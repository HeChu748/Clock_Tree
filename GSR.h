#pragma once

#include <vector>
#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <stdio.h>
#include <algorithm>
#include "constrain.h"
#include "problem.h"
#include "Cluster.h"
#include "kdtree.h"
#include "OverlapBob.h"
#include "./BufferLocationSA.h"
using namespace std;

std::vector<std::vector<int>> cluster_layer_mapping;
double CountRC(double distance, double net_unit_r, double net_unit_c) // distance为曼哈顿距离
{
    double rc = 0.5 * net_unit_r * net_unit_c * distance * distance;
    return rc;
}

// 声明聚类函数
void GSR(
    const ProblemData &P_data,
    const ConstraintsData &C_data,
    std::vector<std::vector<Cluster> *> &AllClusters, // 通过引用传递
    std::vector<Point> &points_ff,
    std::vector<Point> &points_buf,
    int layer);

void GSR(
    const ProblemData &P_data,
    const ConstraintsData &C_data,
    std::vector<std::vector<Cluster> *> &AllClusters, // 通过引用传递
    std::vector<Point> &points_ff,
    std::vector<Point> &points_buf,
    int layer)
{
    // 创建一个新的 vector<Cluster> 实例
    std::vector<Cluster> *Clusters = new std::vector<Cluster>(); // 动态分配内存
    Clusters->emplace_back(0);                                   // 创建c_index=0的Cluster
    Clusters->back().AddElement(0, 0, 0);
    // 将 Clusters 的地址添加到 AllClusters 中
    AllClusters.push_back(Clusters);                             // 直接添加指针
    std::vector<Cluster> *LastClusters = AllClusters[layer - 2]; // 上一级聚类

    // 添加第一个聚类中心
    Clusters->emplace_back(Clusters->size()); // 初始化第一个 Cluster
    Clusters->back().AddElement(1, LastClusters->at(1).GetCX(), LastClusters->at(1).GetCY());

    // 更新聚类中心

    double size_ob = (LastClusters->size() - 1) * P_data.buf_size.x * P_data.buf_size.y; //
    // 计算最大距离
    double max_dis = 2 * sqrt((P_data.die_area[2].x * P_data.die_area[2].y - size_ob) / ((LastClusters->size() - 1) / C_data.max_fanout));
    //cout << "Max Distance: " << max_dis << std::endl;

    Overlap overlap(points_ff, points_buf);

    // 遍历每个 buf
    for (int i = 2; i < LastClusters->size(); i++)
    {
        double min_distance = std::numeric_limits<double>::max(); // 初始化最小距离
        int nearest_cluster = -1;                                 // 最近聚类团索引
        double current_buf_x = LastClusters->at(i).GetCX();       // 缓存当前 buf 位置
        double current_buf_y = LastClusters->at(i).GetCY();

        // 遍历聚类中心，寻找距离最近的聚类中心
        for (size_t j = 1; j < Clusters->size(); ++j) // 直接使用指针
        {
            auto &cluster = (*Clusters)[j];                                                             // 获取当前 Cluster
            double dis = CountDistance(cluster.GetCX(), cluster.GetCY(), current_buf_x, current_buf_y); // 计算当前 duf 至聚类中心的距离
            if (dis < min_distance)
            {
                min_distance = dis;
                nearest_cluster = j;
            }
        }

        // 判断是否满足：最大距离约束和最大扇出约束
        if (!(min_distance < max_dis && (*Clusters)[nearest_cluster].GetFunout() < C_data.max_fanout)) // 不满足
        {
            // 建立该 buf 为新的聚类中心
            Clusters->emplace_back(Cluster(Clusters->size()));            // 创建新的 Cluster 实例
            Clusters->back().AddElement(i, current_buf_x, current_buf_y); // 添加元素
            continue;
        }

        else // 满足，继续判断是否满足最大 rc 约束
        {
            // 将当前 Cluster 复制到 0 索引 Cluster 中
            Clusters->at(0).SetFanout((*Clusters)[nearest_cluster].GetFunout());
            Clusters->at(0).SetCX((*Clusters)[nearest_cluster].GetCX());
            Clusters->at(0).SetCY((*Clusters)[nearest_cluster].GetCY());

            Clusters->at(0).ffs.clear();                            // 清空 ffs 数组
            Clusters->at(0).ffs = (*Clusters)[nearest_cluster].ffs; // 复制 ffs 数组

            // 将当前 ff 添加至 0 索引聚类中心
            Clusters->at(0).AddElement(i, current_buf_x, current_buf_y);

            // 计算加入此 ff 后 0 索引聚类团的 rc
            double net_rc = 0;
            for (const auto &elem : Clusters->at(0).ffs) // 使用范围 for 循环
            {
                net_rc += CountRC(elem.distance, C_data.net_unit_r, C_data.net_unit_c);
            }
            Clusters->at(0).SetRC(net_rc); // 设置新的 rc

            // 判断是否满足最大 rc 约束
            if (Clusters->at(0).GetRC() <= C_data.max_net_rc) // 满足，将此 ff 添加至当前聚类团
            {
                Clusters->at(nearest_cluster).AddElement(i, current_buf_x, current_buf_y);
                Clusters->at(nearest_cluster).SetRC(net_rc); // 将计算好的 rc 直接"复制"到当前聚类团

                // 优化中心位置
                ProblemData::Coordinate best_location = optimizeBufferLocation_SA(Clusters->at(nearest_cluster), *LastClusters);
                // 更新聚类中心位置
                if (best_location.x - P_data.buf_size.x >= P_data.die_area[0].x && best_location.x + P_data.buf_size.x <= P_data.die_area[2].x && best_location.y - P_data.buf_size.y >= P_data.die_area[0].y && best_location.y + P_data.buf_size.y <= P_data.die_area[2].y)
                {
                    // 更新聚类中心位置
                    Clusters->at(0).SetCX(best_location.x); // 更新X坐标
                    Clusters->at(0).SetCY(best_location.y); // 更新Y坐标

                    // 更新 rc(可优化)
                    double net_rc = 0;
                    for (const auto &elem : Clusters->at(0).ffs) // 使用范围 for 循环
                    {
                        net_rc += CountRC(elem.distance, C_data.net_unit_r, C_data.net_unit_c);
                    }
                    Clusters->at(0).SetRC(net_rc); // 更新聚类的 rc 值
                                                   //
                    if (Clusters->at(0).GetRC() <= C_data.max_net_rc)   //判断优化后是否满足rc约束
                    {
                        Clusters->at(nearest_cluster).SetCX(best_location.x); // 返回优化的坐标
                        Clusters->at(nearest_cluster).SetCY(best_location.y);
                        Clusters->at(nearest_cluster).SetRC(net_rc);
                    }
                } // 聚类中心优化结束

                // 移动当前聚类，处理重叠
                Point buf(Clusters->at(nearest_cluster).GetCX(), Clusters->at(nearest_cluster).GetCY(), layer, P_data.buf_size.x, P_data.buf_size.y); // 添加到点集合
                overlap.setCurrentBuf(buf);
                ProblemData::Coordinate optimized_position = overlap.OverlapBob(buf); // 进行重叠处理并获取优化位置

                if ((Clusters->at(nearest_cluster).GetCX() != optimized_position.x) || (Clusters->at(nearest_cluster).GetCY() != optimized_position.y))
                {

                    if (optimized_position.x - P_data.buf_size.x >= P_data.die_area[0].x && optimized_position.x + P_data.buf_size.x <= P_data.die_area[2].x && optimized_position.y - P_data.buf_size.y >= P_data.die_area[0].y && optimized_position.y + P_data.buf_size.y <= P_data.die_area[2].y)
                    {
                        // 更新聚类中心位置
                        Clusters->at(0).SetCX(optimized_position.x); // 更新X坐标
                        Clusters->at(0).SetCY(optimized_position.y); // 更新Y坐标

                        // 更新 rc(可优化)
                        double net_rc = 0;
                        for (const auto &elem : Clusters->at(0).ffs) // 使用范围 for 循环
                        {
                            net_rc += CountRC(elem.distance, C_data.net_unit_r, C_data.net_unit_c);
                        }
                        Clusters->at(0).SetRC(net_rc); // 更新聚类的 rc 值

                        if (Clusters->at(0).GetRC() <= C_data.max_net_rc)   //判断优化后是否满足rc约束
                        {
                            Clusters->at(nearest_cluster).SetCX(optimized_position.x); // 返回优化后坐标
                            Clusters->at(nearest_cluster).SetCY(optimized_position.y);
                            Clusters->at(nearest_cluster).SetRC(net_rc);
                        }
                    }
                } // 重叠优化结束
            }
            else // 不满足，建立该 ff 为新的聚类中心
            {
                Clusters->emplace_back(Cluster(Clusters->size())); // 创建新的聚类
                Clusters->back().AddElement(i, current_buf_x, current_buf_y);
            }
        }
    }


    for (int i = 1; i < Clusters->size(); i++)
    {
        points_buf.emplace_back(Clusters->at(i).GetCX(), Clusters->at(i).GetCY(), layer, P_data.buf_size.x, P_data.buf_size.y); // 添加到点集合
    }
}
