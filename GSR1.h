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
using namespace std;

double CountRC1(double distance, double net_unit_r, double net_unit_c) // distance为曼哈顿距离
{
    double rc = 0.5 * net_unit_r * net_unit_c * distance * distance;
    return rc;
}

// 声明聚类函数
void GSR1(
    const ProblemData &P_data,
    const ConstraintsData &C_data,
    std::vector<std::vector<Cluster> *> &AllClusters, // 通过引用传递
    std::vector<Point> &points_ff,
    std::vector<Point> &points_buf,
    int layer);

void GSR1(
    const ProblemData &P_data,
    const ConstraintsData &C_data,
    std::vector<std::vector<Cluster> *> &AllClusters, // 通过引用传递
    std::vector<Point> &points_ff,
    std::vector<Point> &points_buf,
    int layer)
{
    double size_ob = P_data.ff_count * P_data.ff_size.x * P_data.ff_size.y;

    // 计算最大距离
    double max_dis = 2 * sqrt((P_data.die_area[2].x * P_data.die_area[2].y - size_ob) / (P_data.ff_count / C_data.max_fanout));
    //cout << "Max Distance: " << max_dis << std::endl;

    // 创建一个新的 vector<Cluster> 实例
    std::vector<Cluster> *Clusters = new std::vector<Cluster>(); // 动态分配内存
    Clusters->emplace_back(Clusters->size());                    // 创建c_index=0的Cluster
    Clusters->back().AddElement(0, 0, 0);
    // 将 Clusters 的地址添加到 AllClusters 中
    AllClusters.push_back(Clusters); // 直接添加指针

    // 添加第一个聚类中心
    Clusters->emplace_back(Clusters->size()); // 初始化第一个 Cluster
    Clusters->back().AddElement(1, P_data.ff_positions[0].x, P_data.ff_positions[0].y);
    points_ff.emplace_back(P_data.ff_positions[0].x, P_data.ff_positions[0].y, 0, P_data.ff_size.x, P_data.ff_size.y); // 添加到点集合

    // 更新聚类中心

    // 遍历每个 ff
    for (int i = 2; i <= P_data.ff_count; i++)
    {
        auto current_ff_pos = P_data.ff_positions[i - 1];
        points_ff.emplace_back(current_ff_pos.x, current_ff_pos.y, 0, P_data.ff_size.x, P_data.ff_size.y); // 添加到点集合
    }
    Overlap overlap(points_ff, points_buf); // 创建重叠检测实例
    for (int i = 2; i <= P_data.ff_count; i++)
    {
        //std::cout << "current ff: " << i << std::endl;
        double min_distance = std::numeric_limits<double>::max(); // 初始化最小距离
        int nearest_cluster = -1;                                 // 最近聚类团索引
        auto current_ff_pos = P_data.ff_positions[i - 1];         // 缓存当前 ff 位置

        // 遍历聚类中心，寻找距离最近的聚类中心
        for (size_t j = 1; j < Clusters->size(); ++j) // 直接使用指针
        {
            auto &cluster = (*Clusters)[j]; // 获取当前 Cluster
            // buf.emplace_back(cluster.GetCX(), cluster.GetCY(), layer, P_data.buf_size.x, P_data.buf_size.y); // 添加到点集合
            double dis = CountDistance(cluster.GetCX(), cluster.GetCY(), current_ff_pos.x, current_ff_pos.y); // 计算当前 ff 至聚类中心的距离
            if (dis < min_distance)
            {
                min_distance = dis;
                nearest_cluster = j;
            }
        }

        // 判断是否满足：最大距离约束和最大扇出约束
        if (!(min_distance < max_dis && (*Clusters)[nearest_cluster].GetFunout() < C_data.max_fanout)) // 不满足
        {
            // 建立该 ff 为新的聚类中心
            Clusters->emplace_back(Cluster(Clusters->size()));                  // 创建新的 Cluster 实例
            Clusters->back().AddElement(i, current_ff_pos.x, current_ff_pos.y); // 添加元素
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
            Clusters->at(0).AddElement(i, current_ff_pos.x, current_ff_pos.y);
        }

        // 计算加入此 ff 后 0 索引聚类团的 rc
        double net_rc = 0;
        for (const auto &elem : Clusters->at(0).ffs) // 使用范围 for 循环
        {
            net_rc += CountRC1(elem.distance, C_data.net_unit_r, C_data.net_unit_c);
        }
        Clusters->at(0).SetRC(net_rc); // 设置新的 rc

        // 判断是否满足最大 rc 约束
        if (Clusters->at(0).GetRC() <= C_data.max_net_rc) // 满足，将此 ff 添加至当前聚类团
        {
            Clusters->at(nearest_cluster).AddElement(i, current_ff_pos.x, current_ff_pos.y);
            Clusters->at(nearest_cluster).SetRC(net_rc); 
                                                                                                     // 将计算好的 rc 直接"复制"到当前聚类团
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
                        net_rc += CountRC1(elem.distance, C_data.net_unit_r, C_data.net_unit_c);
                    }
                    Clusters->at(0).SetRC(net_rc); // 更新聚类的 rc 值
                                                   //
                    if (Clusters->at(0).GetRC() <= C_data.max_net_rc)
                    {
                        Clusters->at(nearest_cluster).SetCX(optimized_position.x); // 更新X坐标
                        Clusters->at(nearest_cluster).SetCY(optimized_position.y); // 更新Y坐标
                        Clusters->at(nearest_cluster).SetRC(net_rc);
                    }
                }
            }
        }
        else // 不满足，建立该 ff 为新的聚类中心
        {
            Clusters->emplace_back(Cluster(Clusters->size())); // 创建新的聚类
            Clusters->back().AddElement(i, current_ff_pos.x, current_ff_pos.y);
        }
    }

    for (int i = 1; i < Clusters->size(); i++)
    {
        points_buf.emplace_back(Clusters->at(i).GetCX(), Clusters->at(i).GetCY(), layer, P_data.buf_size.x, P_data.buf_size.y); // 添加到点集合
    }
}
