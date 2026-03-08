#include <iostream>
#include <cmath>
#include <cstring>
#include <vector>
#include <iomanip>
#include <cstdlib>  //包含命令行参数处理库
#include "constrain.h"
#include "problem.h"
#include <stdio.h>
#include <chrono> // 包含时间库
using namespace std;
#include "kdtree.h"
#include "Cluster.h"
#include "GSR.h"
#include "GSR1.h"
#include "./OverlapBob.h"

void printClusterInfo(const std::vector<Cluster> *clusters, int current_buf_num,double buf_width, double buf_length, std::ostream& out)
{
    for (size_t j = 1; j < clusters->size(); ++j)
    {
        double cx = clusters->at(j).GetCX();
        double cy = clusters->at(j).GetCY();

        out << "- BUF" << j - 1 + current_buf_num << " BUF ( " << cx - buf_width / 2 << " " << cy - buf_length / 2 << " ) ;" << std::endl;
    }
}

void printNetInfo(const std::vector<Cluster> *clusters, int lastlayer_buf_num, int current_buf_num, bool isFirstLayer, std::ostream& out)
{
    for (size_t j = 1; j < clusters->size(); ++j)
    {
        out << "- net_" << j - 1 + current_buf_num << " ( BUF" << j - 1 + current_buf_num << " ) ( ";

        const std::vector<int> &children = clusters->at(j).GetChildIndex();
        for (size_t k = 0; k < children.size(); ++k)
        {
            out << (isFirstLayer ? "FF" : "BUF") << (isFirstLayer ?children[k] + lastlayer_buf_num : children[k] + lastlayer_buf_num-1);
            if (k < children.size() - 1)
            {
                out << " ";
            }
        }
        out << " ) ;" << std::endl;
    }
}

int main(int argc, char* argv[]) {
    // 检查命令行参数
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <folder_path> constraints_file=<folder_path>/constraints.txt problem_file=<folder_path>/problem.def" << std::endl;
        return 1;
    }
    // 检查是否需要帮助信息
    if (std::string(argv[1]) == "-h") {
        std::cout << "This program requires one argument for the folder path:" << std::endl;
        std::cout << "Usage: " << argv[0] << " <folder_path>" << std::endl;
        std::cout << "The folder should contain the following files:" << std::endl;
        std::cout << "1. constraints.txt - Path to the constraints file." << std::endl;
        std::cout << "2. problem.def - Path to the problem definition file." << std::endl;
        std::cout << "The output will be saved as solution.def in the same folder." << std::endl;
        return 0;
    }
    // 获取文件夹路径
    std::string folder_path = argv[1];
    
    // 构建文件路径
    std::string constraints_file = folder_path + "/constraints.txt";
    std::string problem_file = folder_path + "/problem.def";

    // 打印文件路径，检查是否正确
    //std::cout << "Constraints file: " << constraints_file << std::endl;
    //std::cout << "Problem file: " << problem_file << std::endl;

    auto start = chrono::high_resolution_clock::now(); // 记录开始时间

    // 创建文件输出流对象并打开文件
    std::ofstream output_file(folder_path + "/solution.def");
    if (!output_file.is_open()) {
        std::cerr << "Error opening file." << std::endl;
        return 1;
    }

    //设置输出格式
    output_file << fixed << setprecision(0);
    ConstraintsData C_data;
    ProblemData P_data;

    // cout << "Constrain Data:" << "\n";
    readConstraints(constraints_file, C_data);
    //C_data.display();

    // cout << "Problem Data:" << "\n";
    readProblemDef(problem_file, P_data);
    //P_data.display();

    output_file << "UNITS DISTANCE MICRONS " <<  P_data.units_distance << " ;" << endl;
    output_file << "DIEAREA " ;
    for (const auto& coord : P_data.die_area){
      output_file << "( " << coord.x << " "<< coord.y << " ) " ;
    }
    output_file << ";" << endl;
    output_file << "FF " << "( " << P_data.ff_size.x << " " << P_data.ff_size.y << " ) ;" << endl;
    output_file << "BUF " << "( " << P_data.buf_size.x << " " << P_data.buf_size.y << " ) ;" << endl;
    output_file << "CLK " << "( " << P_data.clk_position.x << " " << P_data.clk_position.y << " ) ;" << endl;

    int layer = 1;

    // 初始化 AllClusters 和 points
    std::vector<std::vector<Cluster> *> AllClusters; // 定义并初始化 AllClusters

    //std::vector<Point> points;                 // kd树数组
    std::vector<Point> points_ff;                       // kd树数组
    std::vector<Point> points_buf;
    //std::vector<std::vector<int>> cluster_layer_mapping; // 声明 cluster_layer_mapping

    // 第一层聚类
    GSR1(P_data, C_data, AllClusters, points_ff,points_buf, layer);

    std::vector<Cluster> *clusters = AllClusters[AllClusters.size() - 1]; // 获取指向 vector<Cluster> 的指针
   /* while (clusters->size() > C_data.max_fanout)
    {
        layer++;
        GSR(P_data, C_data, AllClusters, points_ff, points_buf, layer );
        clusters = AllClusters[AllClusters.size() - 1]; // 更新指向 vector<Cluster> 的指针
    }*/
    double rc_sum = std::numeric_limits<double>::infinity(); //设为无穷大
    size_t previous_cluster_size = clusters->size();
    while (rc_sum > C_data.max_net_rc )
    {
        layer++;
        GSR(P_data, C_data, AllClusters, points_ff, points_buf, layer );
        clusters = AllClusters[AllClusters.size() - 1]; // 更新指向 vector<Cluster> 的指针
        rc_sum = 0;
        for (const auto& cluster : *clusters){
          double distance = CountDistance(cluster.GetCX(), cluster.GetCY(), P_data.clk_position.x, P_data.clk_position.y);
          double RC = CountRC(distance, C_data.net_unit_r, C_data.net_unit_c);
          rc_sum += RC;   //更新rc_sum
        }
        //cout << "rc_sum:"<< rc_sum <<endl;
        
        //cout << "clusters size" << clusters->size() << endl ;


        //检查是否聚类到极限
        if (clusters->size() == previous_cluster_size){
          break;
        }
        previous_cluster_size = clusters->size();
    }

    // 顶层聚类
    // 确定Topbuf位置
    std::vector<std::vector<ff> *> TopClusters; // 定义并初始化 AllClusters
    ProblemData::Coordinate Topbuf = {0, 0};
    double sum_x = 0;
    double sum_y = 0;
    std::vector<ff> *ffs; // 动态分配内存

    for (int i = 1; i < clusters->size(); i++)
    {
        // 创建一个新的 vector<ff> 实例
        std::vector<ff> *ffs = new std::vector<ff>(); // 动态分配内存
        ffs->emplace_back(ff{0, clusters->at(i).GetCX(), clusters->at(i).GetCY(), 0, 0});
        // 将 ffs 的地址添加到 TopClusters 中
        TopClusters.push_back(ffs); // 直接添加指针
        sum_x += clusters->at(i).GetCX();
        sum_y += clusters->at(i).GetCY();
    }
    Topbuf.x = sum_x / (clusters->size() - 1);
    Topbuf.y = sum_y / (clusters->size() - 1);

    // 向Topbuf"移动"
    double max_distance = static_cast<int>(std::trunc(sqrt(2 * C_data.max_net_rc / (C_data.net_unit_c * C_data.net_unit_r))));
    for (size_t i = 0; i < TopClusters.size(); ++i)
    {
        ffs = TopClusters[i];
        double current_dis = CountDistance(ffs->at(0).ff_x, ffs->at(0).ff_y, Topbuf.x, Topbuf.y);
        double current_rc = CountRC(current_dis, C_data.net_unit_r, C_data.net_unit_c);
        int x = 0; //"移动"到Topbuf还需插入的buf数+1
        x = static_cast<int>(std::trunc(current_dis / max_distance));

        double dis_x = Topbuf.x - ffs->at(0).ff_x;
        double dis_y = Topbuf.y - ffs->at(0).ff_y;
        double flag_x = 1;
        double flag_y = 1;
        if (dis_x < 0)
            flag_x = -1;
        if (dis_y < 0)
            flag_y = -1;

        ProblemData::Coordinate Addbuf = {ffs->at(0).ff_x, ffs->at(0).ff_y};
        for (int j = 0; j < x; j++)
        {
            Addbuf.x += flag_x * (max_distance * (abs(dis_x) / (abs(dis_x) + abs(dis_y))));
            Addbuf.y += flag_y * (max_distance * (abs(dis_y) / (abs(dis_x) + abs(dis_y))));

            ffs->emplace_back(ff{j + 1, Addbuf.x, Addbuf.y, 0, 0});
        }

        if (current_rc - x * C_data.max_net_rc > (C_data.max_net_rc / TopClusters.size()))
        {
            x++;
            Addbuf.x = Topbuf.x - flag_x * ((max_distance / TopClusters.size()) * (abs(dis_x) / (abs(dis_x) + abs(dis_y))));
            Addbuf.y = Topbuf.y - flag_y * ((max_distance / TopClusters.size()) * (abs(dis_y) / (abs(dis_x) + abs(dis_y))));

            ffs->emplace_back(ff{x, Addbuf.x, Addbuf.y, 0, 0});
        }
    }

    // Topbuf聚到clk
    //  创建一个新的 vector<ff> 实例
    ffs = new std::vector<ff>(); // 动态分配内存
    ffs->emplace_back(ff{0, Topbuf.x, Topbuf.y, 0, 0});
    TopClusters.push_back(ffs); // 直接添加指针

    double current_dis = CountDistance(P_data.clk_position.x, P_data.clk_position.y, Topbuf.x, Topbuf.y);
    double current_rc = CountRC(current_dis, C_data.net_unit_r, C_data.net_unit_c);
    int x = 0; //"移动"到clk还需插入的buf数
    x = static_cast<int>(std::trunc(current_dis / max_distance));

    double dis_x = P_data.clk_position.x - Topbuf.x;
    double dis_y = P_data.clk_position.y - Topbuf.y;
    double flag_x = 1;
    double flag_y = 1;
    if (dis_x < 0)
        flag_x = -1;
    if (dis_y < 0)
        flag_y = -1;

    ProblemData::Coordinate Addbuf = {Topbuf.x, Topbuf.y};
    for (int j = 0; j < x; j++)
    {
        Addbuf.x += flag_x * (max_distance * (abs(dis_x) / (abs(dis_x) + abs(dis_y))));
        Addbuf.y += flag_y * (max_distance * (abs(dis_y) / (abs(dis_x) + abs(dis_y))));

        ffs->emplace_back(ff{j + 1, Addbuf.x, Addbuf.y, 0, 0});
    }
   

    // 格式输出
    int lastlayer_buf_num = 0; // 第一层至上一层buf总数
    int current_buf_num = 0;
    
    for (size_t i = 0; i < AllClusters.size(); ++i) {
        clusters = AllClusters[i];
        current_buf_num += clusters->size() - 1;
    }
   for (size_t i = 0; i < TopClusters.size(); ++i)
     {
        ffs = TopClusters[i];
        current_buf_num += ffs->size() - 1;
        if (i == TopClusters.size() - 1)
         {
            if (TopClusters.size() == 2)
                continue;
            else
                current_buf_num++;
         }
     }
    //输出ff和buffer总数
    output_file << "COMPONENTS " << current_buf_num + P_data.ff_count << " ;" << endl;
    //输出ff坐标
    int index = 1;
    for (const auto& pos : P_data.ff_positions){
      output_file << "- FF" << index << " FF " << "( " << pos.x - P_data.ff_size.x/2 << " " << pos.y - P_data.ff_size.y/2 << " ) ;" << endl;
      index++;
    }

    lastlayer_buf_num = 0;
    current_buf_num = 0;
    for (size_t i = 0; i < AllClusters.size(); ++i)
    {
        clusters = AllClusters[i];
        printClusterInfo(clusters, current_buf_num, P_data.buf_size.x, P_data.buf_size.y, output_file);
        current_buf_num += clusters->size() - 1;
    }
    //cout << "current_buf_num: " << current_buf_num << endl ;
     //  输出顶层聚类Cluster
    for (size_t i = 0; i < TopClusters.size(); ++i)
    {
        ffs = TopClusters[i];
        if (i == TopClusters.size() - 1)
        {
            if (TopClusters[0]->at(0).ff_x != Topbuf.x || TopClusters[0]->at(0).ff_y != Topbuf.y)
            {
                output_file << "- BUF" << current_buf_num << " BUF ( " << Topbuf.x << " " << Topbuf.y << " ) ;" << std::endl;
                current_buf_num++;
            }
        }
        if (ffs->size() == 1)
        {
            continue;
        }
        else
        {
            for (size_t j = 1; j < ffs->size(); ++j)
            {
                double cx = ffs->at(j).ff_x;
                double cy = ffs->at(j).ff_y;
                output_file << "- BUF" << j - 1 + current_buf_num << " BUF ( " << cx << " " << cy << " ) ;" << std::endl;
            }
        }
        current_buf_num += ffs->size() - 1;
    }
    output_file << "END COMPONENTS ;" << std::endl;
    output_file << "NETS " << current_buf_num + 1 << " ;" << std::endl;

    lastlayer_buf_num = 0;
    current_buf_num = 0;
    
    for (size_t i = 0; i < AllClusters.size(); ++i)
    {

        if (!i)
        {
            clusters = AllClusters[i];
            printNetInfo(clusters, 0, current_buf_num, i == 0, output_file);
            current_buf_num += clusters->size() - 1; // 更新 current_buf_num
        }
        else
        {
            double tempt = clusters->size() - 1;
            clusters = AllClusters[i];
            printNetInfo(clusters, lastlayer_buf_num, current_buf_num, i == 0, output_file);
            current_buf_num += clusters->size() - 1; // 更新 current_buf_num
            lastlayer_buf_num += tempt;
        }
    }

    // 输出顶层聚类net
    for (size_t i = 0; i < TopClusters.size() - 1; ++i)
    {
        ffs = TopClusters[i];
        for (size_t j = 1; j < ffs->size(); ++j)
        {
            if (j == 1)
            {
                output_file << "- net_" << j - 1 + current_buf_num << " ( BUF" << j - 1 + current_buf_num << " ) ( ";
                output_file << "BUF" << lastlayer_buf_num + i << " ) ;" << std::endl;
                continue;
            }
            output_file << "- net_" << j - 1 + current_buf_num << " ( BUF" << j - 1 + current_buf_num << " ) ( ";
            output_file << "BUF" << j - 1 + current_buf_num - 1 << " ) ;" << std::endl;
        }
        current_buf_num += ffs->size() - 1;
    }
     //cout << "TOPClusters_size" << TopClusters.size() << endl;
    // 输出Topbuf net
    if (TopClusters.size() > 2)
    {
        output_file << "- net_" << current_buf_num << " ( BUF" << current_buf_num << " ) ( ";
        lastlayer_buf_num +=  TopClusters.size() - 1;
        //cout << "lastlayer_buf_num : " << lastlayer_buf_num << endl;

        int last_Top_num = 0;
        for (size_t i = 0; i < TopClusters.size() - 1; ++i)
        {   //cout <<"last_Top_num: " << last_Top_num << endl;
            if(TopClusters[i]->size()==1)
            {
                output_file << "BUF" << lastlayer_buf_num -TopClusters.size() + i +1 << " ";
            }
            else
            {
                last_Top_num += TopClusters[i]->size() - 1 ;
                output_file << "BUF" << lastlayer_buf_num + last_Top_num - 1<< " ";
            }
        }
        
        output_file << ") ;" << std::endl;
        current_buf_num++;
    }

    ffs = TopClusters[TopClusters.size() - 1];
    for (size_t j = 1; j < ffs->size(); ++j)
    {
        output_file << "- net_" << j - 1 + current_buf_num << " ( BUF" << j - 1 + current_buf_num << " ) ( ";
        output_file << "BUF" << j - 1 + current_buf_num - 1 << " ) ;" << std::endl;
    }
    current_buf_num += ffs->size() - 1;

    output_file << "- net_" << current_buf_num << " ( CLK ) ( ";
    output_file << "BUF" << current_buf_num - 1 << " ) ;" << std::endl;
    output_file << "END NETS ;" << std::endl;

    auto end = chrono::high_resolution_clock::now(); // 记录结束时间
    chrono::duration<double> duration = end - start; // 计算用时
    cout << "Time taken: " << duration.count() << " seconds" << std::endl;
    
    return 0;
}
