#ifndef CLUSTER_H
#define CLUSTER_H

#include <iostream>
#include <vector>
#include <cmath>

// 计算曼哈顿距离
double CountDistance(double x1, double y1, double x2, double y2)
{
    return abs(x1 - x2) + abs(y1 - y2);
}

// 定义 ff 结构体
struct ff
{
    int ff_index;     // ff编号（从1开始）
    double ff_x;      // ff X 坐标
    double ff_y;      // ff Y 坐标
    double distance;  // ff和聚类中心的距离
    int parent_index; // 父节点索引（用于最后计算lantancy指标时，查找父节点）

    // 构造函数：初始化所有成员
    ff(int index, double x, double y, double dist, int parent)
        : ff_index(index), ff_x(x), ff_y(y), distance(dist), parent_index(parent)
    {
    }
};

class Cluster
{
public:
    // 构造函数
    Cluster(int index);

    // 添加元素到 vector
    void AddElement(int ff_index, double x, double y);

    // 打印成员变量
    void Print() const;

    // 获取私有变量的值
    int GetFunout() const { return funout; }
    double GetCX() const { return c_x; }
    double GetCY() const { return c_y; }
    int GetCIndex() const { return c_index; }
    double GetRC() const { return rc; }
    
    // 返回 child_index 的引用
    const std::vector<int>& GetChildIndex() const {
        return child_index;
    }
    
    std::vector<ff> ffs; // 动态数组

    // 修改私有变量的值
    void SetRC(double x);  // 设置聚类团当前的rc
    void SetFanout(int x); // 设置聚类团当前的扇出
    void SetCX(double x);  // 设置聚类团当前的x坐标
    void SetCY(double x);  // 设置聚类团当前的y坐标
    void SetCIndex(int x); // 设置聚类团当前的编号

private:
    int funout;  // 聚类团当前扇出
    double c_x;  // 聚类中心X 坐标
    double c_y;  // 聚类中心Y 坐标
    int c_index; // 聚类团编号（从1开始）
    double rc;   // 聚类团当前rc;(main)中更新
    std::vector<int> child_index; // 动态数组
};

// 构造函数
Cluster::Cluster(int index)
    : c_x(0), c_y(0), c_index(index), funout(0), rc(0) {}

// 添加元素到 vector
void Cluster::AddElement(int ff_index, double ff_x, double ff_y)
{
    ff newElement{ff_index, ff_x, ff_y, 0.0, c_index}; // 初始化 ff 结构体
    ffs.push_back(newElement);

    // 更新聚类中心坐标
    funout = ffs.size(); // 更新 funout
    c_x = (c_x * (funout - 1) + ff_x) / funout;
    c_y = (c_y * (funout - 1) + ff_y) / funout;
    child_index.push_back(ff_index);

    // 更新每个 ff 到新的聚类中心的距离
    for (auto& elem : ffs)
    {
        elem.distance = CountDistance(elem.ff_x, elem.ff_y, c_x, c_y);
    }
}

// 设置聚类中心的rc
void Cluster::SetRC(double x)
{
    rc = x;
}

// 设置聚类中心的扇出
void Cluster::SetFanout(int x)
{
    funout = x;
}

// 设置聚类中心的x坐标
void Cluster::SetCX(double x)
{
    c_x = x;
}

// 设置聚类中心的y坐标
void Cluster::SetCY(double x)
{
    c_y = x;
    
    // 更新每个 ff 到新的聚类中心的距离
    for (auto& elem : ffs)
    {
        elem.distance = CountDistance(elem.ff_x, elem.ff_y, c_x, c_y);
    }
}

// 设置聚类中心的编号
void Cluster::SetCIndex(int x)
{
    c_index = x;
}

// 打印成员变量
void Cluster::Print() const
{
    std::cout << "Cluster Coordinates: (" << c_x << ", " << c_y << ")" << std::endl;
    std::cout << "Funout: " << funout << std::endl;
    std::cout << "net rc:" << rc << std::endl;
    std::cout << "Elements in ffs:" << std::endl;

    for (const auto &elem : ffs)
    {
        std::cout << "Index: " << elem.ff_index
                  << ", Coordinates: (" << elem.ff_x << ", " << elem.ff_y
                  << "), Distance: " << elem.distance << std::endl;
    }
}

#endif // CLUSTER_H
