#ifndef __KDTREE_H__
#define __KDTREE_H__

#include <vector>
#include <numeric>
#include <algorithm>
#include <exception>
#include <functional>
#include <memory>
#include <limits>
#include <cmath>

namespace kdt
{
	// k-d 树类模板
	template <class PointT>
	class KDTree
	{
	public:
		// 构造函数
		KDTree() : root_(nullptr) {}												  // 初始化树的根节点为 nullptr
		KDTree(const std::vector<PointT> &points) : root_(nullptr) { build(points); } // 用给定的点构建树

		// 析构函数
		~KDTree() { clear(); } // 清除树的所有节点

		// 重建 k-d 树
		void build(const std::vector<PointT> &points)
		{
			clear();											  // 清空当前树
			points_ = points;									  // 保存输入点
			std::vector<int> indices(points.size());			  // 创建索引数组
			std::iota(std::begin(indices), std::end(indices), 0); // 填充索引数组
			// 递归构建 k-d 树并用智能指针管理内存
			root_ = std::unique_ptr<Node>(buildRecursive(indices.data(), (int)points.size(), 0));
		}

		// 清除 k-d 树
		void clear()
		{
			root_.reset();	 // 重置根节点
			points_.clear(); // 清空点的存储
		}

		// 验证 k-d 树的有效性
		bool validate() const
		{
			try
			{
				validateRecursive(root_.get(), 0); // 从根节点开始验证
			}
			catch (const Exception &)
			{
				return false; // 如果验证失败，则返回 false
			}
			return true; // 验证通过
		}

		// 搜索最近邻
		int nnSearch(const PointT &query, double *minDist = nullptr) const
		{
			int guess;												  // 存储最近邻索引
			double _minDist = std::numeric_limits<double>::max();	  // 初始化最小距离为无穷大
			nnSearchRecursive(query, root_.get(), &guess, &_minDist); // 递归搜索最近邻

			if (minDist)
				*minDist = _minDist; // 如果提供了 minDist，返回最小距离

			return guess; // 返回最近邻索引
		}

		// 搜索 k 个最近邻
		std::vector<int> knnSearch(const PointT &query, int k) const
		{
			KnnQueue queue(k);								  // 创建限制大小的优先队列
			knnSearchRecursive(query, root_.get(), queue, k); // 递归搜索 k 个最近邻

			std::vector<int> indices(queue.size()); // 存储 k 个最近邻的索引
			for (size_t i = 0; i < queue.size(); i++)
				indices[i] = queue[i].second; // 提取索引

			return indices; // 返回索引数组
		}

		// 在给定半径内搜索邻居
		std::vector<int> radiusSearch(const PointT &query, double radius) const
		{
			std::vector<int> indices;									// 存储在半径内的邻居索引
			radiusSearchRecursive(query, root_.get(), indices, radius); // 递归搜索
			return indices;												// 返回邻居索引
		}

	private:
		// k-d 树节点结构
		struct Node
		{
			int idx;					   // 当前节点在点数组中的索引
			std::unique_ptr<Node> next[2]; // 左右子树
			int axis;					   // 划分轴

			Node(int index, int a) : idx(index), axis(a) {} // 节点构造函数
		};

		// k-d 树异常类
		class Exception : public std::exception
		{
			using std::exception::exception; // 继承标准异常
		};

		// 限制优先队列类
		template <class T, class Compare = std::less<T>>
		class BoundedPriorityQueue
		{
		public:
			BoundedPriorityQueue(size_t bound) : bound_(bound) { elements_.reserve(bound + 1); } // 初始化队列

			void push(const T &val) // 添加元素
			{
				// 找到插入位置并插入元素
				auto it = std::find_if(std::begin(elements_), std::end(elements_),
									   [&](const T &element)
									   { return Compare()(val, element); });
				elements_.insert(it, val);

				// 超出限制时，调整队列大小
				if (elements_.size() > bound_)
					elements_.resize(bound_);
			}

			const T &back() const { return elements_.back(); }					 // 获取队尾元素
			const T &operator[](size_t index) const { return elements_[index]; } // 获取指定索引元素
			size_t size() const { return elements_.size(); }					 // 返回队列大小

		private:
			size_t bound_;			  // 队列大小限制
			std::vector<T> elements_; // 存储元素
		};

		using KnnQueue = BoundedPriorityQueue<std::pair<double, int>>; // k 个最近邻队列类型

		std::unique_ptr<Node> root_; // 根节点
		std::vector<PointT> points_; // 存储点的向量

		// 递归构建 k-d 树
		// indices 指向点索引的数组 npoints 当前节点所包含的点的数量 depth 当前树的深度，用于选择划分的轴 return 返回构建的节点指针
		Node *buildRecursive(int *indices, int npoints, int depth)
		{
			if (npoints <= 0)
				return nullptr; // 基本情况：没有节点

			const int axis = depth % PointT::DIM; // 确定划分轴
			const int mid = (npoints - 1) / 2;	  // 找到中位数索引

			// 根据当前轴对索引进行排序
			std::nth_element(indices, indices + mid, indices + npoints,
							 [&](int lhs, int rhs)
							 {
								 return points_[lhs][axis] < points_[rhs][axis];
							 });

			// 创建节点
			auto node = std::make_unique<Node>(indices[mid], axis);
			Node *nodePtr = node.get(); // 获取原始指针

			// 递归构建左右子树
			nodePtr->next[0] = std::unique_ptr<Node>(buildRecursive(indices, mid, depth + 1));
			nodePtr->next[1] = std::unique_ptr<Node>(buildRecursive(indices + mid + 1, npoints - mid - 1, depth + 1));

			return node.release(); // 释放节点所有权
		}

		// 递归验证 k-d 树的结构
		void validateRecursive(const Node *node, int depth) const
		{
			if (node == nullptr)
				return; // 基本情况：节点为空

			const int axis = node->axis;			 // 当前轴
			const Node *node0 = node->next[0].get(); // 左子树
			const Node *node1 = node->next[1].get(); // 右子树

			// 检查节点是否满足 k-d 树性质
			if (node0 && node1)
			{
				if (points_[node->idx][axis] < points_[node0->idx][axis])
					throw Exception(); // 违反性质
				if (points_[node->idx][axis] > points_[node1->idx][axis])
					throw Exception(); // 违反性质
			}

			validateRecursive(node0, depth + 1); // 递归验证左子树
			validateRecursive(node1, depth + 1); // 递归验证右子树
		}

		// 计算两个点之间的曼哈顿距离
		static double distance(const PointT &p, const PointT &q)
		{
			double dist = 0;
			for (size_t i = 0; i < PointT::DIM; i++)
				dist += std::fabs(p[i] - q[i]); // 曼哈顿距离计算
			return dist;						// 返回距离
		}

		// 递归搜索最近邻
		// query 查询的点 node 当前节点 guess 存储找到的最近邻的索引 minDist 存储最小距离
		void nnSearchRecursive(const PointT &query, const Node *node, int *guess, double *minDist) const
		{
			if (node == nullptr)
				return; // 基本情况：节点为空，结束递归

			const PointT &train = points_[node->idx];	// 获取当前节点对应的点
			const double dist = distance(query, train); // 计算查询点与当前点的距离

			if ((dist < *minDist) && dist) // 如果当前距离小于已知最小距离(不含自己)
			{
				*minDist = dist;	// 更新最小距离
				*guess = node->idx; // 更新最近邻的索引
			}

			const int axis = node->axis;									 // 当前节点的划分轴
			const int dir = query[axis] < train[axis] ? 0 : 1;				 // 根据查询点与当前点的比较确定搜索方向
			nnSearchRecursive(query, node->next[dir].get(), guess, minDist); // 递归搜索选定方向的子树

			const double diff = std::fabs(query[axis] - train[axis]);			  // 计算查询点与当前点在当前轴上的差异
			if (diff < *minDist)												  // 如果该差异小于已知最小距离
				nnSearchRecursive(query, node->next[!dir].get(), guess, minDist); // 递归搜索另一个方向的子树
		}

		// 递归搜索 k 个最近邻
		//  query 查询的点  node 当前节点 queue 存储最近邻的队列  k需要返回的最近邻的数量
		void knnSearchRecursive(const PointT &query, const Node *node, KnnQueue &queue, int k) const
		{
			if (node == nullptr)
				return; // 基本情况：节点为空，结束递归

			const PointT &train = points_[node->idx];	 // 获取当前节点对应的点
			const double dist = distance(query, train);	 // 计算查询点与当前点的距离
			queue.push(std::make_pair(dist, node->idx)); // 将距离和当前点索引加入优先队列

			const int axis = node->axis;								// 当前节点的划分轴
			const int dir = query[axis] < train[axis] ? 0 : 1;			// 确定搜索方向
			knnSearchRecursive(query, node->next[dir].get(), queue, k); // 递归搜索选定方向的子树

			const double diff = std::fabs(query[axis] - train[axis]); // 计算差异
			// 如果队列中的元素少于 k 或者差异小于队列中最大距离，则搜索另一个方向的子树
			if ((int)queue.size() < k || diff < queue.back().first)
				knnSearchRecursive(query, node->next[!dir].get(), queue, k);
		}

		// 递归搜索指定半径内的邻居
		// query 查询的点 node 当前节点 indices 存储找到的点的索引 radius 搜索半径
		void radiusSearchRecursive(const PointT &query, const Node *node, std::vector<int> &indices, double radius) const
		{
			if (node == nullptr)
				return; // 基本情况：节点为空，结束递归

			const PointT &train = points_[node->idx];	// 获取当前节点对应的点
			const double dist = distance(query, train); // 计算查询点与当前点的距离

			if (dist < radius)				  // 如果距离小于给定半径
				indices.push_back(node->idx); // 将当前点的索引添加到结果中

			const int axis = node->axis;										  // 当前节点的划分轴
			const int dir = query[axis] < train[axis] ? 0 : 1;					  // 确定搜索方向
			radiusSearchRecursive(query, node->next[dir].get(), indices, radius); // 递归搜索选定方向的子树

			const double diff = std::fabs(query[axis] - train[axis]);				   // 计算差异
			if (diff < radius)														   // 如果差异小于给定半径
				radiusSearchRecursive(query, node->next[!dir].get(), indices, radius); // 递归搜索另一个方向的子树
		}
	};
} // namespace kdt

#endif // !__KDTREE_H__
