
#include "Clustering.h"

#include <lemon/list_graph.h>
#include <lemon/maps.h>
#include <lemon/network_simplex.h>
#include <sys/timeb.h>

#include <algorithm>
#include <array>
#include <cfloat>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <stack>
#include <string>
#include <vector>

#include "../utility/Log.hh"
#include "../utility/Time.hh"
#include "../utility/Usage.hh"


namespace cts::CKMeans {

  using lemon::INVALID;
  using lemon::ListDigraph;
  using lemon::NetworkSimplex;


  Clustering::Clustering(const std::vector<std::pair<float, float>>& sinks,
    std::vector<int>& sinks_idx,
    double r, double c, double max_rc)
    :net_unit_r_(r), net_unit_c_(c), max_net_rc_(max_rc)
  {
    sinks_.reserve(sinks.size());
    for (size_t i = 0; i < sinks.size(); ++i) {
      sinks_.emplace_back(sinks[i].first, sinks[i].second, sinks_idx[i]);
    }
    srand(56);
  }

  Clustering::~Clustering() = default;

  /*** Capacitated K means **************************************************/
  void Clustering::iterKmeans(const unsigned iter,    // 迭代次数
    const unsigned n,       // 簇的数量
    const unsigned cap,     // 每个簇的最大容量
    const unsigned max,     // 控制迭代最大次数的阈值
    const unsigned power,   // 计算距离时的权重因子
    std::vector<std::pair<float, float>>& means)
  {
    std::vector<int> solution(sinks_.size());   // 保存最佳聚类配置的索引，即每个 sink 最终所属的聚类
    float max_silh = -1;    // 跟踪迭代过程中观察到的最大 silhouette 分数
    auto tmp_means = means;
    for (unsigned i = 0; i < iter; ++i) {

      // 调用 Kmeans 函数来执行一次聚类过程，并更新聚类中心点，返回计算得到的 silhouette 分数
      const float silh = Kmeans(n, cap, max, power, tmp_means);

      // 检查这次迭代的 silhouette 分数是否比之前记录的最大分数 max_silh 更高
      if (silh > max_silh) {
        max_silh = silh;
        for (size_t j = 0; j < sinks_.size(); ++j) {
          solution[j] = sinks_[j].cluster_idx;
        }
        // 如果当前 silhouette 分数是新的最大值，则更新聚类中心
        means = tmp_means;
      }
    }

    // 根据 solution 数组中记录的聚类索引来更新每个 sink 的 cluster_idx
    for (size_t i = 0; i < sinks_.size(); ++i) {
      sinks_[i].cluster_idx = solution[i];
    }

  }


  float Clustering::Kmeans(const unsigned n,
    const unsigned cap,
    const unsigned max,
    const unsigned power,
    std::vector<std::pair<float, float>>& means)
  {
    // 为所有sink初始化聚类索引为-1，表示尚未分配到任何聚类
    for (auto& sink : sinks_) {
      sink.cluster_idx = -1;
    }

    // 初始化聚类向量，每个聚类初始时都没有任何sink
    std::vector<std::vector<Sink*>> clusters;

    // 控制聚类迭代的条件变量
    bool stop = false;
    unsigned iter = 1;  // 迭代计数器

    // 聚类优化过程
    while (!stop) {

      // 计算当前所有sink与聚类中心的最小成本流，进行sink到聚类的匹配
      minCostFlow(means, cap, max_net_rc_, power);

      // 清空上一轮的聚类结果
      clusters.clear();
      clusters.resize(n);

      // 先将有簇号的sink放到簇中
      for (Sink& sink : sinks_) {
        // 初始化聚类位置
        int position = 0;
        if (sink.cluster_idx >= 0 && sink.cluster_idx < n) {
          position = sink.cluster_idx;
          // 将sink加入到对应的聚类中
          clusters[position].push_back(&sink);
        }
      }
      // 再处理簇号不满足要求的sink
      for (Sink& sink : sinks_) {
        int position = 0;
        if (sink.cluster_idx < 0 || sink.cluster_idx >= n) {
          // 错误分配检查与处理
          float minimumDist;
          int minimumDistClusterIndex = -1;
          // 遍历所有聚类中心，找到未满员的聚类中心
          for (size_t j = 0; j < means.size(); ++j) {
            if (clusters[j].size() < cap) {
              minimumDist = calcDist(
                std::make_pair(means[j].first, means[j].second), &sink);
              minimumDistClusterIndex = j;
              break;
            }
          }
          // 如果没有找到未满员的聚类中心，则默认分配到第一个聚类
          if (minimumDistClusterIndex == -1) {
            minimumDistClusterIndex = 0;
          }
          else {
            // 找到距离当前sink最近的未满员聚类中心
            for (size_t j = 0; j < means.size(); ++j) {
              if (clusters[j].size() < cap) {
                const float currentDist = calcDist(
                  std::make_pair(means[j].first, means[j].second), &sink);
                if (currentDist < minimumDist) {
                  minimumDist = currentDist;
                  minimumDistClusterIndex = j;
                }
              }
            }
          }
          // 更新聚类位置
          position = minimumDistClusterIndex;
          // 将sink加入到对应的聚类中
          clusters[position].push_back(&sink);
        }
      }

      // 计算新聚类中心并判断是否继续迭代
      float delta = 0;
      for (unsigned i = 0; i < n; ++i) {
        float sum_x = 0, sum_y = 0;
        for (const auto& cluster : clusters[i]) {
          sum_x += cluster->x;
          sum_y += cluster->y;
        }
        const float pre_x = means[i].first;
        const float pre_y = means[i].second;

        if (!clusters[i].empty()) {
          means[i] = std::make_pair(sum_x / clusters[i].size(),
            sum_y / clusters[i].size());
          delta += std::abs(pre_x - means[i].first)
            + std::abs(pre_y - means[i].second);
        }
      }

      // 更新全局聚类结果
      clusters_ = clusters;

      // 判断迭代是否应停止
      if (iter > max || delta < 0.5) {
        stop = true;
      }
      ++iter;
    }
    // 计算并返回silhouette分数
    return calcSilh(means);
  }

  float Clustering::calcSilh(
    const std::vector<std::pair<float, float>>& means) const
  {
    float sum_silh = 0;
    for (const Sink& sink : sinks_) {
      float in_d = 0, out_d = FLT_MAX;
      for (size_t j = 0; j < means.size(); ++j) {
        const float x = means[j].first;
        const float y = means[j].second;
        if (sink.cluster_idx == j) {
          // within the cluster
          in_d = calcDist({ x, y }, &sink);
        }
        else {
          // outside of the cluster
          const float d = calcDist({ x, y }, &sink);
          out_d = std::min(d, out_d);
        }
      }
      const float temp = std::max(out_d, in_d);
      if (temp == 0) {
        if (out_d == 0) {
          sum_silh += -1;
        }
        if (in_d == 0) {
          sum_silh += 1;
        }
      }
      else {
        sum_silh += (out_d - in_d) / temp;
      }
    }
    return sum_silh / sinks_.size();
  }


  void Clustering::minCostFlow(const std::vector<std::pair<float, float>>& means,
    const unsigned cap,
    const double max_net_rc,
    const unsigned power)
  {
    // 创建有向图
    ListDigraph graph;

    // 添加源点和汇点
    ListDigraph::Node src = graph.addNode();
    ListDigraph::Node target = graph.addNode();

    // 存储所有sink节点和聚类中心节点
    std::vector<ListDigraph::Node> sink_nodes, cluster_nodes;

    // 为每个sink和每个聚类中心添加节点
    for (size_t i = 0; i < sinks_.size(); ++i) {
      sink_nodes.push_back(graph.addNode());
    }
    for (size_t i = 0; i < means.size(); ++i) {
      cluster_nodes.push_back(graph.addNode());
    }

    // 从源点到每个sink添加边，代表流动可能性
    ListDigraph::ArcMap<int> edge_capacity(graph);  // 新增：定义边的容量映射
    for (auto& sink : sink_nodes) {
      ListDigraph::Arc edge = graph.addArc(src, sink);
      edge_capacity[edge] = 1;  // 每个sink只能流出一个单位的流量
    }

    // 存储每条边的成本，用于最小成本流计算
    double half_rc = 0.5 * net_unit_r_ * net_unit_c_;
    ListDigraph::ArcMap<double> edge_cost(graph);
    for (size_t i = 0; i < sink_nodes.size(); ++i) {
      for (size_t j = 0; j < means.size(); ++j) {
        double net_rc = calcDist(means[j], &sinks_[i]);
        net_rc = std::pow(net_rc, power);
        net_rc = half_rc * net_rc;
        if (net_rc <= max_net_rc / 2) {
          ListDigraph::Arc edge = graph.addArc(sink_nodes[i], cluster_nodes[j]);
          edge_cost[edge] = net_rc;
          edge_capacity[edge] = 1;
        }
      }
    }

    // 从每个聚类中心到汇点添加边，代表聚类可以接受的总流量
    for (auto& cluster : cluster_nodes) {
      ListDigraph::Arc edge = graph.addArc(cluster, target);
      edge_capacity[edge] = cap; // 确保所有簇的容量都为 cap
    }

    // 配置并运行NetworkSimplex算法
    NetworkSimplex<ListDigraph, int, int> flow(graph);
    flow.costMap(edge_cost);
    flow.upperMap(edge_capacity);
    flow.stSupply(src, target, sinks_.size()); // 设置源点和汇点的总供给量，确保与sink总数相匹配

    if (flow.run() == NetworkSimplex<ListDigraph>::OPTIMAL) {
      ListDigraph::ArcMap<int> solution(graph);  // 新增：存储流量解析结果
      flow.flowMap(solution);

      // 根据流量解析结果更新sink的聚类索引
      for (ListDigraph::OutArcIt a(graph, src); a != INVALID; ++a) {
        if (solution[a] > 0) {
          ListDigraph::Node u = graph.target(a);
          for (ListDigraph::OutArcIt b(graph, u); b != INVALID; ++b) {
            if (solution[b] > 0) {
              int sink_idx = std::distance(sink_nodes.begin(), std::find(sink_nodes.begin(), sink_nodes.end(), u));
              int cluster_idx = std::distance(cluster_nodes.begin(), std::find(cluster_nodes.begin(), cluster_nodes.end(), graph.target(b)));
              sinks_[sink_idx].cluster_idx = cluster_idx;
            }
          }
        }
      }
    }
  }


  void Clustering::getClustersIdx(
    std::vector<std::vector<unsigned>>& newClusters) const
  {
    newClusters.clear();
    newClusters.resize(clusters_.size());
    for (size_t i = 0; i < clusters_.size(); ++i) {
      newClusters[i].resize(clusters_[i].size());
      for (unsigned j = 0; j < clusters_[i].size(); ++j) {
        newClusters[i][j] = clusters_[i][j]->sink_idx;
      }
    }
  }

  /* static */
  float Clustering::calcDist(const std::pair<float, float>& loc, const Sink* sink)
  {
    return calcDist(loc, { sink->x, sink->y });
  }

  /* static */
  float Clustering::calcDist(const std::pair<float, float>& loc1,
    const std::pair<float, float>& loc2)
  {
    return std::abs(loc1.first - loc2.first)
      + std::abs(loc1.second - loc2.second);
  }

}  // namespace cts::CKMeans
