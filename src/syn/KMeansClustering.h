#ifndef KMEANS_CLUSTERING_H
#define KMEANS_CLUSTERING_H

#include <vector>
#include <utility> 
#include "../design/Clock.h"


// KMeans聚类类定义
class KMeansClustering {
public:
  KMeansClustering();
  void run(const std::vector<Point>& data,
    unsigned int numClusters,
    unsigned int maxIterations = 100,
    int cap = -1);

  const std::vector<std::vector<Point>>& getClusters() const {
    return clusters_;
  }  // 获取聚类的点
  const std::vector<std::pair<float, float>>& getCentroids()const {
    return centroids_;
  }  // 获取聚类中心

private:
  std::vector<std::pair<float, float>> centroids_; // 簇的中心点列表
  std::vector<std::vector<Point>> clusters_; // 每个簇的点
  unsigned int numClusters_; // 簇的数量

  void initializeCentroids(const std::vector<Point>& data,
    unsigned int numClusters);

  void reassignClusters(const std::vector<Point>& data,
    bool& changed, int cap);

  int findNearestCentroid(const std::pair<float, float>& point);

  void updateCentroids();

  int findAlternativeCluster(const std::vector<int>& clusterSizes,
    int cap, int exclude);
};

#endif // KMEANS_CLUSTERING_H
