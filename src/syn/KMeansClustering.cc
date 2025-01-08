#include "KMeansClustering.h"
#include <cmath>
#include <limits>
#include <ctime>
#include <algorithm>

KMeansClustering::KMeansClustering() : numClusters_(0) {
}

void KMeansClustering::run(const std::vector<Point>& data, unsigned int numClusters, unsigned int maxIterations, int cap) {
  this->numClusters_ = numClusters;
  if (data.empty()) return;

  centroids_.resize(numClusters);
  initializeCentroids(data, numClusters);

  clusters_.resize(numClusters);
  bool changed = true;
  unsigned int iteration = 0;

  while (changed && iteration < maxIterations) {
    changed = false;
    iteration++;
    reassignClusters(data, changed, cap);
    updateCentroids();
  }
}

void KMeansClustering::initializeCentroids(const std::vector<Point>& data, unsigned int numClusters) {
  for (unsigned int i = 0; i < numClusters; ++i) {
    centroids_[i] = data[i % data.size()].coordinates;
  }
}

void KMeansClustering::reassignClusters(const std::vector<Point>& data, bool& changed, int cap) {
  clusters_.clear();
  clusters_.resize(numClusters_);

  std::vector<int> clusterSizes(numClusters_, 0);

  for (size_t i = 0; i < data.size(); ++i) {
    int bestCluster = findNearestCentroid(data[i].coordinates);
    if (cap > 0 && clusterSizes[bestCluster] >= cap) {
      bestCluster = findAlternativeCluster(clusterSizes, cap, bestCluster);
    }
    clusters_[bestCluster].push_back(data[i]);
    clusterSizes[bestCluster]++;
    changed = true; // Update the flag if any point is reassigned
  }
}

int KMeansClustering::findNearestCentroid(const std::pair<float, float>& point) {
  float minDistance = std::numeric_limits<float>::max();
  int closest = -1;
  for (size_t j = 0; j < centroids_.size(); ++j) {
    float distance = std::pow(point.first - centroids_[j].first, 2) + std::pow(point.second - centroids_[j].second, 2);
    if (distance < minDistance) {
      minDistance = distance;
      closest = j;
    }
  }
  return closest;
}

void KMeansClustering::updateCentroids() {
  std::vector<std::pair<float, float>> newCentroids(numClusters_, { 0, 0 });
  std::vector<int> counts(numClusters_, 0);

  for (int i = 0; i < numClusters_; ++i) {
    for (const auto& point : clusters_[i]) {
      newCentroids[i].first += point.coordinates.first;
      newCentroids[i].second += point.coordinates.second;
      counts[i]++;
    }
    if (counts[i] > 0) {
      centroids_[i].first = newCentroids[i].first / counts[i];
      centroids_[i].second = newCentroids[i].second / counts[i];
    }
  }
}

int KMeansClustering::findAlternativeCluster(const std::vector<int>& clusterSizes, int cap, int exclude) {
  int index = -1, minSize = std::numeric_limits<int>::max();
  for (int i = 0; i < clusterSizes.size(); ++i) {
    if (i != exclude && clusterSizes[i] < minSize && clusterSizes[i] < cap) {
      minSize = clusterSizes[i];
      index = i;
    }
  }
  return index == -1 ? exclude : index;
}
