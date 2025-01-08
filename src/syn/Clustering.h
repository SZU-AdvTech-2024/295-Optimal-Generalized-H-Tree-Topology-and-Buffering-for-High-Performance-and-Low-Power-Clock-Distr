#include <cmath>
#include <string>
#include <vector>

namespace cts::CKMeans {

  struct Sink
  {
    Sink(const float x, const float y, const unsigned idx)
      : x(x), y(y), cluster_idx(-1), sink_idx(idx) {
    };

    const float x, y;
    int cluster_idx;
    const unsigned sink_idx;  // index in sinks_
  };

  class Clustering
  {
  public:
    Clustering(const std::vector<std::pair<float, float>>& sinks,
      std::vector<int>& sinks_idx,double r, double c, double rc);
    ~Clustering();

    void iterKmeans(unsigned iter,
      unsigned n,
      unsigned cap,
      unsigned max,
      unsigned power,
      std::vector<std::pair<float, float>>& means);

    void getClustersIdx(std::vector<std::vector<unsigned>>& newClusters) const;
    std::vector<std::vector<Sink*>> getCluster() { return clusters_; }

  private:
    float Kmeans(unsigned n,
      unsigned cap,
      unsigned max,
      unsigned power,
      std::vector<std::pair<float, float>>& means);
    float calcSilh(const std::vector<std::pair<float, float>>& means) const;
    void minCostFlow(const std::vector<std::pair<float, float>>& means,
      unsigned cap,
      double max_net_rc,
      unsigned power);

    static float calcDist(const std::pair<float, float>& loc, const Sink* sink);
    static float calcDist(const std::pair<float, float>& loc1,
      const std::pair<float, float>& loc2);

    std::vector<Sink> sinks_;
    std::vector<std::vector<Sink*>> clusters_;
    double net_unit_r_;              // 单位长度电阻值r
    double net_unit_c_;              // 电容值c
    double max_net_rc_;
  };

}  // namespace cts::CKMeans
