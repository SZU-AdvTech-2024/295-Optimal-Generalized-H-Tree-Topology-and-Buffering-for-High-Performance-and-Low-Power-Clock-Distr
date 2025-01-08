#include"../design/Clock.h"
#include "KMeansClustering.h"
#include "MinCostFlow.h"
#include "Legalization.h"

#include <cmath>

class HTree {

public:
  HTree(cts::Clock& clock, std::vector<Point>& points, Legalization& legal)
    :clock_(clock), points_(points), legal_(legal) {
  }
  ~HTree() = default;

  cts::Clock buildHTree();
  Point getRootNode() {
    return root_node_;
  }

private:
  void findBoundingBox();

  void findBestTree();

  void findSequence(int direction, int num_level,
    std::vector<std::vector<int>>& w_len,
    std::vector<std::vector<int>>& h_len);

  void tryBuildTree(int direction,
    int lev,
    Point& par_buf,
    const std::vector<int>& seq,
    std::vector<Point>& htree_leafs,
    cts::Clock& tmp_clock,
    Legalization& tmp_legal);

  void findMinLatency(int num_level,
    int direction,
    std::vector<int>& seq,
    std::vector<std::vector<int>> w_len,
    std::vector<std::vector<int>> h_len);

  double calcDist(std::pair<float, float> p1, Point p2);

  void calculateAllPathLatencies(cts::Clock& clock,
    const std::string& bufName,
    double curr_latency,
    int buf_count,
    int curr_level,
    double& total_latency,
    int& path_count);

  void startLegalization(cts::Clock& clock,
    std::pair<float, float>& means,
    Legalization& legalization);

  cts::Clock clock_;
  cts::Clock min_clock_;
  std::vector<Point> points_;
  Legalization legal_;
  Legalization curr_legal_;
  double min_latency_;

  int bounding_box_w_, bounding_box_h_;
  int root_node_x_, root_node_y_;
  Point root_node_;



};