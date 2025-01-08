#include <vector>
#include "../design/Clock.h"
#include "Legalization.h"
#include "../utility/Usage.hh"
#include "../utility/Log.hh"
#include "../utility/Time.hh"
#include "KMeansClustering.h"

class Partition
{

public:
  Partition(cts::Clock &clock) : clock_(clock) {}
  ~Partition() = default;

  void run();

  // getter
  const std::vector<std::vector<std::pair<float, float>>> &getGroups() const
  {
    return groups_;
  }
  const std::vector<std::vector<int>> &getGroupIndexs() const
  {
    return groupIndexs_;
  }
  cts::Clock &getClock()
  {
    return clock_;
  }
  std::vector<std::pair<float, float>> &getBufPoints()
  {
    return buf_points;
  }
  Legalization &getLegalization()
  {
    return legalization_;
  }

private:
  void startLegalization_tmp(cts::Clock &clock,
                             Legalization &tmp_legalization,
                             std::vector<std::pair<float, float>> &means);
  void divideDiearea();

  float calcNetRC(cts::ClockSubNet &subnet);
  int calcDist(cts::ClockInst inst1, cts::ClockInst inst2);

  cts::Clock clock_;
  Legalization legalization_;
  int width_ = 0, height_ = 0;                               // 待划分区域的宽和高
  int row_ = 0, col_ = 0;                                    // 划分的行数和列数
  std::vector<std::vector<std::pair<float, float>>> groups_; // 记录每个组有哪些FF
  std::vector<std::vector<int>> groupIndexs_;                // 表示FF在哪个组中
  std::vector<std::pair<float, float>> buf_points;
};