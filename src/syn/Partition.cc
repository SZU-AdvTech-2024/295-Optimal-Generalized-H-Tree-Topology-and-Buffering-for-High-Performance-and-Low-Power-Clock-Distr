#include "Partition.h"
#include <cmath>
#include "design/Clock.h"

void Partition::run()
{
  divideDiearea();
}

// 单元合法化
void Partition::startLegalization_tmp(cts::Clock &clock,
                                      Legalization &tmp_legalization,
                                      std::vector<std::pair<float, float>> &means)
{

  for (size_t i = 0; i < means.size(); ++i)
  {
    Rectangle tmp_rect{round(means[i].first), round(means[i].second),
                       clock.getBuf_W(), clock.getBuf_H(), "BUF"};
    Rectangle rect = tmp_legalization.legalizeRectangle(tmp_rect);
    tmp_legalization.addFixedRectangle(rect);
    means[i].first = rect.x;
    means[i].second = rect.y;
  }
}

// 版图划分
void Partition::divideDiearea()
{

  issta::Stats stats;

  // 初始化
  legalization_.setDiearea(clock_.getDieArea_X(), clock_.getDieArea_Y());
  std::vector<std::pair<float, float>> data;

  for (cts::ClockInst &ff : clock_.getFFs())
  {
    // 初始化data
    data.emplace_back(ff.getX(), ff.getY());
    // 初始化单元合法化的FF
    Rectangle rect{ff.getX(), ff.getY(),
                   clock_.getFF_W(), clock_.getFF_H(), ff.getName()};
    legalization_.addFixedRectangle(rect);
  }

  col_ = 5;
  bool flag = true;

  do
  {
    flag = true;

    // 大矩形划分成小矩形，初始化信息
    width_ = clock_.getDieArea_X();
    height_ = clock_.getDieArea_Y();
    float ratio = static_cast<float>(height_) / width_;
    row_ = floor(ratio * col_);
    float w = static_cast<float>(width_) / col_;  // 小矩形的宽
    float h = static_cast<float>(height_) / row_; // 小矩形的高
    int numGroups = row_ * col_;                  // 总共划分的组数
    groups_.clear();
    groupIndexs_.clear();
    groups_.resize(numGroups);
    groupIndexs_.resize(numGroups);


    // 对FF进行划分
    for (int i = 0; i < data.size(); ++i)
    {
      int tmp_num_col, tmp_num_row, idx;
      // 计算每个FF应该放在哪个小矩形中
      if (data[i].first == width_)
      {
        tmp_num_row = data[i].second / h;
        idx = (tmp_num_row + 1) * col_ - 1;
      }
      else if (data[i].second == height_)
      {
        tmp_num_col = data[i].first / w;
        idx = (row_ - 1) * col_ + tmp_num_col;
      }
      else
      {
        tmp_num_col = data[i].first / w;
        tmp_num_row = data[i].second / h;
        idx = tmp_num_row * col_ + tmp_num_col;
      }
      // FF放入对应的矩形中
      groups_[idx].push_back(data[i]);
      groupIndexs_[idx].push_back(i);
    }

    std::vector<std::vector<std::vector<Point>>> res_clusters(groups_.size());
    std::vector<std::vector<std::pair<float, float>>> res_centroids(groups_.size());

    int cnt = 1;
    Legalization tmp_legalization = legalization_;

    // 对每个组进行试探性聚类，检查这样划分net_rc是否会超过max_rc
    for (size_t i = 0; i < groups_.size(); ++i)
    {
     
      cts::Clock tmp_clock;
      tmp_clock = clock_;

      // 初始化数据
      std::vector<Point> tmp_data;
      tmp_data.resize(groups_[i].size());

      for (size_t j = 0; j < groups_[i].size(); ++j)
      {
        tmp_data[j].coordinates.first = groups_[i][j].first;
        tmp_data[j].coordinates.second = groups_[i][j].second;
        tmp_data[j].index = groupIndexs_[i][j];
      }

      std::vector<std::vector<Point>> tmp_clusters;
      std::vector<std::pair<float, float>> tmp_centroids;

      if (tmp_data.size() != 0)
      {
        // 找出聚类中心
        KMeansClustering kmeans;
        double num = static_cast<double>(groups_[i].size()) / clock_.getMaxFanout();
        unsigned int numClusters = static_cast<unsigned int>(ceil(num));
        kmeans.run(tmp_data, numClusters, 100, clock_.getMaxFanout());
        tmp_clusters = kmeans.getClusters();
        tmp_centroids = kmeans.getCentroids();

        // 单元合法化
        startLegalization_tmp(tmp_clock, tmp_legalization, tmp_centroids);

        // 检查每个簇的net_rc是否超过max_rc
        for (int j = 0; j < tmp_clusters.size(); ++j)
        {

          // 创建子网
          cts::ClockSubNet tmp_subnet("tmp_subnet");
          cts::ClockInst tmp_newBuf("BUF", cts::CLOCK_BUFFER, tmp_centroids[j].first,
                                    tmp_centroids[j].second, clock_.getBuf_W(), clock_.getBuf_H());
          tmp_subnet.addInst(tmp_newBuf);

          for (int k = 0; k < tmp_clusters[j].size(); ++k)
          {
            cts::ClockInst tmp_FF("FF", cts::CLOCK_FF, tmp_clusters[j][k].coordinates.first,
                                  tmp_clusters[j][k].coordinates.second, clock_.getFF_W(), clock_.getFF_H());
            tmp_subnet.addInst(tmp_FF);
            
          }

          // 计算子网的net_rc
          float net_rc = calcNetRC(tmp_subnet);

          // 如果超出max_rc，重新划分
          if (net_rc > clock_.getMaxNetRC())
          {
            flag = false;
            break;
          }
        }
      }
      // 保存每个组的聚类结果
      res_clusters[i] = tmp_clusters;
      res_centroids[i] = tmp_centroids;
    }

    // 如果rc都满足条件，则保存结果，更新clock
    if (flag == true)
    {

      // 保存合法化
      legalization_ = tmp_legalization;

      int buf_idx = 1;
      // 保存每个组的结果
      for (size_t i = 0; i < groups_.size(); ++i)
      {
        for (size_t j = 0; j < res_clusters[i].size(); ++j)
        {
          // 创建一个子网
          int numBuffer = clock_.getNumBuffer();
          numBuffer += 1;
          clock_.setNumBuffer(numBuffer);
          std::string subnetName = "net_buf" + std::to_string(numBuffer);
          cts::ClockSubNet &subnet = clock_.addSubNet(subnetName);

          // 先往子网中放入一个簇的中心点
          cts::ClockInst tmp_newBuf("BUF" + std::to_string(numBuffer),
                                    cts::CLOCK_BUFFER, res_centroids[i][j].first,
                                    res_centroids[i][j].second, clock_.getBuf_W(), clock_.getBuf_H());
          subnet.addInst(tmp_newBuf);
          clock_.addClockBuffer(tmp_newBuf);

          // 收集中心点，作为下一层聚类的数据
          std::pair<float, float> buf_point;
          buf_point.first = res_centroids[i][j].first;
          buf_point.second = res_centroids[i][j].second;
          buf_points.push_back(buf_point);
          // Point buf_point;
          // buf_point.coordinates.first = res_centroids[i][j].first;
          // buf_point.coordinates.second = res_centroids[i][j].second;
          // buf_point.index = buf_idx;
          // buf_idx++;
          // buf_points.push_back(buf_point);

          // 再往子网中放入这个簇的其他点
          for (size_t k = 0; k < res_clusters[i][j].size(); ++k)
          {
            cts::ClockInst tmp_FF("FF" + std::to_string(res_clusters[i][j][k].index + 1),
                                  cts::CLOCK_FF, res_clusters[i][j][k].coordinates.first,
                                  res_clusters[i][j][k].coordinates.second, clock_.getFF_W(), clock_.getFF_H());
            subnet.addInst(tmp_FF);
          }
        }

      }
    }
    else
    {
      col_++;
    }

  } while (flag == false);
}

float Partition::calcNetRC(cts::ClockSubNet &subnet)
{
  std::deque<cts::ClockInst> clockInst = subnet.getClockInst();
  cts::ClockInst buffer = clockInst[0];
  double half_rc = 0.5 * clock_.getNetUnitR() * clock_.getNetUnitC();
  double net_rc = 0.0;
  for (int i = 1; i < clockInst.size(); ++i)
  {
    double rc = half_rc * pow(calcDist(buffer, clockInst[i]), 2);
    net_rc += rc;
  }
  return net_rc;
}

int Partition::calcDist(cts::ClockInst inst1, cts::ClockInst inst2)
{
  // 计算曼哈顿距离
  return std::abs((inst1.getX() + inst1.getW() / 2) - (inst2.getX() + inst2.getW() / 2)) + std::abs((inst1.getY() + inst1.getH() / 2) - (inst2.getY() + inst2.getH() / 2));
}