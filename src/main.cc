#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "syn/Clustering.h"
#include "syn/Legalization.h"
#include "syn/KMeansClustering.h"
#include "syn/GHTreeBuilder.h"
#include "syn/HTreeBuilder.h"
#include "Parser/ArgParser.h"
#include "design/Clock.h"
#include "syn/Partition.h"
#include "utility/Log.hh"
#include "utility/Time.hh"
#include "utility/Usage.hh"

using namespace cts::CKMeans;

// 读文件
void readInputFile(const std::string& filename1,
  const std::string& filename2, cts::Clock& clock)
{

  std::ifstream file1(filename1); // 创建文件流对象，用于读取文件
  std::ifstream file2(filename2);

  if (!file1)
  { // 检查文件是否成功打开
    std::cout << "Error: problem.def not found!" << std::endl;
    return; // 文件未成功打开，终止整个函数
  }
  if (!file2)
  { // 检查文件是否成功打开
    std::cout << "Error: constraint.txt not found!" << std::endl;
    return; // 文件未成功打开，终止整个函数
  }

  std::string line;     // 用来存储从文件中读取的行
  std::size_t startPos; // 被截取字符串的第一个位置
  std::size_t endPos;   // 被截取字符串的最后一个位置
  std::string data;     // 使用字符串流来转换字符串为数字

  int unit_distance_micron, die_area_x, die_area_y, max_fanout;
  int ff_w, ff_h, buf_w, buf_h, clk_x, clk_y;
  double net_unit_r, net_unit_c, max_net_rc, buffer_delay;

  // 读UNIT_DISTANCE_MICRON
  std::getline(file1, line);
  startPos = line.rfind('S');
  data = line.substr(startPos + 1);
  std::istringstream iss1(data);
  iss1 >> unit_distance_micron;
  clock.setUnitsDistanceMicrons(unit_distance_micron);

  // 读DIEAREA
  std::getline(file1, line);
  startPos = line.find('(');
  startPos = line.find('(', startPos + 1);
  startPos = line.find('(', startPos + 1);
  endPos = line.find(')', startPos);
  data = line.substr(startPos + 1, endPos - startPos - 1);
  std::istringstream iss2(data);
  iss2 >> die_area_x >> die_area_y;
  clock.setDieArea(die_area_x, die_area_y);

  // 读出第3行FF的宽和高
  std::getline(file1, line);
  startPos = line.find('(');
  startPos++;
  endPos = line.find(')', startPos);
  data = line.substr(startPos, endPos - startPos);
  std::istringstream iss3(data);
  iss3 >> ff_w >> ff_h;
  clock.setFFInfo(ff_w, ff_h);

  // 读BUF的宽和高
  std::getline(file1, line);
  startPos = line.find('(');
  startPos++;
  endPos = line.find(')', startPos);
  data = line.substr(startPos, endPos - startPos);
  std::istringstream iss4(data);
  iss4 >> buf_w >> buf_h;
  clock.setBufferInfo(buf_w, buf_h);

  // 读clk的坐标
  std::getline(file1, line);
  startPos = line.find('(');
  startPos++;
  endPos = line.find(')', startPos);
  data = line.substr(startPos, endPos - startPos);
  std::istringstream iss5(data);
  iss5 >> clk_x >> clk_y;
  clock.setClockPin(clk_x, clk_y);

  // 跳过第6行
  std::getline(file1, line);

  // 从第7行开始读取FF的坐标
  float ff_x, ff_y;
  std::string ff_name;
  while (std::getline(file1, line))
  {
    startPos = line.find("FF");
    if (startPos != std::string::npos)
    {
      endPos = line.find(' ', startPos);
      ff_name = line.substr(startPos, endPos - startPos);

      startPos = line.find('(', endPos) + 1;
      data = line.substr(startPos, line.find(')', startPos) - startPos);
      std::istringstream iss7(data);
      iss7 >> ff_x >> ff_y;
      cts::ClockInst tmp_FF(ff_name, cts::CLOCK_FF, ff_x, ff_y, ff_w, ff_h);
      clock.addFF(tmp_FF);
    }
  }

  // 读net_unit_r
  std::getline(file2, line);
  startPos = line.rfind('=');
  data = line.substr(startPos + 1);
  std::istringstream iss8(data);
  iss8 >> net_unit_r;
  clock.setNetUnitR(net_unit_r);

  // 读net_unit_c
  std::getline(file2, line);
  startPos = line.rfind('=');
  data = line.substr(startPos + 1);
  std::istringstream iss9(data);
  iss9 >> net_unit_c;
  clock.setNetUnitC(net_unit_c);

  // 读max_net_rc
  std::getline(file2, line);
  startPos = line.rfind('=');
  data = line.substr(startPos + 1);
  std::istringstream iss10(data);
  iss10 >> max_net_rc;
  clock.setMaxNetRC(max_net_rc);

  // 读maxfanout
  std::getline(file2, line);
  startPos = line.rfind('=');
  data = line.substr(startPos + 1);
  std::istringstream iss11(data);
  iss11 >> max_fanout;
  clock.setMaxFanout(max_fanout);

  // 读buffer_delay
  std::getline(file2, line);
  startPos = line.rfind('=');
  data = line.substr(startPos + 1);
  std::istringstream iss12(data);
  iss12 >> buffer_delay;
  clock.setBufferDelay(buffer_delay);

  file1.close();
  file2.close();
}

// 单元合法化
void startLegalization(cts::Clock& clock,
  std::pair<float, float>& means,
  Legalization& legalization)
{
  Rectangle tmp_rect{ means.first, means.second,
                     clock.getBuf_W(), clock.getBuf_H(), "BUF" };
  Rectangle rect = legalization.legalizeRectangle(tmp_rect);
  legalization.addFixedRectangle(rect);
  means.first = rect.x;
  means.second = rect.y;
}

// 单元合法化
void startLegalization(cts::Clock& clock,
  std::vector<std::pair<float, float>>& means,
  Legalization& legalization)
{
  for (size_t i = 0; i < means.size(); ++i)
  {
    Rectangle tmp_rect{ means[i].first, means[i].second,
                       clock.getBuf_W(), clock.getBuf_H(), "BUF" };
    Rectangle rect = legalization.legalizeRectangle(tmp_rect);
    legalization.addFixedRectangle(rect);
    means[i].first = rect.x;
    means[i].second = rect.y;
  }
}

void connectClkToRoot(cts::Clock& clock, const Point& root_node, Legalization& legal2) {
  // 起点(CLK)和终点(root)的坐标
  float x1 = clock.getClockPinX();
  float y1 = clock.getClockPinY();
  float x2 = root_node.coordinates.first;
  float y2 = root_node.coordinates.second;

  // 使用曼哈顿距离
  double L = std::abs(x2 - x1) + std::abs(y2 - y1);
  double d = clock.getBufferDelay();
  double r = clock.getNetUnitR();
  double c = clock.getNetUnitC();

  int current_buffer_count = clock.getNumBuffer();

  if (d >= 0.69 * r * c * L * L) {
    // 直接连接CLK和根节点
    cts::ClockSubNet& clk_net = clock.addSubNet("net_clk");
    cts::ClockInst tmp_clk("CLK");
    cts::ClockInst tmp_rootBuf("BUF" + std::to_string(root_node.index),
      cts::CLOCK_BUFFER,
      root_node.coordinates.first,
      root_node.coordinates.second,
      clock.getBuf_W(),
      clock.getBuf_H());
    clk_net.addInst(tmp_clk);
    clk_net.addInst(tmp_rootBuf);
  }
  else {
    // 计算理论最优buffer数量
    double optimal_n_float = L * std::sqrt(0.69 * r * c / d) - 1;
    int floor_n = std::floor(optimal_n_float);
    int ceil_n = std::ceil(optimal_n_float);

    // 计算延迟的lambda函数
    auto calcLatency = [&](int n) {
      return (n + 1) * d + 0.69 * r * c * L * L / (n + 1);
      };

    // 计算向上向下取整的延迟
    double floor_latency = calcLatency(floor_n);
    double ceil_latency = calcLatency(ceil_n);

    // 选择延迟最小的方案
    int final_n = (floor_latency <= ceil_latency) ? floor_n : ceil_n;
    double final_latency = (floor_latency <= ceil_latency) ? floor_latency : ceil_latency;

    // 计算x和y方向上的步长
    float delta_x = (x2 - x1) / (final_n + 1);
    float delta_y = (y2 - y1) / (final_n + 1);

    // 插入buffer并连接
    std::vector<cts::ClockInst> buffers;
    for (int i = 0; i < final_n; i++) {
      // 按比例计算每个buffer的坐标
      float x = x1 + (i + 1) * delta_x;
      float y = y1 + (i + 1) * delta_y;

      // 进行合法化
      std::pair<float, float> coords = { x, y };
      startLegalization(clock, coords, legal2);
      x = coords.first;
      y = coords.second;

      // buffer编号从大到小
      int buf_index = current_buffer_count + final_n - i;

      cts::ClockInst new_buf("BUF" + std::to_string(buf_index),
        cts::CLOCK_BUFFER,
        x, y,
        clock.getBuf_W(),
        clock.getBuf_H());
      buffers.push_back(new_buf);
      clock.addClockBuffer(new_buf);
    }

    // 更新clock中的buffer总数
    clock.setNumBuffer(current_buffer_count + final_n);

    // 创建CLK到第一个buffer的连接
    cts::ClockSubNet& clk_net = clock.addSubNet("net_clk");
    cts::ClockInst tmp_clk("CLK");
    clk_net.addInst(tmp_clk);
    clk_net.addInst(buffers[0]);

    // 连接中间的buffers
    for (size_t i = 0; i < buffers.size() - 1; i++) {
      std::string subnet_name = "net_buf" + std::to_string(current_buffer_count + final_n - i);
      cts::ClockSubNet& subnet = clock.addSubNet(subnet_name);
      subnet.addInst(buffers[i]);
      subnet.addInst(buffers[i + 1]);
    }

    // 连接最后一个buffer到root
    std::string last_subnet_name = "net_buf" + std::to_string(current_buffer_count + 1);
    cts::ClockSubNet& last_subnet = clock.addSubNet(last_subnet_name);
    last_subnet.addInst(buffers.back());
    cts::ClockInst tmp_rootBuf("BUF" + std::to_string(root_node.index),
      cts::CLOCK_BUFFER,
      root_node.coordinates.first,
      root_node.coordinates.second,
      clock.getBuf_W(),
      clock.getBuf_H());
    last_subnet.addInst(tmp_rootBuf);
  }
}

// 聚类
void startClustering(cts::Clock& clock)
{
  issta::Stats stats;

  // 底层FF聚类
  Partition partition(clock);
  partition.run();
  clock = partition.getClock();

  // 初始化下一层聚类的数据
  std::vector<std::pair<float, float>> points;
  points = partition.getBufPoints();
  std::vector<int> points_idx;
  int idx = 1;
  for (auto& point : points)
  {
    points_idx.push_back(idx++);
  }

  // 初始化单元合法化
  Legalization legalization2(clock.getDieArea_X(), clock.getDieArea_Y());
  for (cts::ClockInst& ff : clock.getFFs())
  {
    Rectangle rect{ ff.getX(), ff.getY(),
                   clock.getFF_W(), clock.getFF_H(), "FF" };
    legalization2.addFixedRectangle(rect);
  }
  for (auto& point : points)
  {
    Rectangle rect{ point.first, point.second,
                   clock.getBuf_W(), clock.getBuf_H(), "BUF" };
    legalization2.addFixedRectangle(rect);
  }

  // 不断聚类直到只剩一个簇
  while (points.size() > 1)
  {
    std::vector<std::vector<Sink*>> clusters;
    std::vector<std::pair<float, float>> means;

    
    if (points.size() < 256)
    {
      std::vector<Point> data;
      data.resize(points.size());
      for (size_t i = 0; i < points.size(); i++) {
        data[i].coordinates.first = points[i].first;
        data[i].coordinates.second = points[i].second;
        data[i].index = points_idx[i];
      }

      // 初始化GHTree
      HTree htree(clock, data, legalization2);

      // 构建GHTree
      clock = htree.buildHTree();

      // 获取H树的根节点
      Point root_node = htree.getRootNode();

      // 将rootnode和clk相连
      connectClkToRoot(clock, root_node, legalization2);

      break;
    }
    else
    {
      int fanout = clock.getMaxFanout();
      bool flag = true;

      // 初始化 Clustering
      Clustering clustering(points, points_idx,
        clock.getNetUnitR(), clock.getNetUnitC(), clock.getMaxNetRC());

      do
      {
        flag = true;

        // 设置簇的数量
        double num = static_cast<double>(points.size()) / fanout;
        unsigned numClusters = static_cast<unsigned>(ceil(num));

        // 初始化聚类中心点
        means.clear();
        for (unsigned i = 0; i < numClusters; ++i)
        {
          means.push_back(points[i % points.size()]); // 选择初始中心点
        }

        // 执行聚类
        clustering.iterKmeans(3, numClusters, fanout, 100, 2,
          means); // 使用 max_fanout 限制每个簇的最大容量

        // 保存聚类结果
        clusters.clear();
        clusters = clustering.getCluster();

        // 单元合法化
        startLegalization(clock, means, legalization2);

        // 计算根据当前fanout聚类，rc会不会超
        for (size_t i = 0; i < means.size(); ++i)
        {
          cts::ClockSubNet tmp_subnet("tmp_subnet");
          cts::ClockInst tmp_newBuf("new_BUF", cts::CLOCK_BUFFER, means[i].first,
            means[i].second, clock.getBuf_W(), clock.getBuf_H());
          tmp_subnet.addInst(tmp_newBuf);
          for (Sink* sink : clusters[i])
          {
            cts::ClockInst tmp_oldBuf("old_BUF", cts::CLOCK_BUFFER,
              sink->x, sink->y, clock.getBuf_W(), clock.getBuf_H());
            tmp_subnet.addInst(tmp_oldBuf);
          }
          float net_rc = clock.calcNetRC(tmp_subnet);
          if (net_rc > clock.getMaxNetRC())
          {
            flag = false;
            fanout--;
            means.clear();
            clusters.clear();
            break;
          }
        }
      } while (flag == false);

      // 保存结果
      for (size_t i = 0; i < means.size(); ++i)
      {
        // 创建一个子网
        int numBuffer = clock.getNumBuffer();
        numBuffer += 1;
        clock.setNumBuffer(numBuffer);
        std::string subnetName = "net_buf" + std::to_string(numBuffer);
        cts::ClockSubNet& subnet = clock.addSubNet(subnetName);

        // 先往子网中放入一个簇的中心点
        cts::ClockInst tmp_newBuf("BUF" + std::to_string(numBuffer),
          cts::CLOCK_BUFFER, means[i].first,
          means[i].second, clock.getBuf_W(), clock.getBuf_H());
        subnet.addInst(tmp_newBuf);
        clock.addClockBuffer(tmp_newBuf);

        // 再往子网中放入这个簇的其他点
        for (Sink* sink : clusters[i])
        {
          cts::ClockInst tmp_Buf("BUF" + std::to_string(sink->sink_idx),
            cts::CLOCK_BUFFER, sink->x, sink->y, clock.getBuf_W(), clock.getBuf_H());
          subnet.addInst(tmp_Buf);
        }

        // clock.calcNetRC(subnetName);
      }

      // 更新下一轮聚类的初始数据
      points = means;
      points_idx.clear();
      for (int i = clock.getNumBuffer() - means.size(); i < clock.getNumBuffer(); i++)
      {
        points_idx.emplace_back(i + 1);
      }
      means.clear();
    }
  }
}

// 输出文件
void printOutputFile(const std::string& inputPath,
  const std::string& outputPath, cts::Clock& clock)
{

  std::ifstream inputFile(inputPath);
  std::vector<std::string> lines;
  std::string line;

  // 复制输入文件
  while (getline(inputFile, line))
  {
    lines.push_back(line);
  }
  inputFile.close();

  // 修改component的数量
  lines[5] = "COMPONENTS " +
    std::to_string(clock.getFFs().size() + clock.getNumBuffer()) + " ;";

  // 添加buffer
  std::vector<std::string> bufLines;
  for (auto& buffer : clock.getBuffers())
  {
    std::string tmp_line = "- " + buffer.getName() + " BUF ( " + std::to_string(buffer.getX()) + " " + std::to_string(buffer.getY()) + " ) ;";
    bufLines.push_back(tmp_line);
  }
  lines.insert(lines.begin() + 6 + clock.getFFs().size(),
    bufLines.begin(), bufLines.end());

  // 添加net
  lines.push_back("NETS " + std::to_string(clock.getSubNets().size()) + " ;");
  for (auto& subnet : clock.getSubNets())
  {
    std::deque<cts::ClockInst> instances = subnet.getClockInst();
    std::string tmp_line = "- " + subnet.getName() + " ( " + instances[0].getName() + " ) (";
    for (int i = 1; i < instances.size(); i++)
    {
      tmp_line = tmp_line + " " + instances[i].getName();
    }
    tmp_line += " ) ;";
    lines.push_back(tmp_line);
  }
  lines.push_back("END NETS ;");

  // 输出文件
  std::ofstream outputFile(outputPath);
  for (int i = 0; i < lines.size() - 1; i++)
  {
    outputFile << lines[i] << "\n";
  }
  outputFile << lines[lines.size() - 1];
  outputFile.close();
}

int main(int argc, char** argv)
{
  LOG_INFO("start run iCTS at %s", GetCurrWallTime().c_str());
  issta::Stats stats;
  ArgParser arg_parser;
  arg_parser.parseCmd(argc, argv);
  std::string problem_file = arg_parser.get_problem_file();
  std::string constraint_file = arg_parser.get_constraint_file();
  std::string output_file = arg_parser.get_output_file();

  cts::Clock clock("clk", -1, -1);

  // // 读赛题大样例
  readInputFile(problem_file, constraint_file, clock);

  // readInputFile("/home/sujianrong/edacontest-cts/exp/pre_submit/case2/problem.def",
  //   "/home/sujianrong/edacontest-cts/exp/pre_submit/case2/constraints.txt",
  //   clock);


  startClustering(clock);

  // 输出赛题def
  printOutputFile(problem_file, output_file, clock);

  // printOutputFile("/home/sujianrong/edacontest-cts/exp/pre_submit/case2/problem.def",
  //   "/home/sujianrong/edacontest-cts/exp/pre_submit/case2/solution.def", clock);


  LOG_INFO("finish run iCTS at %s", GetCurrWallTime().c_str());
  double memory_delta = stats.memoryDelta();
  LOG_INFO("memory usage: %f", memory_delta);
  double time_delta = stats.elapsedRunTime();
  LOG_INFO("time elapsed: %f", time_delta);

  return 0;
}
