#include "HTreeBuilder.h"

cts::Clock HTree::buildHTree()
{
  findBoundingBox();
  findBestTree();
  return clock_;
}

void HTree::findBoundingBox()
{
  // 记录叶子结点的坐标
  int min_x, min_y, max_x, max_y;
  min_x = max_x = points_[0].coordinates.first;
  min_y = max_y = points_[0].coordinates.second;

  // 找bounding_box
  for (size_t i = 1; i < points_.size(); ++i)
  {
    if (points_[i].coordinates.first < min_x)
    {
      min_x = points_[i].coordinates.first;
    }
    if (points_[i].coordinates.second < min_y)
    {
      min_y = points_[i].coordinates.second;
    }

    if (points_[i].coordinates.first > max_x)
    {
      max_x = points_[i].coordinates.first;
    }
    if (points_[i].coordinates.second > max_y)
    {
      max_y = points_[i].coordinates.second;
    }
  }
  bounding_box_w_ = max_x - min_x;
  bounding_box_h_ = max_y - min_y;

  // 找根节点
  root_node_x_ = min_x + bounding_box_w_ / 2;
  root_node_y_ = min_y + bounding_box_h_ / 2;
}

void HTree::findBestTree()
{
  int num_level = floor(log2(points_.size())) - 1;
  min_latency_ = std::numeric_limits<double>::max();

  // 存储每一层的宽和高
  std::vector<std::vector<int>> w_len(num_level, std::vector<int>(5));
  std::vector<std::vector<int>> h_len(num_level, std::vector<int>(5));

  // H树构建序列
  std::vector<int> seq; // debug

  // 从横向开始找w的长度
  findSequence(0, num_level, w_len, h_len);
  findMinLatency(num_level, 0, seq, w_len, h_len);

  // 从纵向开始找h的长度
  w_len.clear();
  h_len.clear();
  seq.clear(); // 清空序列
  w_len.resize(num_level, std::vector<int>(5));
  h_len.resize(num_level, std::vector<int>(5));

  findSequence(1, num_level, w_len, h_len);
  findMinLatency(num_level, 1, seq, w_len, h_len);

  clock_ = min_clock_;
  legal_ = curr_legal_;
}

void HTree::findSequence(int direction, int num_level,
                         std::vector<std::vector<int>> &w_len,
                         std::vector<std::vector<int>> &h_len)
{
  // 动态更新bounding_box的w和h
  std::vector<int> tmp_box_w;
  std::vector<int> tmp_box_h;

  tmp_box_w.push_back(bounding_box_w_);
  tmp_box_h.push_back(bounding_box_h_);

  int cnt_w = 0, cnt_h = 0;
  for (int i = 0; i < num_level; i++)
  {
    direction++;
    // 横向w
    if (direction % 2 == 1)
    {
      // 只使用一种长度：当前box宽度的一半
      w_len[cnt_w][0] = tmp_box_w[cnt_w] / 2;
      // 更新bounding_box的w
      cnt_w++;
      tmp_box_w.push_back(w_len[cnt_w - 1][0]);
    }
    // 纵向h
    else
    {
      // 只使用一种长度：当前box高度的一半
      h_len[cnt_h][0] = tmp_box_h[cnt_h] / 2;
      // 更新bounding_box的h
      cnt_h++;
      tmp_box_h.push_back(h_len[cnt_h - 1][0]);
    }
  }
}

void HTree::tryBuildTree(int direction,
                         int lev,
                         Point &par_buf,
                         const std::vector<int> &seq,
                         std::vector<Point> &htree_leafs,
                         cts::Clock &tmp_clock,
                         Legalization &tmp_legal)
{

  if (lev == seq.size())
  {
    // 添加叶子节点到htree_leafs列表
    htree_leafs.push_back(par_buf);

    // 创建叶子节点的buffer
    std::string leafBufName = "BUF" + std::to_string(par_buf.index);
    cts::ClockInst leaf_buf(leafBufName,
                            cts::CLOCK_BUFFER,
                            par_buf.coordinates.first,
                            par_buf.coordinates.second,
                            tmp_clock.getBuf_W(),
                            tmp_clock.getBuf_H());

    tmp_clock.addClockBuffer(leaf_buf);
    return;
  }

  // check H tree branch opt
  double L = seq[lev] / 2;
  double d = tmp_clock.getBufferDelay();
  double r = tmp_clock.getNetUnitR();
  double c = tmp_clock.getNetUnitC();
  bool h_opt = d < 1.38 * r * c * L * L;
  size_t n = h_opt ? std::round(L * std::sqrt(1.38 * r * c / d) - 1) : 1;
  if (h_opt)
  {
    n = n == 0 ? 1 : n;
  }

  // 横向
  if (direction % 2 == 0)
  {
    // 对父节点进行合法化，使用临时的legalization对象
    std::pair<float, float> par_coord = {par_buf.coordinates.first, par_buf.coordinates.second};
    startLegalization(tmp_clock, par_coord, tmp_legal);
    par_buf.coordinates.first = par_coord.first;
    par_buf.coordinates.second = par_coord.second;

    // 创建子网
    std::string subnetName = "net_buf" + std::to_string(par_buf.index);
    cts::ClockSubNet &subnet = tmp_clock.addSubNet(subnetName);
    cts::ClockInst tmp_newBuf("BUF" + std::to_string(par_buf.index),
                              cts::CLOCK_BUFFER, par_buf.coordinates.first,
                              par_buf.coordinates.second, tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
    subnet.addInst(tmp_newBuf);
    tmp_clock.addClockBuffer(tmp_newBuf);

    // 将父节点buffer添加到临时legalization中
    Rectangle par_rect{static_cast<int>(par_buf.coordinates.first),
                       static_cast<int>(par_buf.coordinates.second),
                       tmp_clock.getBuf_W(),
                       tmp_clock.getBuf_H(),
                       "BUF"};
    tmp_legal.addFixedRectangle(par_rect);

    // h tree branch opt
    if (h_opt)
    {
      // insert n buffer from tmp_newBuf to buf1
      auto pre_loc1 = par_buf;
      auto *pre_sub_net1 = &subnet;

      Point buf1;
      buf1.coordinates.first = par_buf.coordinates.first - seq[lev] / 2;
      buf1.coordinates.second = par_buf.coordinates.second;

      auto delta_x = (buf1.coordinates.first - pre_loc1.coordinates.first) / (n + 1);
      auto delta_y = (buf1.coordinates.second - pre_loc1.coordinates.second) / (n + 1);
      for (size_t i = 0; i < n; ++i)
      {
        Point temp_opt_buf_loc;
        temp_opt_buf_loc.coordinates.first = pre_loc1.coordinates.first + delta_x;
        temp_opt_buf_loc.coordinates.second = pre_loc1.coordinates.second + delta_y;
        temp_opt_buf_loc.index = tmp_clock.getNumBuffer() + 1;
        tmp_clock.setNumBuffer(temp_opt_buf_loc.index);
        // leagalize
        std::pair<float, float> temp_opt_coord = {temp_opt_buf_loc.coordinates.first, temp_opt_buf_loc.coordinates.second};
        startLegalization(tmp_clock, temp_opt_coord, tmp_legal);
        cts::ClockInst temp_opt_Buf("BUF" + std::to_string(temp_opt_buf_loc.index),
                                    cts::CLOCK_BUFFER, temp_opt_coord.first, temp_opt_coord.second,
                                    tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
        tmp_legal.addFixedRectangle(Rectangle{static_cast<int>(temp_opt_buf_loc.coordinates.first),
                                              static_cast<int>(temp_opt_buf_loc.coordinates.second),
                                              tmp_clock.getBuf_W(),
                                              tmp_clock.getBuf_H(),
                                              "BUF"});
        pre_sub_net1->addInst(temp_opt_Buf);
        // update pre*
        pre_loc1 = temp_opt_buf_loc;
        auto temp_sub_net_name = "net_buf" + std::to_string(temp_opt_buf_loc.index);
        auto &new_net = tmp_clock.addSubNet(temp_sub_net_name);
        pre_sub_net1 = &new_net;
        pre_sub_net1->addInst(temp_opt_Buf);
        tmp_clock.addClockBuffer(temp_opt_Buf);
      }
      // connect buf1
      buf1.index = tmp_clock.getNumBuffer() + 1;
      tmp_clock.setNumBuffer(buf1.index);
      // leagalize
      std::pair<float, float> buf1_coord = {buf1.coordinates.first, buf1.coordinates.second};
      startLegalization(tmp_clock, buf1_coord, tmp_legal);
      cts::ClockInst tmp_Buf1("BUF" + std::to_string(buf1.index),
                              cts::CLOCK_BUFFER, buf1.coordinates.first, buf1.coordinates.second,
                              tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      tmp_legal.addFixedRectangle(Rectangle{static_cast<int>(buf1.coordinates.first),
                                            static_cast<int>(buf1.coordinates.second),
                                            tmp_clock.getBuf_W(),
                                            tmp_clock.getBuf_H(),
                                            "BUF"});
      pre_sub_net1->addInst(tmp_Buf1);

      // insert n buffer from tmp_newBuf to buf2
      auto pre_loc2 = par_buf;
      auto *pre_sub_net2 = &subnet;
      Point buf2;
      buf2.coordinates.first = par_buf.coordinates.first + seq[lev] / 2;
      buf2.coordinates.second = par_buf.coordinates.second;

      delta_x = (buf2.coordinates.first - pre_loc2.coordinates.first) / (n + 1);
      delta_y = (buf2.coordinates.second - pre_loc2.coordinates.second) / (n + 1);
      for (size_t i = 0; i < n; ++i)
      {
        Point temp_opt_buf_loc;
        temp_opt_buf_loc.coordinates.first = pre_loc2.coordinates.first + delta_x;
        temp_opt_buf_loc.coordinates.second = pre_loc2.coordinates.second + delta_y;
        temp_opt_buf_loc.index = tmp_clock.getNumBuffer() + 1;
        tmp_clock.setNumBuffer(temp_opt_buf_loc.index);
        // leagalize
        std::pair<float, float> temp_opt_coord = {temp_opt_buf_loc.coordinates.first, temp_opt_buf_loc.coordinates.second};
        startLegalization(tmp_clock, temp_opt_coord, tmp_legal);
        cts::ClockInst temp_opt_Buf("BUF" + std::to_string(temp_opt_buf_loc.index),
                                    cts::CLOCK_BUFFER, temp_opt_coord.first, temp_opt_coord.second,
                                    tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
        tmp_legal.addFixedRectangle(Rectangle{static_cast<int>(temp_opt_buf_loc.coordinates.first),
                                              static_cast<int>(temp_opt_buf_loc.coordinates.second),
                                              tmp_clock.getBuf_W(),
                                              tmp_clock.getBuf_H(),
                                              "BUF"});
        pre_sub_net2->addInst(temp_opt_Buf);
        // update pre*
        pre_loc2 = temp_opt_buf_loc;
        auto temp_sub_net_name = "net_buf" + std::to_string(temp_opt_buf_loc.index);
        auto &new_net = tmp_clock.addSubNet(temp_sub_net_name);
        pre_sub_net2 = &new_net;
        pre_sub_net2->addInst(temp_opt_Buf);
        tmp_clock.addClockBuffer(temp_opt_Buf);
      }
      // connect buf2
      buf2.index = tmp_clock.getNumBuffer() + 1;
      tmp_clock.setNumBuffer(buf2.index);
      // leagalize
      std::pair<float, float> buf2_coord = {buf2.coordinates.first, buf2.coordinates.second};
      startLegalization(tmp_clock, buf2_coord, tmp_legal);
      cts::ClockInst tmp_Buf2("BUF" + std::to_string(buf2.index),
                              cts::CLOCK_BUFFER, buf2.coordinates.first, buf2.coordinates.second,
                              tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      tmp_legal.addFixedRectangle(Rectangle{static_cast<int>(buf2.coordinates.first),
                                            static_cast<int>(buf2.coordinates.second),
                                            tmp_clock.getBuf_W(),
                                            tmp_clock.getBuf_H(),
                                            "BUF"});
      pre_sub_net2->addInst(tmp_Buf2);

      // recursive H tree
      tryBuildTree(direction + 1, lev + 1, buf1, seq, htree_leafs, tmp_clock, tmp_legal);
      tryBuildTree(direction + 1, lev + 1, buf2, seq, htree_leafs, tmp_clock, tmp_legal);
    }
    else
    { // 计算坐标
      Point buf1, buf2;
      buf1.coordinates.first = par_buf.coordinates.first - seq[lev] / 2;
      buf1.coordinates.second = par_buf.coordinates.second;
      buf1.index = tmp_clock.getNumBuffer() + 1;
      tmp_clock.setNumBuffer(buf1.index);
      buf2.coordinates.first = par_buf.coordinates.first + seq[lev] / 2;
      buf2.coordinates.second = par_buf.coordinates.second;
      buf2.index = tmp_clock.getNumBuffer() + 1;
      tmp_clock.setNumBuffer(buf2.index);

      // 对两个子节点进行合法化，使用临时的legalization对象
      std::pair<float, float> buf1_coord = {buf1.coordinates.first, buf1.coordinates.second};
      std::pair<float, float> buf2_coord = {buf2.coordinates.first, buf2.coordinates.second};
      startLegalization(tmp_clock, buf1_coord, tmp_legal);
      startLegalization(tmp_clock, buf2_coord, tmp_legal);
      buf1.coordinates.first = buf1_coord.first;
      buf1.coordinates.second = buf1_coord.second;
      buf2.coordinates.first = buf2_coord.first;
      buf2.coordinates.second = buf2_coord.second;

      // 添加到子网中
      cts::ClockInst tmp_Buf1("BUF" + std::to_string(buf1.index),
                              cts::CLOCK_BUFFER, buf1.coordinates.first, buf1.coordinates.second,
                              tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      subnet.addInst(tmp_Buf1);

      cts::ClockInst tmp_Buf2("BUF" + std::to_string(buf2.index),
                              cts::CLOCK_BUFFER, buf2.coordinates.first, buf2.coordinates.second,
                              tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      subnet.addInst(tmp_Buf2);

      // 将子节点buffer添加到临时legalization中
      Rectangle buf1_rect{static_cast<int>(buf1.coordinates.first),
                          static_cast<int>(buf1.coordinates.second),
                          tmp_clock.getBuf_W(),
                          tmp_clock.getBuf_H(),
                          "BUF"};
      tmp_legal.addFixedRectangle(buf1_rect);

      Rectangle buf2_rect{static_cast<int>(buf2.coordinates.first),
                          static_cast<int>(buf2.coordinates.second),
                          tmp_clock.getBuf_W(),
                          tmp_clock.getBuf_H(),
                          "BUF"};
      tmp_legal.addFixedRectangle(buf2_rect);

      // 递归构建H树
      tryBuildTree(direction + 1, lev + 1, buf1, seq, htree_leafs, tmp_clock, tmp_legal);
      tryBuildTree(direction + 1, lev + 1, buf2, seq, htree_leafs, tmp_clock, tmp_legal);
    }
  }
  // 纵向
  else
  {
    // 对父节点进行合法化
    std::pair<float, float> par_coord = {par_buf.coordinates.first, par_buf.coordinates.second};
    startLegalization(tmp_clock, par_coord, tmp_legal);
    par_buf.coordinates.first = par_coord.first;
    par_buf.coordinates.second = par_coord.second;

    // 创建子网
    std::string subnetName = "net_buf" + std::to_string(par_buf.index);
    cts::ClockSubNet &subnet = tmp_clock.addSubNet(subnetName);
    cts::ClockInst tmp_newBuf("BUF" + std::to_string(par_buf.index),
                              cts::CLOCK_BUFFER, par_buf.coordinates.first,
                              par_buf.coordinates.second, tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
    subnet.addInst(tmp_newBuf);
    tmp_clock.addClockBuffer(tmp_newBuf);

    // 将父节点buffer添加到临时legalization中
    Rectangle par_rect{static_cast<int>(par_buf.coordinates.first),
                       static_cast<int>(par_buf.coordinates.second),
                       tmp_clock.getBuf_W(),
                       tmp_clock.getBuf_H(),
                       "BUF"};
    tmp_legal.addFixedRectangle(par_rect);

    if (h_opt)
    {
      // insert n buffer from tmp_newBuf to buf1
      auto pre_loc1 = par_buf;
      auto *pre_sub_net1 = &subnet;
      Point buf1;
      buf1.coordinates.first = par_buf.coordinates.first;
      buf1.coordinates.second = par_buf.coordinates.second - seq[lev] / 2;

      auto delta_x = (buf1.coordinates.first - pre_loc1.coordinates.first) / (n + 1);
      auto delta_y = (buf1.coordinates.second - pre_loc1.coordinates.second) / (n + 1);
      for (size_t i = 0; i < n; ++i)
      {
        Point temp_opt_buf_loc;
        temp_opt_buf_loc.coordinates.first = pre_loc1.coordinates.first + delta_x;
        temp_opt_buf_loc.coordinates.second = pre_loc1.coordinates.second + delta_y;
        temp_opt_buf_loc.index = tmp_clock.getNumBuffer() + 1;
        tmp_clock.setNumBuffer(temp_opt_buf_loc.index);
        // leagalize
        std::pair<float, float> temp_opt_coord = {temp_opt_buf_loc.coordinates.first, temp_opt_buf_loc.coordinates.second};
        startLegalization(tmp_clock, temp_opt_coord, tmp_legal);
        cts::ClockInst temp_opt_Buf("BUF" + std::to_string(temp_opt_buf_loc.index),
                                    cts::CLOCK_BUFFER, temp_opt_coord.first, temp_opt_coord.second,
                                    tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
        tmp_legal.addFixedRectangle(Rectangle{static_cast<int>(temp_opt_buf_loc.coordinates.first),
                                              static_cast<int>(temp_opt_buf_loc.coordinates.second),
                                              tmp_clock.getBuf_W(),
                                              tmp_clock.getBuf_H(),
                                              "BUF"});
        pre_sub_net1->addInst(temp_opt_Buf);

        // update pre*
        pre_loc1 = temp_opt_buf_loc;
        auto temp_sub_net_name = "net_buf" + std::to_string(temp_opt_buf_loc.index);
        auto &new_net = tmp_clock.addSubNet(temp_sub_net_name);
        pre_sub_net1 = &new_net;
        pre_sub_net1->addInst(temp_opt_Buf);
        tmp_clock.addClockBuffer(temp_opt_Buf);
      }
      // connect buf1
      buf1.index = tmp_clock.getNumBuffer() + 1;
      tmp_clock.setNumBuffer(buf1.index);
      // leagalize
      std::pair<float, float> buf1_coord = {buf1.coordinates.first, buf1.coordinates.second};
      startLegalization(tmp_clock, buf1_coord, tmp_legal);
      cts::ClockInst tmp_Buf1("BUF" + std::to_string(buf1.index),
                              cts::CLOCK_BUFFER, buf1.coordinates.first, buf1.coordinates.second,
                              tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      tmp_legal.addFixedRectangle(Rectangle{static_cast<int>(buf1.coordinates.first),
                                            static_cast<int>(buf1.coordinates.second),
                                            tmp_clock.getBuf_W(),
                                            tmp_clock.getBuf_H(),
                                            "BUF"});
      pre_sub_net1->addInst(tmp_Buf1);

      // insert n buffer from tmp_newBuf to buf2
      auto pre_loc2 = par_buf;
      auto *pre_sub_net2 = &subnet;
      Point buf2;
      buf2.coordinates.first = par_buf.coordinates.first;
      buf2.coordinates.second = par_buf.coordinates.second + seq[lev] / 2;

      delta_x = (buf2.coordinates.first - pre_loc2.coordinates.first) / (n + 1);
      delta_y = (buf2.coordinates.second - pre_loc2.coordinates.second) / (n + 1);
      for (size_t i = 0; i < n; ++i)
      {
        Point temp_opt_buf_loc;
        temp_opt_buf_loc.coordinates.first = pre_loc2.coordinates.first + delta_x;
        temp_opt_buf_loc.coordinates.second = pre_loc2.coordinates.second + delta_y;
        temp_opt_buf_loc.index = tmp_clock.getNumBuffer() + 1;
        tmp_clock.setNumBuffer(temp_opt_buf_loc.index);
        // leagalize
        std::pair<float, float> temp_opt_coord = {temp_opt_buf_loc.coordinates.first, temp_opt_buf_loc.coordinates.second};
        startLegalization(tmp_clock, temp_opt_coord, tmp_legal);
        cts::ClockInst temp_opt_Buf("BUF" + std::to_string(temp_opt_buf_loc.index),
                                    cts::CLOCK_BUFFER, temp_opt_coord.first, temp_opt_coord.second,
                                    tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
        tmp_legal.addFixedRectangle(Rectangle{static_cast<int>(temp_opt_buf_loc.coordinates.first),
                                              static_cast<int>(temp_opt_buf_loc.coordinates.second),
                                              tmp_clock.getBuf_W(),
                                              tmp_clock.getBuf_H(),
                                              "BUF"});
        pre_sub_net2->addInst(temp_opt_Buf);

        // update pre*
        pre_loc2 = temp_opt_buf_loc;
        auto temp_sub_net_name = "net_buf" + std::to_string(temp_opt_buf_loc.index);
        auto &new_net = tmp_clock.addSubNet(temp_sub_net_name);
        pre_sub_net2 = &new_net;
        pre_sub_net2->addInst(temp_opt_Buf);
        tmp_clock.addClockBuffer(temp_opt_Buf);
      }
      // connect buf2
      buf2.index = tmp_clock.getNumBuffer() + 1;
      tmp_clock.setNumBuffer(buf2.index);
      // leagalize
      std::pair<float, float> buf2_coord = {buf2.coordinates.first, buf2.coordinates.second};
      startLegalization(tmp_clock, buf2_coord, tmp_legal);
      cts::ClockInst tmp_Buf2("BUF" + std::to_string(buf2.index),
                              cts::CLOCK_BUFFER, buf2.coordinates.first, buf2.coordinates.second,
                              tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      tmp_legal.addFixedRectangle(Rectangle{static_cast<int>(buf2.coordinates.first),
                                            static_cast<int>(buf2.coordinates.second),
                                            tmp_clock.getBuf_W(),
                                            tmp_clock.getBuf_H(),
                                            "BUF"});
      pre_sub_net2->addInst(tmp_Buf2);

      // recursive H tree
      tryBuildTree(direction + 1, lev + 1, buf1, seq, htree_leafs, tmp_clock, tmp_legal);
      tryBuildTree(direction + 1, lev + 1, buf2, seq, htree_leafs, tmp_clock, tmp_legal);
    }
    else
    {
      // 计算坐标
      Point buf1, buf2;
      buf1.coordinates.first = par_buf.coordinates.first;
      buf1.coordinates.second = par_buf.coordinates.second + seq[lev] / 2;
      buf1.index = tmp_clock.getNumBuffer() + 1;
      tmp_clock.setNumBuffer(buf1.index);
      buf2.coordinates.first = par_buf.coordinates.first;
      buf2.coordinates.second = par_buf.coordinates.second - seq[lev] / 2;
      buf2.index = tmp_clock.getNumBuffer() + 1;
      tmp_clock.setNumBuffer(buf2.index);

      // 对两个子节点进行合法化
      std::pair<float, float> buf1_coord = {buf1.coordinates.first, buf1.coordinates.second};
      std::pair<float, float> buf2_coord = {buf2.coordinates.first, buf2.coordinates.second};
      startLegalization(tmp_clock, buf1_coord, tmp_legal);
      startLegalization(tmp_clock, buf2_coord, tmp_legal);
      buf1.coordinates.first = buf1_coord.first;
      buf1.coordinates.second = buf1_coord.second;
      buf2.coordinates.first = buf2_coord.first;
      buf2.coordinates.second = buf2_coord.second;

      // 添加到子网中
      cts::ClockInst tmp_Buf1("BUF" + std::to_string(buf1.index),
                              cts::CLOCK_BUFFER, buf1.coordinates.first, buf1.coordinates.second,
                              tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      subnet.addInst(tmp_Buf1);

      cts::ClockInst tmp_Buf2("BUF" + std::to_string(buf2.index),
                              cts::CLOCK_BUFFER, buf2.coordinates.first, buf2.coordinates.second,
                              tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      subnet.addInst(tmp_Buf2);

      // 将子节点buffer添加到临时legalization中
      Rectangle buf1_rect{static_cast<int>(buf1.coordinates.first),
                          static_cast<int>(buf1.coordinates.second),
                          tmp_clock.getBuf_W(),
                          tmp_clock.getBuf_H(),
                          "BUF"};
      tmp_legal.addFixedRectangle(buf1_rect);

      Rectangle buf2_rect{static_cast<int>(buf2.coordinates.first),
                          static_cast<int>(buf2.coordinates.second),
                          tmp_clock.getBuf_W(),
                          tmp_clock.getBuf_H(),
                          "BUF"};
      tmp_legal.addFixedRectangle(buf2_rect);

      // 递归构建H树
      tryBuildTree(direction + 1, lev + 1, buf1, seq, htree_leafs, tmp_clock, tmp_legal);
      tryBuildTree(direction + 1, lev + 1, buf2, seq, htree_leafs, tmp_clock, tmp_legal);
    }
  }
}

void HTree::findMinLatency(int num_level,
                           int direction,
                           std::vector<int> &seq,
                           std::vector<std::vector<int>> w_len,
                           std::vector<std::vector<int>> h_len)
{
  auto temp_direction = direction;
  for (size_t lev = 0; lev < num_level; ++lev)
  {
    temp_direction++;
    // 横向
    if (temp_direction % 2 == 1)
    {
      seq.push_back(w_len[lev / 2][0]);
    }
    // 纵向
    else
    {
      seq.push_back(h_len[lev / 2][0]);
    }
  }

  // 建时钟树
  cts::Clock tmp_clock = clock_;
  Legalization tmp_legal = legal_;

  Point root_node;
  root_node.coordinates.first = root_node_x_;
  root_node.coordinates.second = root_node_y_;

  // 对根节点进行合法化，使用临时的legalization对象
  std::pair<float, float> root_coord = {root_node.coordinates.first, root_node.coordinates.second};
  startLegalization(tmp_clock, root_coord, tmp_legal);
  root_node.coordinates.first = root_coord.first;
  root_node.coordinates.second = root_coord.second;

  // 将根节点添加到临时legalization中
  Rectangle root_rect{static_cast<int>(root_node.coordinates.first),
                      static_cast<int>(root_node.coordinates.second),
                      tmp_clock.getBuf_W(),
                      tmp_clock.getBuf_H(),
                      "BUF"};
  tmp_legal.addFixedRectangle(root_rect);

  root_node.index = tmp_clock.getNumBuffer() + 1;
  tmp_clock.setNumBuffer(root_node.index);
  std::vector<Point> htree_leafs; // H树的叶子节点

  tryBuildTree(direction, 0, root_node, seq, htree_leafs, tmp_clock, tmp_legal);

  // Kmeans++找初始中心点
  // KMeansPlusPlus kmeanspp;

  // unsigned int numClusters = std::pow(2, num_level);
  // kmeanspp.run(points_, numClusters, 100);
  // std::vector<std::pair<float, float>> tmp_centroids;
  // tmp_centroids = kmeanspp.getCentroids();

  // 使用MinCostFlow聚类
  MinCostFlow<Point> mcf;
  for (size_t i = 0; i < points_.size(); ++i)
  {
    mcf.add_node(points_[i].coordinates.first,
                 points_[i].coordinates.second, points_[i]);
  }
  // 使用H树叶节点的坐标作为中心点
  for (const auto &leaf : htree_leafs)
  {
    mcf.add_center(leaf.coordinates.first, leaf.coordinates.second);
  }
  // for (size_t i = 0;i < tmp_centroids.size();++i) {
  //   mcf.add_center(tmp_centroids[i].first, tmp_centroids[i].second);
  // }

  double mcf_num = static_cast<double>(points_.size()) / std::pow(2, num_level);
  int mcf_fanout = static_cast<int>(ceil(mcf_num));
  std::vector<std::vector<Point>> clusters = mcf.run(mcf_fanout);

  // 初始化H树叶节点，一开始都可用
  std::vector<std::pair<Point, int>> htree_free_leafs(htree_leafs.size());
  for (size_t i = 0; i < htree_free_leafs.size(); ++i)
  {
    htree_free_leafs[i].first = htree_leafs[i];
    htree_free_leafs[i].second = 1;
  }

  // 每个簇找最近的H树的叶子节点
  std::vector<std::pair<std::vector<Point>, Point>> cluster_leaf_mapping;
  for (size_t i = 0; i < clusters.size(); ++i)
  {
    // 计算簇中所有点的坐标平均值
    float sum_x = 0, sum_y = 0;
    for (const auto &point : clusters[i])
    {
      sum_x += point.coordinates.first;
      sum_y += point.coordinates.second;
    }
    // 计算中心点坐标
    std::pair<float, float> centroid = {
        sum_x / clusters[i].size(),
        sum_y / clusters[i].size()};

    double minDist = std::numeric_limits<double>::max();
    int minIdx = -1;

    for (size_t j = 0; j < htree_free_leafs.size(); ++j)
    {
      if (htree_free_leafs[j].second == 1)
      {
        double tmp_dist = calcDist(centroid, htree_free_leafs[j].first);
        if (tmp_dist < minDist)
        {
          minDist = tmp_dist;
          minIdx = j;
        }
      }
    }

    if (minIdx != -1)
    {
      htree_free_leafs[minIdx].second = 0; // 标记为已占用
      // 保存簇和叶子节点的映射关系
      cluster_leaf_mapping.push_back({clusters[i], htree_free_leafs[minIdx].first});
    }
    else
    {
      std::cerr << "Error: No available leaf nodes for cluster " << i << std::endl;
    }
  }

  std::vector<cts::ClockSubNet> nets;
  // 把每个簇接到H树的叶子节点上
  for (auto &mapping : cluster_leaf_mapping)
  {
    auto &cluster = mapping.first;          // 簇中的点
    auto &htree_leaf_node = mapping.second; // H树叶子节点

    // 为这个H树叶子节点创建子网
    std::string subnetName = "net_buf" + std::to_string(htree_leaf_node.index);
    cts::ClockSubNet &subnet = tmp_clock.addSubNet(subnetName);

    // 添加H树叶子节点buffer到子网
    cts::ClockInst htree_leaf_buf("BUF" + std::to_string(htree_leaf_node.index),
                                  cts::CLOCK_BUFFER, htree_leaf_node.coordinates.first,
                                  htree_leaf_node.coordinates.second, tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
    subnet.addInst(htree_leaf_buf);

    // 将簇中的所有点连接到这个H树叶子节点
    for (auto &point_buf : cluster)
    {
      cts::ClockInst sink_buf("BUF" + std::to_string(point_buf.index),
                              cts::CLOCK_BUFFER, point_buf.coordinates.first,
                              point_buf.coordinates.second, tmp_clock.getBuf_W(), tmp_clock.getBuf_H());
      subnet.addInst(sink_buf);
    }

    nets.push_back(subnet);
  }
  // 计算这棵树的average_latency
  double tmp_latency = 0.0;
  double total_latency = 0.0;
  int path_count = 0;

  std::string rootBufName = "BUF" + std::to_string(root_node.index);
  // 递归计算所有路径的延迟
  calculateAllPathLatencies(tmp_clock, rootBufName, 0, 1, 0, total_latency, path_count);
  // 计算平均延迟
  tmp_latency = total_latency / path_count;

  if (tmp_latency < min_latency_)
  {
    min_clock_ = tmp_clock;
    curr_legal_ = tmp_legal;
    min_latency_ = tmp_latency;
    root_node_ = root_node;
  }
}

double HTree::calcDist(std::pair<float, float> p1, Point p2)
{
  return std::abs(p1.first - p2.coordinates.first) + std::abs(p1.second - p2.coordinates.second);
}

void HTree::calculateAllPathLatencies(cts::Clock &clock,
                                      const std::string &bufName,
                                      double curr_latency,
                                      int buf_count,
                                      int curr_level,
                                      double &total_latency,
                                      int &path_count)
{

  double BUF_DELAY = clock.getBufferDelay();

  // 如果到达第5层，累加当前路径的延迟
  if (curr_level == floor(log2(points_.size())))
  {
    total_latency += curr_latency;
    path_count++;
    return;
  }

  // 找到以该buffer为driver的子网
  for (auto &subnet : clock.getSubNets())
  {
    std::deque<cts::ClockInst> instances = subnet.getClockInst();
    if (instances[0].getName() == bufName)
    {
      // 遍历该buffer驱动的所有sink
      for (int i = 1; i < instances.size(); i++)
      {
        // 计算当前net的rc
        double rc = 0.5 * clock.getNetUnitR() * clock.getNetUnitC() *
                    pow(clock.calcDist(instances[0], instances[i]), 2);

        // 计算net delay
        double net_delay = 0.69 * rc;

        // 计算当前buffer的延迟
        double curr_path_latency = curr_latency;

        // 递归计算每条路径的延迟
        calculateAllPathLatencies(clock,
                                  instances[i].getName(),
                                  curr_path_latency + net_delay + BUF_DELAY,
                                  buf_count + 1,
                                  curr_level + 1,
                                  total_latency,
                                  path_count);
      }
      break;
    }
  }
}

void HTree::startLegalization(cts::Clock &clock, std::pair<float, float> &means, Legalization &legalization)
{

  Rectangle tmp_rect{means.first, means.second,
                     clock.getBuf_W(), clock.getBuf_H(), "BUF"};
  Rectangle rect = legalization.legalizeRectangle(tmp_rect);
  legalization.addFixedRectangle(rect);
  means.first = rect.x;
  means.second = rect.y;
}
