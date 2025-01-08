
#pragma once

#include <iostream>
#include <cassert>
#include <cstdint>
#include <deque>
#include <functional>
#include <unordered_map>
#include <vector>
#include <string>
#include <set>
#include <cmath>

// 定义一个结构体，用来存储每个数据点的索引和坐标
struct Point {
  Point() : index(0), coordinates{ 0.0, 0.0 } {}
  Point(int idx, std::pair<float, float> coords)
    : index(idx), coordinates(coords) {
  }

  int index;                           // 数据点的原始索引
  std::pair<float, float> coordinates; // 数据点的坐标

};

namespace cts {

  enum InstType : uint8_t
  {
    CLOCK_BUFFER,
    CLOCK_FF,
    CLOCK_CLK
  };

  class ClockInst
  {
  public:
    ClockInst(const std::string& name) :name_(name) {}

    ClockInst(const std::string& name,
      InstType type,
      int x,
      int y,
      int w,
      int h)
      : name_(name),
      type_(type),
      location_x(x),
      location_y(y),
      inst_w(w),
      inst_h(h)
    {
    }

    std::string getName() const { return name_; }
    int getX() const { return location_x; }
    int getY() const { return location_y; }
    int getW() const { return inst_w; }
    int getH() const { return inst_h; }
    bool isClockBuffer() const { return type_ == CLOCK_BUFFER; }


  private:
    std::string name_;
    InstType type_;
    int location_x;
    int location_y;
    int inst_w;
    int inst_h;
  };

  class ClockSubNet
  {
  private:
    std::string name_;
    std::deque<ClockInst> instances_;
    std::unordered_map<ClockInst*, unsigned> mapInstToIdx_;
    bool leafLevel_ = false;

  public:
    explicit ClockSubNet(const std::string& name) : name_(name) {}

    // setter
    void setLeafLevel(bool isLeaf) { leafLevel_ = isLeaf; }

    // adder
    void addInst(ClockInst& inst)
    {
      instances_.push_back(inst);
      mapInstToIdx_[&inst] = instances_.size() - 1;
    }

    // getter
    std::string getName() const { return name_; }
    int getNumFFs() const { return instances_.size() - 1; }
    std::deque<ClockInst> getClockInst() { return instances_; }

    bool isLeafLevel() const { return leafLevel_; }
    unsigned findIndex(ClockInst* inst) const { return mapInstToIdx_.at(inst); }

  };

  class Clock
  {
  private:
    std::string netName_;
    int clockPinX_;                 // clk的横坐标
    int clockPinY_;                 // clk的纵坐标
    int units_distance_microns_;    // 缩放比例
    int die_area_x_, die_area_y_;   // 版图的宽和高
    int buf_w_, buf_h_;             // BUFFER的宽和高
    int ff_w_, ff_h_;               // FF的宽和高
    int max_fanout_;                // 最大扇出
    int numBuffer_ = 0;             // Buffer的数量
    double net_unit_r_;             // 单位长度电阻值r
    double net_unit_c_;             // 电容值c
    double max_net_rc_;
    double buffer_delay_;

    std::deque<ClockInst> FFs_;
    std::deque<ClockInst> clockBuffers_;
    std::deque<ClockSubNet> subNets_;
    std::unordered_map<std::string, ClockInst*> mapNameToInst_;
    std::unordered_map<std::string, ClockSubNet*> mapNameToSubNet_;

    unsigned numLevels_ = 0;

  public:
    Clock() {
      netName_ = "no name";
      clockPinX_ = -1;
      clockPinY_ = -1;
    }

    Clock(const std::string& netName,
      int clockPinX,
      int clockPinY)
      :netName_(netName),
      clockPinX_(clockPinX),
      clockPinY_(clockPinY) {
    }


    // getter
    int getUnitsDistanceMicrons() { return units_distance_microns_; }
    int getDieArea_X() { return die_area_x_; }
    int getDieArea_Y() { return die_area_y_; }
    int getClockPinX() { return clockPinX_; }
    int getClockPinY() { return clockPinY_; }
    int getBuf_W() { return buf_w_; }
    int getBuf_H() { return buf_h_; }
    int getFF_W() { return ff_w_; }
    int getFF_H() { return ff_h_; }
    int getMaxFanout() { return max_fanout_; }
    int getNumBuffer() { return numBuffer_; }
    double getNetUnitR() { return net_unit_r_; }
    double getNetUnitC() { return net_unit_c_; }
    double getMaxNetRC() { return max_net_rc_; }
    double getBufferDelay() { return buffer_delay_; }
    std::string getName() const { return netName_; }
    unsigned getNumFFs() const { return FFs_.size(); }
    unsigned getMaxLevel() const { return numLevels_; }
    std::deque<ClockInst>& getFFs() { return FFs_; }
    std::deque<ClockInst>& getBuffers() { return clockBuffers_; }
    std::deque<ClockSubNet>& getSubNets() { return subNets_; }


    // setter
    void setUnitsDistanceMicrons(int units_distance_microns) {
      units_distance_microns_ = units_distance_microns;
    }
    void setDieArea(int die_area_x, int die_area_y) {
      die_area_x_ = die_area_x;
      die_area_y_ = die_area_y;
    }
    void setBufferInfo(int buf_w, int buf_h) {
      buf_w_ = buf_w;
      buf_h_ = buf_h;
    }
    void setFFInfo(int ff_w, int ff_h) {
      ff_w_ = ff_w;
      ff_h_ = ff_h;
    }
    void setMaxFanout(int max_fanout) {
      max_fanout_ = max_fanout;
    }
    void setMaxLevel(unsigned level) { numLevels_ = level; }
    void setClockPin(int clockPinX, int clockPinY) {
      clockPinX_ = clockPinX;
      clockPinY_ = clockPinY;
    }
    void setNumBuffer(int numBuffer) {
      numBuffer_ = numBuffer;
    }
    void setNetUnitR(double net_unit_r) {
      net_unit_r_ = net_unit_r;
    }
    void setNetUnitC(double net_unit_c) {
      net_unit_c_ = net_unit_c;
    }
    void setMaxNetRC(double max_net_rc) {
      max_net_rc_ = max_net_rc;
    }
    void setBufferDelay(double buffer_delay){
      buffer_delay_ = buffer_delay;
    }

    // adder
    ClockSubNet& addSubNet(const std::string& name) {
      subNets_.emplace_front(name);   // 头插法
      mapNameToSubNet_[name] = &subNets_.front();
      return subNets_.front();
    }
    void addFF(ClockInst& inst) {
      FFs_.emplace_back(inst);
    }
    void addClockBuffer(ClockInst& inst) {
      clockBuffers_.emplace_back(inst);
    }

    // other function
    ClockInst* findClockInstByName(const std::string& name) {
      if (mapNameToInst_.find(name) == mapNameToInst_.end()) {
        return nullptr;
      }
      return mapNameToInst_.at(name);
    }

    ClockSubNet* findClockSubNetByName(const std::string& name) {
      if (mapNameToSubNet_.find(name) == mapNameToSubNet_.end()) {
        return nullptr;
      }
      return mapNameToSubNet_.at(name);
    }

    void printFF() {
      for (ClockInst& ff : FFs_) {
        std::cout << ff.getName() << " " << ff.getX()
          << " " << ff.getY() << std::endl;
      }
    }

    void printResult() {
      for (auto& buffer : clockBuffers_) {
        std::cout << "- " << buffer.getName() << " BUF ( "
          << buffer.getX() << " " << buffer.getY() << " ) ;" << std::endl;
      }
      std::cout << "NETS " << subNets_.size() << " ;" << std::endl;
      for (auto& subnet : subNets_) {
        std::deque<ClockInst> instances = subnet.getClockInst();
        std::cout << "- " << subnet.getName() << " ( "
          << instances[0].getName() << " ) (";
        for (int i = 1;i < instances.size();i++) {
          std::cout << " " << instances[i].getName();
        }
        std::cout << " ) ;" << std::endl;
      }
      std::cout << "END NETS" << std::endl;
    }

    int calcDist(ClockInst inst1, ClockInst inst2) {
      // 计算曼哈顿距离
      return std::abs((inst1.getX() + inst1.getW() / 2) - (inst2.getX() + inst2.getW() / 2))
        + std::abs((inst1.getY() + inst1.getH() / 2) - (inst2.getY() + inst2.getH() / 2));
    }

    void calcNetRC(const std::string& net_name) {
      ClockSubNet subnet = *findClockSubNetByName(net_name);
      std::deque<ClockInst> clockInst = subnet.getClockInst();
      ClockInst buffer = clockInst[0];
      double half_rc = 0.5 * net_unit_r_ * net_unit_c_;
      double net_rc = 0.0;

      std::cout << net_name << " information: " << std::endl;
      for (int i = 1; i < clockInst.size(); ++i) {
        double rc = half_rc * pow(calcDist(buffer, clockInst[i]), 2);
        std::cout << buffer.getName() << " to "
          << clockInst[i].getName() << " rc: " << rc << std::endl;
        net_rc += rc;
      }
      std::cout << net_name << " total net_rc is: " << net_rc << std::endl;
    }

    float calcNetRC(cts::ClockSubNet& subnet) {
      std::deque<cts::ClockInst> clockInst = subnet.getClockInst();
      cts::ClockInst buffer = clockInst[0];
      double half_rc = 0.5 * net_unit_r_ * net_unit_c_;
      double net_rc = 0.0;
      for (int i = 1; i < clockInst.size(); ++i) {
        double rc = half_rc * pow(calcDist(buffer, clockInst[i]), 2);
        net_rc += rc;
      }
      return net_rc;
    }



  };

}  // namespace cts
