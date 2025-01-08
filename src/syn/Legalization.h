#ifndef LEGALIZATION_H
#define LEGALIZATION_H

#include <iostream>
#include <vector>
#include <queue>
#include <string>
#include <unordered_set>
#include <cmath>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>

// Namespace aliases for convenience
namespace bg = boost::geometry;
namespace bgi = boost::geometry::index;

// Rectangle structure
struct Rectangle
{
  int x, y, width, height;
  std::string type;

  // Method to check if two rectangles overlap
  bool overlaps(const Rectangle &other) const
  {
    return !(x + width <= other.x || other.x + other.width <= x ||
             y + height <= other.y || other.y + other.height <= y);
  }

  // Method to check if a rectangle is fully contained within this rectangle
  bool contains(const Rectangle &other) const
  {
    return (x <= other.x && x + width >= other.x + other.width &&
            y <= other.y && y + height >= other.y + other.height);
  }

  // Equality operator for unordered_set
  bool operator==(const Rectangle &other) const
  {
    return x == other.x &&
           y == other.y &&
           width == other.width &&
           height == other.height &&
           type == other.type;
  }
};

// Hash function for Rectangle
namespace std
{
  template <>
  struct hash<Rectangle>
  {
    std::size_t operator()(const Rectangle &rect) const
    {
      std::size_t h1 = std::hash<int>{}(rect.x);
      std::size_t h2 = std::hash<int>{}(rect.y);
      std::size_t h3 = std::hash<int>{}(rect.width);
      std::size_t h4 = std::hash<int>{}(rect.height);
      // std::size_t h5 = std::hash<std::string>{}(rect.type);
      return h1 ^ (h2 << 1) ^ (h3 << 2) ^ (h4 << 3);
    }
  };
} // namespace std

// Define Boost Geometry types
typedef bg::model::point<int, 2, bg::cs::cartesian> BoostPoint;
typedef bg::model::box<BoostPoint> BoostBox;
typedef std::pair<BoostBox, Rectangle> RTreeValue;

// Comparator for priority queue based on total delta distance
struct CompareRectangles
{
  bool operator()(const std::pair<Rectangle, int> &a, const std::pair<Rectangle, int> &b) const
  {
    return a.second > b.second; // Min-heap based on delta distance
  }
};

// Legalization class
class Legalization
{
public:
  Legalization() : dieWidth_(0), dieHeight_(0) {}

  Legalization(int dieWidth, int dieHeight)
      : dieWidth_(dieWidth), dieHeight_(dieHeight)
  {
    setDiearea(dieWidth, dieHeight);
  }

  // Copy assignment operator
  void operator=(const Legalization &other)
  {
    dieWidth_ = other.dieWidth_;
    dieHeight_ = other.dieHeight_;
    fixedRectanglesRTree_ = other.fixedRectanglesRTree_;
  }

  // Set die area dimensions
  void setDiearea(int dieWidth, int dieHeight)
  {
    dieWidth_ = dieWidth;
    dieHeight_ = dieHeight;
  }

  // Add a fixed rectangle to the R-tree
  void addFixedRectangle(const Rectangle &rect)
  {
    BoostBox box = rectToBoostBox(rect);
    fixedRectanglesRTree_.insert(std::make_pair(box, rect));
  }

  // Legalize a rectangle by finding a non-overlapping position
  Rectangle legalizeRectangle(const Rectangle &current);

private:
  int dieWidth_, dieHeight_;

  // R-tree for fixed rectangles
  bgi::rtree<RTreeValue, bgi::quadratic<16>> fixedRectanglesRTree_;

  // Convert Rectangle to BoostBox
  BoostBox rectToBoostBox(const Rectangle &rect) const
  {
    BoostPoint lower_left(rect.x, rect.y);
    BoostPoint upper_right(rect.x + rect.width, rect.y + rect.height);
    return BoostBox(lower_left, upper_right);
  }

  // Query overlapping rectangles using R-tree
  std::vector<Rectangle> queryOverlappingRectangles(const Rectangle &rect) const
  {
    BoostBox queryBox = rectToBoostBox(rect);
    std::vector<RTreeValue> result_s;
    fixedRectanglesRTree_.query(bgi::intersects(queryBox), std::back_inserter(result_s));

    std::vector<Rectangle> overlappingRects;
    for (const auto &value : result_s)
    {
      if (rect.overlaps(value.second))
      {
        overlappingRects.push_back(value.second);
      }
    }
    return overlappingRects;
  }

  // Check if rectangle is within die area
  bool isWithinDiearea(const Rectangle &rect) const
  {
    return rect.x >= 0 &&
           rect.y >= 0 &&
           (rect.x + rect.width) <= dieWidth_ &&
           (rect.y + rect.height) <= dieHeight_;
  }
};
#endif