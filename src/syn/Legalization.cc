#include "Legalization.h"

Rectangle Legalization::legalizeRectangle(const Rectangle &current)
{
  // If current position is already non-overlapping, return it
  auto overlapping = queryOverlappingRectangles(current);
  if (overlapping.empty())
  {
    return current;
  }

  // Priority queue for BFS, prioritized by total delta distance
  std::priority_queue<
      std::pair<Rectangle, int>,
      std::vector<std::pair<Rectangle, int>>,
      CompareRectangles>
      pq;

  // Initialize with the current rectangle and zero distance
  pq.emplace(current, 0);

  int origin_x = current.x;
  int origin_y = current.y;

  // Visited set to avoid revisiting the same positions
  std::unordered_set<Rectangle> visited;
  visited.insert(current);

  while (!pq.empty())
  {
    auto [rect, distance] = pq.top();
    pq.pop();

    // Query overlapping rectangles
    auto overlaps = queryOverlappingRectangles(rect);
    if (overlaps.empty())
    {
      // Found a non-overlapping position
      return rect;
    }

    // For each overlapping rectangle, calculate possible shifts
    for (const auto &overlapRect : overlaps)
    {
      // Calculate deltas in four directions
      // Shift Left: Move rect to the left of overlapRect
      int shiftLeft = overlapRect.x - rect.width;
      if (shiftLeft >= 0)
      {
        Rectangle newRect = rect;
        newRect.x = shiftLeft;
        if (isWithinDiearea(newRect) && visited.find(newRect) == visited.end())
        {
          int newDistance = std::abs(origin_x - newRect.x) + std::abs(origin_y - newRect.y);
          pq.emplace(newRect, newDistance);
          visited.insert(newRect);
        }
      }

      // Shift Right: Move rect to the right of overlapRect
      int shiftRight = overlapRect.x + overlapRect.width;
      if (shiftRight + rect.width <= dieWidth_)
      {
        Rectangle newRect = rect;
        newRect.x = shiftRight;
        if (isWithinDiearea(newRect) && visited.find(newRect) == visited.end())
        {
          int newDistance = std::abs(origin_x - newRect.x) + std::abs(origin_y - newRect.y);
          pq.emplace(newRect, newDistance);
          visited.insert(newRect);
        }
      }

      // Shift Down: Move rect below overlapRect
      int shiftDown = overlapRect.y - rect.height;
      if (shiftDown >= 0)
      {
        Rectangle newRect = rect;
        newRect.y = shiftDown;
        if (isWithinDiearea(newRect) && visited.find(newRect) == visited.end())
        {
          int newDistance = std::abs(origin_x - newRect.x) + std::abs(origin_y - newRect.y);
          pq.emplace(newRect, newDistance);
          visited.insert(newRect);
        }
      }

      // Shift Up: Move rect above overlapRect
      int shiftUp = overlapRect.y + overlapRect.height;
      if (shiftUp + rect.height <= dieHeight_)
      {
        Rectangle newRect = rect;
        newRect.y = shiftUp;
        if (isWithinDiearea(newRect) && visited.find(newRect) == visited.end())
        {
          int newDistance = std::abs(origin_x - newRect.x) + std::abs(origin_y - newRect.y);
          pq.emplace(newRect, newDistance);
          visited.insert(newRect);
        }
      }
    }
  }

  std::cout << "Error: No non-overlapping position found for rectangle" << std::endl;
  return current;
}
