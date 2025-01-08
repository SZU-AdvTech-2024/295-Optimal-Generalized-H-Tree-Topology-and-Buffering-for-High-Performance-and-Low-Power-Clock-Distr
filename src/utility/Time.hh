#pragma once

#include <chrono>
#include <string>

inline std::string GetCurrWallTime() {
  using namespace std::chrono;

  auto now = system_clock::now();
  auto time_now = system_clock::to_time_t(now);

  std::string curr_wall_time_str = ctime(&time_now);
  return curr_wall_time_str;
}