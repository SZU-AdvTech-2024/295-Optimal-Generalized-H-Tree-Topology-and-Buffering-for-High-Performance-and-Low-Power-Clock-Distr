#pragma once
#include <algorithm>
#include <cctype>
#include <string>
#include <string_view>

namespace issta {
class Str {
 public:
  static bool contain(std::string& str, std::string sub_str) {
    return str.find(sub_str) != std::string::npos;
  }

  static bool containsSubstringIgnoreCase(std::string_view str,
                                          std::string_view substring) {
    auto it = std::search(str.begin(), str.end(), substring.begin(),
                          substring.end(), [](char ch1, char ch2) {
                            return std::tolower(ch1) == std::tolower(ch2);
                          });

    return (it != str.end());
  }
};
}  // namespace issta