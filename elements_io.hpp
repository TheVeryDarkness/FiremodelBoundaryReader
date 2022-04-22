#pragma once
#include "types.hpp"
#include <fstream>
#include <iostream>
#include <string>
#include <string_view>
#include <tuple>
#include <vector>

using std::cerr;
using std::char_traits;
using std::endl;
using std::getline;
using std::ifstream;
using std::isblank;
using std::isspace;
using std::istream;
using std::string;
using std::tuple;
using std::vector;
using std::ws;

/// @brief Read nodes data from stream associated with standard file.
/// @param in: Stream associated with standard file.
/// @return Nodes coordinates.
static inline vector<float> read_nodes(istream &in) {
  u32 count = 0;
  in >> count;
  string s;
  getline(in, s);
  s.clear();
  vector<float> res;
  for (u32 i = 0; i < count; ++i) {
    float x = 0, y = 0, z = 0;
    in >> x >> y >> z;
    res.push_back(x);
    res.push_back(y);
    res.push_back(z);
  }
  return res;
}

/// @brief Read elements data from stream associated with standard file.
/// @param in: Stream associated with standard file.
/// @return Indices of Element vertices, count of element vertices.
static inline tuple<vector<u32>, vector<u32>> read_elements(istream &in) {
  tuple<vector<u32>, vector<u32>> res;
  auto &[indices, sizes] = res;
  while (in.peek() != char_traits<char>::eof() && in) {
    u32 elements_count = 0, vertices_count = 0;
    in >> elements_count >> vertices_count;
    for (size_t i = 0; i < elements_count; ++i) {
      for (size_t j = 0; j < vertices_count; ++j) {
        u32 index = 0;
        in >> index;
        indices.push_back(index);
      }
      sizes.push_back(vertices_count);
    }
    in >> ws;
  }
  return res;
}
