#pragma once
#include "types.hpp"
#include <map>
using std::initializer_list;
using std::map;

static inline const map<u32, vector<u32>> &get_element_frames() {
  static map<u32, vector<u32>> element_frames;
  if (element_frames.empty()) {
    element_frames.emplace(8, vector<u32>{0, 1, 1, 2, 2, 3, 3, 0, 0, 4, 1, 5,
                                          2, 6, 3, 7, 4, 5, 5, 6, 6, 7, 7, 4});
    element_frames.emplace(
        20,
        vector<u32>{0, 8,  8,  1, 1, 9,  9,  2, 2, 10, 10, 3, 3, 11, 11, 0,
                    0, 16, 16, 4, 1, 17, 17, 5, 2, 18, 18, 6, 3, 19, 19, 7,
                    4, 12, 12, 5, 5, 13, 13, 6, 6, 14, 14, 7, 7, 15, 15, 4});
    element_frames.emplace(10, vector<u32>{0, 4, 4, 1, 1, 5, 5, 2, 2, 6, 6, 0,
                                           0, 7, 7, 3, 1, 8, 8, 3, 2, 9, 9, 3});
  }
  return element_frames;
}

static inline const map<u32, vector<u32>> &get_element_triangles() {
  static map<u32, vector<u32>> element_frames;
  if (element_frames.empty()) {
    // 0 1 2 3
    // 4 5 6 7
    // 0 1 5 4
    // 2 3 7 6
    // 1 2 6 5
    // 3 0 4 7
    element_frames.emplace(8, vector<u32>{0, 1, 2, 0, 2, 3, 4, 5, 6, 4, 6, 7,
                                          0, 1, 5, 0, 5, 4, 2, 3, 7, 2, 7, 6,
                                          1, 2, 6, 1, 6, 5, 3, 0, 4, 3, 4, 7});
    //  0  8  1  9  2 10  3 11
    //  4 12  5 13  6 14  7 15
    //  0  8  1 17  5 12  4 16
    //  1  9  2 18  6 13  5 17
    //  2 10  3 19  7 14  6 18
    //  3 11  0 16  4 15  7 19
    initializer_list<u32> l20 = {
        0, 8,  11, 8,  1, 9,  2, 10, 9,  10, 3, 11, 8,  9,  11, 9,  10, 11, //
        4, 12, 15, 12, 5, 13, 6, 14, 13, 14, 7, 15, 12, 13, 15, 13, 14, 15, //
        0, 8,  16, 8,  1, 17, 5, 12, 17, 12, 4, 16, 8,  17, 16, 17, 12, 16, //
        1, 9,  17, 9,  2, 18, 6, 13, 18, 13, 5, 17, 9,  18, 17, 18, 13, 17, //
        2, 10, 18, 10, 3, 19, 7, 14, 19, 14, 6, 18, 10, 19, 18, 19, 14, 18, //
        3, 11, 19, 11, 0, 16, 4, 15, 16, 15, 7, 19, 11, 16, 19, 16, 15, 19  //
    };
    element_frames.emplace(20, vector<u32>{l20});
    // 0 4 1 8 3 7
    // 1 5 2 9 3 8
    // 2 6 0 7 3 9
    // 0 6 2 5 1 4
    element_frames.emplace(10, vector<u32>{0, 4, 7, 4, 1, 8, 7, 8, 3, 4, 8, 7,
                                           1, 5, 8, 5, 2, 9, 8, 9, 3, 5, 9, 8,
                                           2, 6, 9, 6, 0, 7, 9, 7, 3, 6, 7, 9,
                                           0, 6, 4, 6, 2, 5, 4, 5, 1, 6, 5, 4});
  }
  return element_frames;
}

static inline tuple<tuple<vector<float>>, vector<GLuint>>
from_elements(const vector<float> &nodes, const vector<u32> &vertex_count,
              const vector<u32> &elements, bool wireframe, u32 sum,
              float ratio = 1) {
  tuple<tuple<vector<float>>, vector<GLuint>> res;
  auto &[vertices, indices] = res;
  auto &[position] = vertices;
  position = nodes;
  for (auto &p : position)
    p *= ratio;

  indices.reserve(3 * sum);

  u32 index = 0;

  set<u32> unrecognized;

  size_t out = 0;
  if (wireframe) {
    for (const auto &vc : vertex_count) {
      auto &map = get_element_frames();
      auto iter = map.find(vc);
      if (iter != map.cend()) {
        auto &[c, order] = *iter;
        for (auto i : order) {
          assert(i < c);
          auto _i = elements.at(index + i);
          if (_i >= sum)
            ++out;
          indices.push_back(_i);
        }
      } else {
        unrecognized.insert(vc);

        for (auto p = elements.cbegin() + index;
             p < elements.cbegin() + index + vc - 1; ++p) {
          indices.push_back(*p);
          indices.push_back(*(p + 1));
        }
      }
      index += vc;
    }
  } else {
    for (const auto &vc : vertex_count) {
      auto &map = get_element_triangles();
      auto iter = map.find(vc);
      if (iter != map.cend()) {
        auto &[c, order] = *iter;
        for (auto i : order) {
          assert(i < c);
          auto _i = elements.at(index + i);
          if (_i >= sum)
            ++out;
          indices.push_back(_i);
        }
      } else {
        unrecognized.insert(vc);

        for (auto p = elements.cbegin() + index;
             p < elements.cbegin() + index + vc - 2; ++p) {
          indices.push_back(*p);
          indices.push_back(*(p + 1));
          indices.push_back(*(p + 2));
        }
      }
      index += vc;
    }
  }

  if (!unrecognized.empty()) {
    cout << "Element with ";
    for (auto vc : unrecognized)
      cout << vc << ' ';
    cout << "nodes are not recoginzed." << endl;
  }

  if (out > 0)
    clog << out << " indices out of nodes count are detected." << endl;

  return res;
}
