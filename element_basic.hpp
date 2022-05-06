#pragma once
#include "types.hpp"
#include <algorithm>
#include <iostream>
#include <map>
#include <numeric>
#include <set>
using std::accumulate;
using std::cerr;
using std::conditional_t;
using std::initializer_list;
using std::map;
using std::set;
using std::tuple;

// Thanks to https://www.mm.bme.hu/
// I J K L M N O P Q R  S  T  U  V  W  X  Y  Z  A  B
// 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
// 4
// face 1 (J-I-K)
// face 2 (I-J-L)
// face 3 (J-K-L)
// face 4 (K-I-L)
// 1 0 2
// 0 1 3
// 1 2 3
// 2 0 3
// 8
// face 1 (J-I-L-K)
// face 2 (I-J-N-M)
// face 3 (J-K-O-N)
// face 4 (K-L-P-O)
// face 5 (L-I-M-P)
// face 6 (M-N-O-P)
// 1 0 3 2
// 0 1 5 4
// 1 2 6 5
// 2 3 7 6
// 3 0 4 7
// 4 5 6 7
// 10
// face 1 (J-M-I-O-K-N)
// face 2 (I-M-J-Q-L-P)
// face 3 (J-N-K-R-L-Q)
// face 4 (K-O-I-P-L-R)
// 1 4 0 6 2 5
// 0 4 1 8 3 7
// 1 5 2 9 3 8
// 2 6 0 7 3 9
// 20
// face 1 (J-Q-I-T-L-S-K-R)
// face 2 (I-Q-J-Z-N-U-M-Y)
// face 3 (J-R-K-A-O-V-N-Z)
// face 4 (K-S-L-B-P-W-O-A)
// face 5 (L-T-I-Y-M-X-P-B)
// face 6 (M-U-N-V-O-W-P-X)
// 1  8 0 11 3 10 2  9
// 0  8 1 17 5 12 4 16
// 1  9 2 18 6 13 5 17
// 2 10 3 19 7 14 6 18
// 3 11 0 16 4 15 7 19
// 4 12 5 13 6 14 7 15

static inline const map<u32, tuple<u32, vector<u32>>> &get_element_polygons() {
  static map<u32, tuple<u32, vector<u32>>> element_polygons;

  if (element_polygons.empty()) {
    initializer_list<u32> l8 = {1, 0, 3, 2, //
                                0, 1, 5, 4, //
                                1, 2, 6, 5, //
                                2, 3, 7, 6, //
                                3, 0, 4, 7, //
                                4, 5, 6, 7};
    element_polygons.emplace(8, make_tuple(4, vector<u32>{l8}));
    initializer_list<u32> l20 = {1, 8,  0, 11, 3, 10, 2, 9,  //
                                 0, 8,  1, 17, 5, 12, 4, 16, //
                                 1, 9,  2, 18, 6, 13, 5, 17, //
                                 2, 10, 3, 19, 7, 14, 6, 18, //
                                 3, 11, 0, 16, 4, 15, 7, 19, //
                                 4, 12, 5, 13, 6, 14, 7, 15};
    element_polygons.emplace(20, make_tuple(8, vector<u32>{l20}));
    initializer_list<u32> l10 = {1, 4, 0, 6, 2, 5, //
                                 0, 4, 1, 8, 3, 7, //
                                 1, 5, 2, 9, 3, 8, //
                                 2, 6, 0, 7, 3, 9};
    element_polygons.emplace(10, make_tuple(6, vector<u32>{l10}));
  }
  return element_polygons;
}

static inline const map<u32, vector<u32>> &get_element_frames() {
  static map<u32, vector<u32>> element_frames;
  if (element_frames.empty()) {
    element_frames.emplace(8, vector<u32>{0, 1, 1, 2, 2, 3, 3, 0, //
                                          0, 4, 1, 5, 2, 6, 3, 7, //
                                          4, 5, 5, 6, 6, 7, 7, 4});
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
    element_frames.emplace(8, vector<u32>{0, 1, 2, 0, 2, 3, //
                                          4, 5, 6, 4, 6, 7, //
                                          0, 1, 5, 0, 5, 4, //
                                          2, 3, 7, 2, 7, 6, //
                                          1, 2, 6, 1, 6, 5, //
                                          3, 0, 4, 3, 4, 7});
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

/// @brief
/// @param sizes
/// @param indices
/// @param number_map
/// @return Polygon vertex counts, polygon vertex indices, element polygons
/// counts, node numbers, element numbers, surface numbers
template <bool withElementSize>
static inline tuple<vector<u32>, vector<u32>,
                    conditional_t<withElementSize, vector<u8>, tuple<>>,
                    vector<u32>, vector<u32>, vector<u8>>
get_polygon(const vector<u32> &sizes, const vector<u32> &indices,
            const vector<u32> &node_number_map,
            const vector<u32> &element_number_map) {
  const bool withNodeNumber = !node_number_map.empty();
  const bool withElementNumber = !element_number_map.empty();
  auto count = accumulate(sizes.begin(), sizes.end(), 0);
  assert(count == indices.size());
  tuple<vector<u32>, vector<u32>,
        conditional_t<withElementSize, vector<u8>, tuple<>>, vector<u32>,
        vector<u32>, vector<u8>>
      res;
  auto &[polygon_sizes, polygon_indices, element_sizes, node_numbers,
         element_numbers, surface_numbers] = res;
  auto p = indices.begin();
  const auto e = indices.end();
  auto &map = get_element_polygons();

  u32 i_element = 0;
  for (auto sz : sizes) {
    assert(p < e);
    auto &[polygon_size, polygon_vertex_indices] = map.at(sz);
    for (size_t i = 0; i < polygon_vertex_indices.size() / polygon_size; ++i) {
      polygon_sizes.push_back(polygon_size);
      if (withElementNumber)
        element_numbers.push_back(element_number_map[i_element]);

      surface_numbers.push_back(i);
    }
    for (auto polygon_vertex_index : polygon_vertex_indices) {
      assert(polygon_vertex_index < sz);
      polygon_indices.push_back(*(p + polygon_vertex_index));
      if (withNodeNumber)
        node_numbers.push_back(node_number_map[*(p + polygon_vertex_index)]);
    }
    if constexpr (withElementSize) {
      assert(polygon_vertex_indices.size() <
             polygon_size * (numeric_limits<u8>::max() + 1U));
      element_sizes.push_back(u8(polygon_vertex_indices.size() / polygon_size));
    }
    p += sz;
    ++i_element;
  }
  assert(p == e);
  return res;
}

static inline vector<u32> from_polygons(const vector<u32> &polygon_sizes,
                                        const vector<u32> &polygon_indices,
                                        const bool wireframe) {
  auto sum = accumulate(polygon_sizes.begin(), polygon_sizes.end(), 0);
  assert(sum == polygon_indices.size());
  auto p = polygon_indices.begin();
  auto e = polygon_indices.end();
  vector<u32> res;

  for (auto size : polygon_sizes) {
    assert(p != e);

    if (wireframe) {
      for (size_t i = 1; i < size; ++i) {
        res.push_back(*(p + i - 1));
        res.push_back(*(p + i));
      }
      res.push_back(*p);
      res.push_back(*(p + size - 1));
    } else {
      for (size_t i = 2; i < size; ++i) {
        res.push_back(*p);
        res.push_back(*(p + i - 1));
        res.push_back(*(p + i));
      }
    }
    p += size;
  }
  assert(p == e);
  return res;
}

/// @brief Return real indices of nodes.
/// @param nodes: Nodes coordinates.
/// @param vertex_count: Count of nodes in elements
/// @param elements: Indices of nodes in element.
/// @param wireframe
/// @return Real indices of nodes in primitives.
static inline vector<u32> get_primitives(const vector<float> &nodes,
                                         const vector<u32> &vertex_count,
                                         const vector<u32> &elements,
                                         const bool wireframe) {
  vector<u32> indices;

  u32 sum = accumulate(vertex_count.begin(), vertex_count.end(), 0U);
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
    cerr << "Element with ";
    for (auto vc : unrecognized)
      cerr << vc << ' ';
    cerr << "nodes are not recoginzed." << endl;
  }

  return indices;
}

static inline tuple<tuple<vector<float>>, vector<u32>>
from_elements(const vector<float> &nodes, const vector<u32> &vertex_count,
              const vector<u32> &elements, const bool wireframe,
              const float ratio) {
  tuple<tuple<vector<float>>, vector<u32>> res;
  auto &[vertices, indices] = res;
  auto &[position] = vertices;
  position = nodes;
  for (auto &p : position)
    p *= ratio;

  indices = get_primitives(nodes, vertex_count, elements, wireframe);
  return res;
}

template <typename Ty> static inline bool set_contains(const set<Ty> &s, Ty v) {
  return s.find(v) != s.end();
}
template <typename Ty>
static inline bool set_contains_all(const set<Ty> &s, initializer_list<Ty> l) {
  for (auto &v : l)
    if (!set_contains(s, v))
      return false;
  return true;
}
template <typename Ty>
static inline bool set_contains_all(const set<Ty> &s, const vector<Ty> &l) {
  for (auto &v : l)
    if (!set_contains(s, v))
      return false;
  return true;
}

static inline tuple<vector<float>, vector<float>, vector<float>>
coordinates(const vector<float> &vec, const vector<u32> &indices) {
  tuple<vector<float>, vector<float>, vector<float>> res;
  auto &[x, y, z] = res;
  x.reserve(indices.size());
  y.reserve(indices.size());
  z.reserve(indices.size());
  for (auto i : indices) {
    x.push_back(vec[3 * i]);
    y.push_back(vec[3 * i + 1]);
    z.push_back(vec[3 * i + 2]);
  }
  return res;
}

/// @brief Return real indices of nodes.
/// @param nodes: Nodes coordinates.
/// @param vertex_count: Count of nodes in elements
/// @param elements: Indices of nodes in element.
/// @param wireframe
/// @return Real indices of nodes in primitives.
static inline vector<float>
polygon_as_centroid(const vector<float> &nodes,
                    const vector<u32> &polygon_sizes,
                    const vector<u32> &polygon_indices) {
  vector<float> res;

  res.reserve(3 * polygon_sizes.size());

  assert(accumulate(polygon_sizes.begin(), polygon_sizes.end(), 0) ==
         polygon_indices.size());

  auto p = polygon_indices.begin();
  auto e = polygon_indices.end();
  for (const auto &sz : polygon_sizes) {
    assert(p < e);
    float x = 0;
    float y = 0;
    float z = 0;
    for (size_t i = 0; i < sz; ++i) {
      auto I = *p++;
      x += nodes[3 * I];
      y += nodes[3 * I + 1];
      z += nodes[3 * I + 2];
    }
    res.push_back(x / sz);
    res.push_back(y / sz);
    res.push_back(z / sz);
  }

  return res;
}