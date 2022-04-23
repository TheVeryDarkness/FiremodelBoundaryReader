#pragma once

#include "element_basic.hpp"
#include "fds_basic.hpp"
#include <cassert>

using std::make_tuple;

static inline tuple<tuple<vector<float>>, vector<u32>>
from_patches_and_elements(const vector<float> &nodes,
                          const vector<u32> &vertex_count,
                          const vector<u32> &elements,
                          const vector<patch_info> &patches, bool wireframe,
                          u32 sum) {
  auto &&[_, i1] = from_patches<false>(patches, wireframe);
  auto &[p1, __] = _;
  auto &&[___, i2] =
      from_elements(nodes, vertex_count, elements, sum, wireframe);
  auto &[p2] = ___;

  p2.reserve(p1.size() + p2.size());
  assert(p2.size() / 2 < numeric_limits<u32>::max());
  u32 base = u32(p2.size() / 3);

  for (auto p = p1.cbegin(); p != p1.cend();) {
    auto i = *p++;
    auto j = *p++;
    auto k = *p++;
    p2.push_back(float(i * mesh.cell_size + mesh.x0));
    p2.push_back(float(j * mesh.cell_size + mesh.y0));
    p2.push_back(float(k * mesh.cell_size + mesh.z0));
  }
  for (auto i : i1)
    i2.push_back(base + i);
  return make_tuple(make_tuple(std::move(p2)), std::move(i2));
}

static inline tuple<tuple<vector<float>>, vector<GLuint>>
from_matrix_data(const vector<float> &data, u32 n) {
  tuple<tuple<vector<float>>, vector<GLuint>> res;
  auto &[vertices, indices] = res;
  auto &[position] = vertices;

  position.reserve(data.size() * 3);
  indices.reserve(wireframe ? indices.size() * 4 : indices.size() * 6);

  u32 index = 0;
  for (const auto &p : data) {
    auto i = index / n, j = index % n;
    if (j != 0)
      if (wireframe) {
        // 1
        // 2
        indices.push_back(i * n + j);
        indices.push_back(i * n + j - 1);
      }
    if (i != 0) {
      if (wireframe) {
        // 2 1
        indices.push_back(i * n + j);
        indices.push_back((i - 1) * n + j);
      } else {
        //   2
        // 3 1
        indices.push_back(i * n + j);
        indices.push_back(i * n + j + 1);
        indices.push_back((i - 1) * n + j);
        // 2 1
        // 3
        indices.push_back(i * n + j + 1);
        indices.push_back((i - 1) * n + j + 1);
        indices.push_back((i - 1) * n + j);
      }
    }

    position.push_back((float)i);
    position.push_back((float)j);
    position.push_back((float)p);

    ++index;
  }

  return res;
}