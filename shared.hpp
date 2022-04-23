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
from_matrix_data(const vector<float> &data, u32 n, bool wireframe) {
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

static inline set<u32> node_on_boundary(const vector<patch_info> &patches,
                                        const vector<float> &nodes) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  set<u32> res;

  u32 i = 0;

  float xmin = nodes[0], xmax = nodes[0], ymin = nodes[1], ymax = nodes[1],
        zmin = nodes[2], zmax = nodes[2];
  for (auto p = nodes.begin(), end = nodes.end(); p != end;) {
    float x = *p++;
    float y = *p++;
    float z = *p++;
    xmin = min(x, xmin);
    xmax = max(x, xmax);
    ymin = min(y, ymin);
    ymax = max(y, ymax);
    zmin = min(z, zmin);
    zmax = max(z, zmax);
  }
  auto &f = patches.front();
  u32 Imin = f.I1, Imax = f.I2, Jmin = f.J1, Jmax = f.J2, Kmin = f.K1,
      Kmax = f.K2;
  for (auto &p : patches) {
    Imin = min(p.I1, Imin);
    Imax = max(p.I2, Imax);
    Jmin = min(p.J1, Jmin);
    Jmax = max(p.J2, Jmax);
    Kmin = min(p.K1, Kmin);
    Kmax = max(p.K2, Kmax);
  }

  float Xmin = mesh.x0 + Imin * mesh.cell_size;
  float Xmax = mesh.x0 + Imax * mesh.cell_size;
  float Ymin = mesh.y0 + Jmin * mesh.cell_size;
  float Ymax = mesh.y0 + Jmax * mesh.cell_size;
  float Zmin = mesh.z0 + Kmin * mesh.cell_size;
  float Zmax = mesh.z0 + Kmax * mesh.cell_size;
  // cin >> cell_size >> x0 >> y0 >> z0;

  constexpr float r = .1f;
  for (auto p = nodes.begin(), end = nodes.end(); p != end;) {
    float x = *p++;
    float y = *p++;
    float z = *p++;
    for (const auto &patch : patches) {
      bool X = patch.IOR == 1 || patch.IOR == -1;
      bool Y = patch.IOR == 2 || patch.IOR == -2;
      bool Z = patch.IOR == 3 || patch.IOR == -3;
      float x1 = mesh.x0 + patch.I1 * mesh.cell_size - mesh.cell_size * r * X;
      float x2 = mesh.x0 + patch.I2 * mesh.cell_size + mesh.cell_size * r * X;
      float y1 = mesh.y0 + patch.J1 * mesh.cell_size - mesh.cell_size * r * Y;
      float y2 = mesh.y0 + patch.J2 * mesh.cell_size + mesh.cell_size * r * Y;
      float z1 = mesh.z0 + patch.K1 * mesh.cell_size - mesh.cell_size * r * Z;
      float z2 = mesh.z0 + patch.K2 * mesh.cell_size + mesh.cell_size * r * Z;

      if (x1 <= x && x <= x2 && y1 <= y && y <= y2 && z1 <= z && z <= z2) {
        res.insert(i);
        break;
      }
    }
    ++i;
  }
  return res;
}