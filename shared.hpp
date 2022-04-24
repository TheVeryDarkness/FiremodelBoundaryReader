#pragma once

#include "element_basic.hpp"
#include "fds_basic.hpp"
#include <cassert>

using std::make_tuple;

static inline tuple<tuple<vector<float>>, vector<u32>>
from_patches_and_elements(const vector<float> &nodes,
                          const vector<u32> &vertex_count,
                          const vector<u32> &elements,
                          const vector<patch_info> &patches,
                          const bool wireframe, const u32 sum,
                          const float ratio) {
  auto [_, i1] = from_patches<false>(patches, wireframe);
  auto &[p1, __] = _;
  auto [___, i2] =
      from_elements(nodes, vertex_count, elements, wireframe, sum, ratio);
  auto &[p2] = ___;

  p2.reserve(p1.size() + p2.size());
  assert(p2.size() / 3 < numeric_limits<u32>::max());
  assert(p2.size() % 3 == 0);
  u32 base = u32(p2.size() / 3);

  for (auto p = p1.cbegin(); p != p1.cend();) {
    auto i = *p++;
    auto j = *p++;
    auto k = *p++;
    p2.push_back(float(i * mesh.cell_size + mesh.x0) * ratio);
    p2.push_back(float(j * mesh.cell_size + mesh.y0) * ratio);
    p2.push_back(float(k * mesh.cell_size + mesh.z0) * ratio);
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

static inline bool in_box(float x, float y, float z, float x1, float x2,
                          float y1, float y2, float z1, float z2) {
  return x1 <= x && x <= x2 && y1 <= y && y <= y2 && z1 <= z && z <= z2;
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

  cout << "Nodes:  [" << xmin << ", " << xmax << "]*[" << ymin << ", " << ymax
       << "]*[" << zmin << ", " << zmax << "]" << endl;
  cout << "Meshes: [" << Xmin << ", " << Xmax << "]*[" << Ymin << ", " << Ymax
       << "]*[" << Zmin << ", " << Zmax << "]" << endl;

  constexpr float r = .01f;
  for (auto p = nodes.begin(), end = nodes.end(); p != end;) {
    float x = *p++;
    float y = *p++;
    float z = *p++;
    for (const auto &patch : patches) {
      float x1 = mesh.x0 + patch.I1 * mesh.cell_size - mesh.cell_size * r;
      float x2 = mesh.x0 + patch.I2 * mesh.cell_size + mesh.cell_size * r;
      float y1 = mesh.y0 + patch.J1 * mesh.cell_size - mesh.cell_size * r;
      float y2 = mesh.y0 + patch.J2 * mesh.cell_size + mesh.cell_size * r;
      float z1 = mesh.z0 + patch.K1 * mesh.cell_size - mesh.cell_size * r;
      float z2 = mesh.z0 + patch.K2 * mesh.cell_size + mesh.cell_size * r;

      if (in_box(x, y, z, x1, x2, y1, y2, z1, z2)) {
        res.insert(i);
        break;
      }
    }
    ++i;
  }
  return res;
}

static inline vector<set<u32>>
node_on_each_patch_boundary(const vector<patch_info> &patches,
                            const vector<float> &nodes) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  vector<set<u32>> res;
  res.resize(patches.size());

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

  cout << "Nodes:  [" << xmin << ", " << xmax << "]*[" << ymin << ", " << ymax
       << "]*[" << zmin << ", " << zmax << "]" << endl;
  cout << "Meshes: [" << Xmin << ", " << Xmax << "]*[" << Ymin << ", " << Ymax
       << "]*[" << Zmin << ", " << Zmax << "]" << endl;

  constexpr float r = .01f;
  for (auto p = nodes.begin(), end = nodes.end(); p != end;) {
    float x = *p++;
    float y = *p++;
    float z = *p++;
    size_t i_patch = 0;
    for (const auto &patch : patches) {
      float x1 = mesh.x0 + patch.I1 * mesh.cell_size - mesh.cell_size * r;
      float x2 = mesh.x0 + patch.I2 * mesh.cell_size + mesh.cell_size * r;
      float y1 = mesh.y0 + patch.J1 * mesh.cell_size - mesh.cell_size * r;
      float y2 = mesh.y0 + patch.J2 * mesh.cell_size + mesh.cell_size * r;
      float z1 = mesh.z0 + patch.K1 * mesh.cell_size - mesh.cell_size * r;
      float z2 = mesh.z0 + patch.K2 * mesh.cell_size + mesh.cell_size * r;

      if (in_box(x, y, z, x1, x2, y1, y2, z1, z2)) {
        res[i_patch].insert(i);
        break;
      }
      ++i_patch;
    }
    ++i;
  }
  return res;
}

static inline vector<u32>
primitive_on_boundary(const vector<patch_info> &patches,
                      const vector<float> &nodes, const vector<u32> &elements,
                      const bool wireframe) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  vector<u32> res;

  auto _nodes = node_on_each_patch_boundary(patches, nodes);

  constexpr float r = .01f;

  for (const auto &set : _nodes) {
    if (wireframe)
      for (auto p = elements.begin(), end = elements.end(); p != end;) {
        u32 i0 = *p++;
        u32 i1 = *p++;
        assert(i0 < nodes.size());
        assert(i1 < nodes.size());
        if (set_contains_all(set, {i0, i1})) {
          res.push_back(i0);
          res.push_back(i1);
          break;
        }
      }
    else
      for (auto p = elements.begin(), end = elements.end(); p != end;) {
        u32 i0 = *p++;
        u32 i1 = *p++;
        u32 i2 = *p++;
        assert(i0 < nodes.size());
        assert(i1 < nodes.size());
        assert(i2 < nodes.size());
        if (set_contains_all(set, {i0, i1, i2})) {
          res.push_back(i0);
          res.push_back(i1);
          res.push_back(i2);
          break;
        }
      }
  }
  return res;
}

static inline tuple<vector<u32>, vector<u32>>
polygon_on_boundary(const vector<patch_info> &patches,
                    const vector<float> &nodes,
                    const vector<u32> &polygon_sizes,
                    const vector<u32> &polygon_indices, const bool wireframe) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  auto sum = accumulate(polygon_sizes.begin(), polygon_sizes.end(), 0);
  assert(sum == polygon_indices.size());
  tuple<vector<u32>, vector<u32>> res;
  auto &[sizes, indices] = res;

  auto _nodes = node_on_each_patch_boundary(patches, nodes);

  constexpr float r = .01f;

  auto &polygon_map = get_element_polygons();

  for (const auto &set : _nodes) {

    auto p = polygon_indices.begin();
    auto e = polygon_indices.end();
    for (auto sz : polygon_sizes) {
      assert(p < e);
      vector<u32> polygon_vertex_indices = {p, p + sz};

      bool in = true;
      for (auto index : polygon_vertex_indices) {
        if (!set_contains(set, index))
          in = false;
      }

      if (in) {
        for (auto i : polygon_vertex_indices)
          indices.push_back(i);
        sizes.push_back(sz);
      }
      p += sz;
    }
  }
  return res;
}