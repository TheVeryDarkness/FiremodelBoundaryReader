#pragma once

#include "element_basic.hpp"
#include "fds_basic.hpp"
#include <cassert>

using std::all_of;
using std::int_fast64_t;
using std::make_tuple;
using std::minmax_element;
using std::uint_fast8_t;

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
      from_elements(nodes, vertex_count, elements, wireframe, ratio);
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

static inline tuple<tuple<vector<float>>, vector<u32>>
from_matrix_data(const vector<float> &data, u32 n, bool wireframe) {
  tuple<tuple<vector<float>>, vector<u32>> res;
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
      } else if (j != 0) {
        //   1
        // 2 3
        indices.push_back(i * n + j);
        indices.push_back((i - 1) * n + j - 1);
        indices.push_back(i * n + j - 1);
        // 2 1
        // 3
        indices.push_back(i * n + j);
        indices.push_back((i - 1) * n + j);
        indices.push_back((i - 1) * n + j - 1);
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

static const float tolerance = .01f;

static inline void check_border(const vector<patch_info> &patches,
                                const vector<float> &nodes) {
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
}

static inline set<u32> node_on_boundary(const vector<patch_info> &patches,
                                        const vector<float> &nodes) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  set<u32> res;

  check_border(patches, nodes);

  u32 i = 0;

  for (auto p = nodes.begin(), end = nodes.end(); p != end;) {
    float x = *p++;
    float y = *p++;
    float z = *p++;
    for (const auto &patch : patches) {
      float x1 =
          mesh.x0 + patch.I1 * mesh.cell_size - mesh.cell_size * tolerance;
      float x2 =
          mesh.x0 + patch.I2 * mesh.cell_size + mesh.cell_size * tolerance;
      float y1 =
          mesh.y0 + patch.J1 * mesh.cell_size - mesh.cell_size * tolerance;
      float y2 =
          mesh.y0 + patch.J2 * mesh.cell_size + mesh.cell_size * tolerance;
      float z1 =
          mesh.z0 + patch.K1 * mesh.cell_size - mesh.cell_size * tolerance;
      float z2 =
          mesh.z0 + patch.K2 * mesh.cell_size + mesh.cell_size * tolerance;

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

  check_border(patches, nodes);

  for (auto p = nodes.begin(), end = nodes.end(); p != end;) {
    float x = *p++;
    float y = *p++;
    float z = *p++;
    size_t i_patch = 0;
    for (const auto &patch : patches) {
      float x1 =
          mesh.x0 + patch.I1 * mesh.cell_size - mesh.cell_size * tolerance;
      float x2 =
          mesh.x0 + patch.I2 * mesh.cell_size + mesh.cell_size * tolerance;
      float y1 =
          mesh.y0 + patch.J1 * mesh.cell_size - mesh.cell_size * tolerance;
      float y2 =
          mesh.y0 + patch.J2 * mesh.cell_size + mesh.cell_size * tolerance;
      float z1 =
          mesh.z0 + patch.K1 * mesh.cell_size - mesh.cell_size * tolerance;
      float z2 =
          mesh.z0 + patch.K2 * mesh.cell_size + mesh.cell_size * tolerance;

      if (in_box(x, y, z, x1, x2, y1, y2, z1, z2)) {
        res[i_patch].insert(i);
      }
      ++i_patch;
    }
    ++i;
  }
  return res;
}

static inline vector<set<u32>>
node_on_each_connected_patch_boundary(const vector<patch_info> &patches,
                                      const vector<float> &nodes) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  vector<set<u32>> res;
  res.resize(patches.size());

  check_border(patches, nodes);

  auto regions = merge(patches);

  return regions.filter_for_each_connected_region(nodes);
}

static inline vector<u32> primitive_on_boundary(
    const vector<patch_info> &patches, const vector<float> &nodes,
    const vector<u32> &primitive_indices, const bool wireframe) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  vector<u32> res;

  auto _nodes = node_on_each_connected_patch_boundary(patches, nodes);

  // Loop for nodes on each boundary
  for (const auto &set : _nodes) {
    if (wireframe)
      for (auto p = primitive_indices.begin(), end = primitive_indices.end();
           p != end;) {
        u32 i0 = *p++;
        u32 i1 = *p++;
        assert(i0 < nodes.size() / 3);
        assert(i1 < nodes.size() / 3);
        if (set_contains_all(set, {i0, i1})) {
          res.push_back(i0);
          res.push_back(i1);
        }
      }
    else
      for (auto p = primitive_indices.begin(), end = primitive_indices.end();
           p != end;) {
        u32 i0 = *p++;
        u32 i1 = *p++;
        u32 i2 = *p++;
        assert(i0 < nodes.size() / 3);
        assert(i1 < nodes.size() / 3);
        assert(i2 < nodes.size() / 3);
        if (set_contains_all(set, {i0, i1, i2})) {
          res.push_back(i0);
          res.push_back(i1);
          res.push_back(i2);
        }
      }
  }
  return res;
}

template <bool withElementNumber>
static inline tuple<vector<u32>, vector<u32>, vector<u32>> polygon_on_boundary(
    const vector<patch_info> &patches, const vector<float> &nodes,
    const vector<u32> &polygon_sizes, const vector<u32> &polygon_indices,
    const vector<u32> &element_index_map) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  auto sum = accumulate(polygon_sizes.begin(), polygon_sizes.end(), 0);
  assert(sum == polygon_indices.size());
  assert(polygon_sizes.size() == element_index_map.size());
  tuple<vector<u32>, vector<u32>, vector<u32>> res;
  auto &[sizes, indices, numbers] = res;

  auto _nodes = node_on_each_connected_patch_boundary(patches, nodes);

  vector<u32> polygon_vertex_indices;

  for (const auto &set : _nodes) {
    auto p = polygon_indices.begin();
    auto e = polygon_indices.end();

    u32 i_polygon = 0;
    // For each polygon
    for (auto sz : polygon_sizes) {
      assert(p < e);
      polygon_vertex_indices.clear();
      polygon_vertex_indices.insert(polygon_vertex_indices.cend(), p, p + sz);

      bool in = true;
      // For each vertex of current polygon
      for (auto index : polygon_vertex_indices) {
        if (!set_contains(set, index)) {
          in = false;
          break;
        }
      }

      if (in) {
        for (auto i : polygon_vertex_indices) {
          assert(i < nodes.size() / 3);
          indices.push_back(i);
        }
        if constexpr (withElementNumber)
          numbers.push_back(element_index_map[i_polygon]);
        sizes.push_back(sz);
      }

      p += sz;
      ++i_polygon;
    }
  }
  return res;
}

static inline int_fast8_t sign(int_fast64_t i) {
  return i > 0 ? 1 : i < 0 ? -1 : 0;
}

static inline bool inside(const vector<u32> &I, const vector<u32> &J, u32 i,
                          u32 j) {
  using i64f = int_fast64_t;
  uint_fast8_t last = [&I, &J, i, j]() {
    i64f i1 = I.back();
    i64f i2 = I.front();
    i64f j1 = J.back();
    i64f j2 = J.front();
    auto dir = (i2 - i1) * (j - j1) - (j2 - j1) * (i - i1);
    return sign(dir);
  }();
  //  Pi+1Pi cross PPi
  bool res = true;
  for (size_t index = 0; index + 1 < I.size(); ++index) {
    i64f i1 = I[index];
    i64f i2 = I[index + 1];
    i64f j1 = J[index];
    i64f j2 = J[index + 1];
    auto dir = (i2 - i1) * (j - j1) - (j2 - j1) * (i - i1);
    res = res && (sign(dir) * last >= 0);
  }
  return res;
}

static inline tuple<vector<u32>, vector<u32>, vector<u32>>
mesh_coordinates(const vector<float> &vec, const vector<u32> &indices) {
  tuple<vector<u32>, vector<u32>, vector<u32>> res;
  auto &[x, y, z] = res;
  for (auto i : indices) {
    x.push_back((vec[3 * i] - mesh.x0) / mesh.cell_size);
    y.push_back((vec[3 * i + 1] - mesh.y0) / mesh.cell_size);
    z.push_back((vec[3 * i + 2] - mesh.z0) / mesh.cell_size);
  }
  return res;
}

template <size_t dim1, size_t dim2>
static inline void find(const u32 i1, const u32 i2, const u32 j1, const u32 j2,
                        const vector<frame> &frames, u32 i_patch,
                        const patch_info &patch, const vector<u32> &I,
                        const vector<u32> &J, vector<float> &sum, u32 &count) {
  assert(count == 0);
  assert(sum.size() == frames.size());
  auto [I1, I2] = patch.border<dim1>();
  auto [J1, J2] = patch.border<dim2>();
  auto imin = max(I1, i1);
  auto imax = min(I2, i2);
  auto jmin = max(J1, j1);
  auto jmax = min(J2, j2);
  for (u32 i = imin; i <= imax; ++i) {
    for (u32 j = jmin; j <= jmax; ++j) {
      if (inside(I, J, i, j)) {
        for (size_t index = 0; index < frames.size(); ++index) {
          sum[index] += frames[index]
                            .data[i_patch]
                            .data[i - I1 + patch.length<dim1>() * (j - J1)];
        }
        count += 1;
      }
    }
  }
}

static inline float average(const vector<u32> &vec) {
  return (accumulate(vec.begin(), vec.end(), 0.f) / vec.size());
}

/// @return Polygon surface number and polygon average
static inline tuple<vector<u32>, vector<float>, vector<float>> polygon_average(
    const vector<patch_info> &patches, const vector<float> &nodes,
    const vector<u8> &element_sizes, const vector<u32> &polygon_sizes,
    const vector<u32> &polygon_indices, const vector<frame> &frames) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  size_t none = 0;
  size_t degenerated = 0;
  auto polygon_indices_count =
      accumulate(polygon_sizes.begin(), polygon_sizes.end(), 0);
  assert(polygon_indices_count == polygon_indices.size());
  tuple<vector<u32>, vector<float>, vector<float>> res;
  auto &[numbers, positions, data] = res;

  auto _nodes = node_on_each_patch_boundary(patches, nodes);

  constexpr float r = .01f;

  // data.reserve(frames.size() * polygon_indices.size());

  // Take variables outside.
  vector<float> sum;
  vector<u32> polygon_vertex_indices;

  u32 i_patch = 0;
  for (const auto &set : _nodes) {
    auto &patch = patches[i_patch];
    auto p = polygon_indices.begin();
    auto e = polygon_indices.end();
    auto P = element_sizes.begin();
    auto E = element_sizes.end();

    u8 i_surface = 0;
    for (auto sz : polygon_sizes) {
      assert(p < e);
      assert(P < E);
      polygon_vertex_indices.clear();
      polygon_vertex_indices.insert(polygon_vertex_indices.end(), p, p + sz);

      bool in = true;
      for (auto index : polygon_vertex_indices) {
        if (!set_contains(set, index))
          in = false;
      }

      if (in) {
        const auto &[I, J, K] = mesh_coordinates(nodes, polygon_vertex_indices);

        numbers.push_back(i_surface);
        positions.push_back(average(I));
        positions.push_back(average(J));
        positions.push_back(average(K));

        const auto [_i1, _i2] = minmax_element(I.cbegin(), I.cend());
        const auto [_j1, _j2] = minmax_element(J.cbegin(), J.cend());
        const auto [_k1, _k2] = minmax_element(K.cbegin(), K.cend());
        auto i1 = *_i1;
        auto i2 = *_i2;
        auto j1 = *_j1;
        auto j2 = *_j2;
        auto k1 = *_k1;
        auto k2 = *_k2;

        sum.clear();
        sum.resize(frames.size(), 0.f);
        u32 count = 0;

        switch (patches[i_patch].IOR) {
        case 1:
        case -1:
          find<1, 2>(j1, j2, k1, k2, frames, i_patch, patch, J, K, sum, count);
          break;
        case 2:
        case -2:
          find<0, 2>(i1, i2, k1, k2, frames, i_patch, patch, I, K, sum, count);
          break;
        case 3:
        case -3:
          find<0, 1>(i1, i2, j1, j2, frames, i_patch, patch, I, J, sum, count);
          break;
        default:
          assert(false);
          break;
        }

        if (count > 0) {
          for (auto &s : sum)
            s /= count;
        } else {
          if (all_same(I) || all_same(J) || all_same(K))
            ++degenerated;
          else
            ++none;
        }

        data.insert(data.end(), sum.cbegin(), sum.cend());
      }

      p += sz;
      ++i_surface;
      if (i_surface == *P) {
        i_surface = 0;
        ++P;
      }
    }

    ++i_patch;
  }
  // data.shrink_to_fit();
  if (degenerated)
    clog << degenerated << " surfaces seem to be degenerated.\n";
  if (none)
    cerr << none
         << " non-degenerated surfaces on boundary don't contain data "
            "point.\n";
  return res;
}