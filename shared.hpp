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

/// @brief
/// @param patches
/// @param nodes
/// @param polygon_sizes
/// @param polygon_indices
/// @param node_numbers_map
/// @param element_numbers_map
/// @param polygon_surface_numbers
/// @return Polygons vertex counts, polygon vertex indices, node numbers,
/// element numbers, surface numbers
static inline tuple<vector<u32>, vector<u32>, vector<u32>, vector<u32>,
                    vector<u8>>
polygon_on_boundary(const vector<patch_info> &patches,
                    const vector<float> &nodes,
                    const vector<u32> &polygon_sizes,
                    const vector<u32> &polygon_indices,
                    const vector<u32> &node_numbers_map,
                    const vector<u32> &element_numbers_map,
                    const vector<u8> &polygon_surface_numbers) {
  const bool withNodeNumber = !node_numbers_map.empty();
  const bool withElementNumber = !element_numbers_map.empty();

  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  assert(accumulate(polygon_sizes.begin(), polygon_sizes.end(), 0) ==
         polygon_indices.size());
  if (withNodeNumber)
    assert(polygon_indices.size() == node_numbers_map.size());
  if (withElementNumber)
    assert(polygon_sizes.size() == element_numbers_map.size());
  assert(polygon_surface_numbers.size() == polygon_sizes.size());
  tuple<vector<u32>, vector<u32>, vector<u32>, vector<u32>, vector<u8>> res;
  auto &[sizes, indices, node_numbers, elem_numbers, load_keys] = res;

  auto _nodes = node_on_each_connected_patch_boundary(patches, nodes);

  vector<u32> polygon_vertex_indices;

  for (const auto &set : _nodes) {
    auto p = polygon_indices.begin();
    auto e = polygon_indices.end();
    auto P = polygon_surface_numbers.begin();
    auto E = polygon_surface_numbers.end();

    u32 i_polygon = 0;
    // For each polygon
    for (auto sz : polygon_sizes) {
      assert(p < e);
      polygon_vertex_indices.clear();
      polygon_vertex_indices.insert(polygon_vertex_indices.cend(), p, p + sz);

      // For each vertex of current polygon
      const bool in = set_contains_all(set, polygon_vertex_indices);

      if (in) {
        for (auto i : polygon_vertex_indices) {
          assert(i < nodes.size() / 3);
          indices.push_back(i);
          if (withNodeNumber)
            node_numbers.push_back(node_numbers_map[i_polygon]);
        }
        if (withElementNumber)
          elem_numbers.push_back(element_numbers_map[i_polygon]);
        sizes.push_back(sz);
        load_keys.push_back(*P);
      }

      p += sz;
      ++P;
      ++i_polygon;
    }
    assert(i_polygon == polygon_sizes.size());
    assert(p == e);
    assert(P == E);
  }
  assert(accumulate(sizes.cbegin(), sizes.cend(), 0) == indices.size());
  if (withNodeNumber)
    assert(node_numbers.size() == indices.size());
  if (withElementNumber)
    assert(elem_numbers.size() == sizes.size());
  return res;
}

static inline int_fast8_t sign(float i) { return i > 0 ? 1 : i < 0 ? -1 : 0; }

static inline bool inside(const vector<float> &I, const vector<float> &J, u32 i,
                          u32 j) {
  using T = float;
  uint_fast8_t last = [&I, &J, i, j]() {
    T i1 = I.back();
    T i2 = I.front();
    T j1 = J.back();
    T j2 = J.front();
    auto dir = (i2 - i1) * (j - j1) - (j2 - j1) * (i - i1);
    return sign(dir);
  }();
  //  PiPi+1 cross PiP
  bool res = true;
  for (size_t index = 0; index + 1 < I.size(); ++index) {
    T i1 = I[index];
    T i2 = I[index + 1];
    T j1 = J[index];
    T j2 = J[index + 1];
    auto dir = (i2 - i1) * (j - j1) - (j2 - j1) * (i - i1);
    res = res && (sign(dir) * last >= 0);
  }
  return res;
}

static inline tuple<vector<float>, vector<float>, vector<float>>
mesh_coordinates(const vector<float> &vec, const vector<u32> &indices) {
  tuple<vector<float>, vector<float>, vector<float>> res;
  auto &[x, y, z] = res;
  for (auto i : indices) {
    x.push_back((vec[3 * i] - mesh.x0) / mesh.cell_size);
    y.push_back((vec[3 * i + 1] - mesh.y0) / mesh.cell_size);
    z.push_back((vec[3 * i + 2] - mesh.z0) / mesh.cell_size);
  }
  return res;
}

static inline tuple<vector<float>, vector<float>, vector<float>>
round_all(const vector<float> &I, const vector<float> &J,
          const vector<float> &K) {
  tuple<vector<float>, vector<float>, vector<float>> res;
  auto &[i, j, k] = res;
  for (size_t index = 0; index < I.size(); ++index) {
    i.push_back(floorf(I[index]));
    j.push_back(floorf(J[index]));
    k.push_back(floorf(K[index]));
  }
  return res;
}

[[nodiscard]] static inline vector<float>
calculate_node(const vector<frame> &frames, const vector<patch_info> &patches,
               u32 i, u32 j, u32 k) {
  vector<float> res;
  for (size_t i_patch = 0; i_patch < patches.size(); ++i_patch) {
    auto &patch = patches[i_patch];
    auto [I1, I2] = patch.border<0>();
    auto [J1, J2] = patch.border<1>();
    auto [K1, K2] = patch.border<2>();
    if (I1 <= i && i <= I2 && J1 <= j && j <= J2 && K1 <= k && k <= K2) {
      for (auto &frame : frames)
        res.push_back(
            frame.data[i_patch].data[(i - I1) + patch.I() * (j - J1) +
                                     patch.I() * patch.J() * (k - K1)]);
      return res;
    }
  }
  assert(false);
  return res;
}

template <size_t dim0, size_t dim1, size_t dim2>
[[nodiscard]] static inline vector<float>
find(const vector<frame> &frames, const vector<patch_info> &patches,
     const vector<float> &I, const vector<float> &J, const vector<float> &K) {
  assert(I.size() == J.size());
  vector<float> res;
  for (const auto &frame : frames) {
    float sum = 0;
    for (size_t index = 0; index < I.size(); ++index) {
      u32 i = lroundf(I[index]);
      u32 j = lroundf(J[index]);
      u32 k = lroundf(K[index]);
      for (size_t i_patch = 0; i_patch < patches.size(); ++i_patch) {
        auto &patch = patches[i_patch];
        auto [I1, I2] = patch.border<dim1>();
        auto [J1, J2] = patch.border<dim2>();
        auto [K1, K2] = patch.border<dim0>();
        assert(K1 == K1);
        if (k == K1 && I1 <= i && i <= I2 && J1 <= j && j <= J2) {
          for (auto &frame : frames)
            sum += frame.data[i_patch]
                       .data[(i - I1) + patch.length<dim1>() * (j - J1)];
          goto NEXT;
        }
      }
      assert(false);
    NEXT:;
    }
    res.push_back(sum / I.size());
  }
  return res;
}

static inline float average(const vector<float> &vec) {
  return (accumulate(vec.begin(), vec.end(), 0.f) / vec.size());
}

static inline void log_points(const vector<float> &I, const vector<float> &J,
                              const vector<float> &K) {
  assert(I.size() == J.size());
  assert(I.size() == K.size());
  for (size_t i = 0; i < I.size(); ++i)
    clog << setw(8) << I[i] << ' ' << setw(8) << J[i] << ' ' << setw(8) << K[i]
         << '\n';
}

/// @return Polygon surface number, polygon surface centroid positions and
/// polygon average
static inline tuple<vector<u8>, vector<float>, vector<float>>
polygon_average(const vector<patch_info> &patches, const vector<float> &nodes,
                const vector<u8> &polygon_surface_numbers,
                const vector<u32> &on_boundary_polygon_sizes,
                const vector<u32> &on_boundary_polygon_indices,
                const vector<frame> &frames) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  size_t none = 0;
  size_t degenerated = 0;
  auto polygon_indices_count = accumulate(on_boundary_polygon_sizes.begin(),
                                          on_boundary_polygon_sizes.end(), 0);
  assert(polygon_indices_count == on_boundary_polygon_indices.size());
  tuple<vector<u8>, vector<float>, vector<float>> res;
  auto &[numbers, positions, data] = res;

  // data.reserve(frames.size() * polygon_indices.size());

  // Take variables outside.
  vector<float> sum;
  vector<u32> polygon_vertex_indices;

  // For the set of nodes on each connected region.
  auto p = on_boundary_polygon_indices.begin();
  auto e = on_boundary_polygon_indices.end();
  auto P = polygon_surface_numbers.begin();
  auto E = polygon_surface_numbers.end();

  // For each polygon
  for (auto sz : on_boundary_polygon_sizes) {
    assert(p < e);
    assert(P < E);
    polygon_vertex_indices.clear();
    polygon_vertex_indices.insert(polygon_vertex_indices.end(), p, p + sz);

    const auto &[I, J, K] = mesh_coordinates(nodes, polygon_vertex_indices);
    const auto &[i, j, k] = round_all(I, J, K);

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

    for (size_t index = 0; index < I.size(); ++index) {
      auto _i = I[index];
      auto _j = J[index];
      auto _k = K[index];
      auto n = calculate_node(frames, patches, lroundf(_i), lroundf(_j),
                              lroundf(_k));
      for (size_t t = 0; t < frames.size(); ++t)
        sum[t] += n[t];
    }
    for (size_t t = 0; t < frames.size(); ++t)
      sum[t] /= I.size();

    numbers.push_back(*P);
    positions.push_back(average(I));
    positions.push_back(average(J));
    positions.push_back(average(K));
    data.insert(data.end(), sum.cbegin(), sum.cend());
    p += sz;
    ++P;
  }
  assert(p == e);
  assert(P == E);

  // data.shrink_to_fit();
  if (degenerated)
    clog << degenerated << " surfaces seem to be degenerated.\n";
  if (none)
    cerr << none
         << " non-degenerated surfaces on boundary don't contain data "
            "point.\n";
  return res;
}
