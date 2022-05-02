#pragma once
#include "elements_io.hpp"
#include "fs.hpp"
#include "graphics.hpp"
#include "io.hpp"
#include "proxy.hpp"
#include "types.hpp"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <vector>

using std::conditional_t;
using std::cout;
using std::decay_t;
using std::endl;
using std::ifstream;
using std::ios;
using std::ofstream;
using std::ostream;
using std::set;
using std::setprecision;
using std::stringstream;
using std::tuple;
using std::filesystem::create_directory;
using std::filesystem::exists;

static inline u16 default_precision = 0;
static inline u16 input_precision() {
  if (default_precision)
    return default_precision;
  u16 precision = 0;
  cout << "Precision:";
  cin >> precision;
  return precision;
}
static inline void set_default_precision(u16 precision) {
  default_precision = precision;
}

static inline void as_csv(ostream &o, const vector<float> &data, u32 m, u32 n,
                          const char *sep = ",") {
  u32 i = 0;
  for (auto v : data) {
    o << v;
    ++i;
    if (i % n == 0) {
      o << endl;
    } else
      o << sep;
  }
  o << endl;
}

template <size_t sz = 2>
static inline tuple<vector<float>, u32, u32>
from_csv(istream &in, const char (&sep)[sz] = ",") {
  tuple<vector<float>, u32, u32> res;
  auto &[vec, m, n] = res;

  while (in) {
    string line;
    getline(in, line);
    stringstream ss(line);

    size_t i = 0;
    float e;
    ss >> e;
    vec.push_back(e);
    while (ss) {
      check(in, sep);
      ss >> e;
      vec.push_back(e);
      ++i;
    }
    if (i != n)
      throw file_error();
  }

  return res;
}

static inline void print_patch(const vector<float> &data, u32 m, u32 n) {
  cout << endl;
  as_csv(cout, data, m, n, " ");
}

static inline void save_patch_as_binary(const vector<float> &data,
                                        const vector<u32> &sizes,
                                        ofstream &fout) {
  write_vector_binary(fout, sizes);
  write_binary(fout, u32(0));
  write_vector_binary(fout, data);
}

static inline tuple<u32, u32> get_matrix_size(const vector<u32> &sizes) {

  if (sizes.size() > 2) {
    cout << "Data dimension " << sizes.size() << " > 2." << endl;
    return {0, 0};
  } else if (sizes.size() == 0) {
    return {0, 1};
  } else if (sizes.size() == 2) {
    auto m = sizes[0], n = sizes[1];
    return {m, n};
  } else {
    auto m = (u32)1, n = sizes[0];
    return {m, n};
  }
}

static inline void save_patch_as_csv_text(const vector<float> &data, u32 m,
                                          u32 n, ofstream &fout) {
  as_csv(fout, data, m, n);
}

static inline void reduce_dimension(vector<u32> &sizes, size_t target) {
  for (size_t i = 0; i < sizes.size() && sizes.size() > target; ++i) {
    if (sizes[i] == 1) {
      sizes.erase(sizes.cbegin() + i);
      --i;
    }
  }
}

/// @brief Get data offset in the flattened multi-dimension array.
static inline u32 array_offset(const vector<u32> &sizes,
                               const vector<u32> &index) {
  u32 offset = 0;
  auto s1 = sizes.begin(), s2 = sizes.end();
  auto i1 = index.begin(), i2 = index.end();
  for (; s1 != s2 && i1 != i2; ++s1, ++i1)
    offset *= *s1, offset += *i1;
  return offset;
}

/// @return Stride on the dimension
static inline u32 array_stride(const vector<u32> &sizes, size_t d) {
  u32 res = 1;
  auto p = sizes.begin() + d, end = sizes.end();
  for (; p != end; ++p)
    res *= *p;
  return res;
}

static inline vector<float> average(const vector<float> &data,
                                    vector<u32> &sizes, size_t dimension) {

  const auto m = sizes[dimension], n = array_stride(sizes, dimension + 1);

  vector<float> new_result;
  new_result.resize(data.size() / m);

  const auto &_data = data;
  const u32 l = (u32)data.size() / m / n;
  assert(data.size() == l * m * n);
#define SIMPLE_AVERAGE
#ifdef SIMPLE_AVERAGE
  // Will this be optimized well?
  for (size_t i = 0; i < l; ++i) {
    for (size_t k = 0; k < n; ++k) {
      for (size_t j = 0; j < m; ++j) {
        new_result[i * n + k] += data[i * m * n + j * n + k];
      }
    }
  }
#else
  for (const float *p0 = new_result.data(), // [i][]
           *const end0 = new_result.data() + new_result.size(),
                   *_p0 = data.data(); // [i][][]
       p0 != end0; p0 += n, _p0 += m * n) {
#ifndef NDEBUG
    assert(p0 >= new_result.data());
    const u32 i = (u32)(p0 - new_result.data()) / n;
    const u32 di = (u32)(p0 - new_result.data()) % n;
    assert(di == 0);
    assert(_p0 == data.data() + i * m * n);
#endif // !NDEBUG
    for (auto *p1 = p0, // [i][k]
             *const end1 = p0 + n,
              *_p1 = _p0; // [i][][k]
         p1 != end1; ++p1, ++_p1) {
#ifndef NDEBUG
      const u32 k = (u32)(p1 - p0);
      assert(0 <= k && k < n);
      assert(p1 == new_result.data() + i * n + k);
      assert(_p1 == data.data() + i * m * n + k);
#endif // !NDEBUG
      for (auto *_p2 = _p1, *const _end2 = _p1 + m * n; _p2 != _end2;
           _p2 += n) {
#ifndef NDEBUG
        const u32 j = (u32)(_p2 - _p1) / n;
        const u32 dj = (u32)(_p2 - _p1) % n;
        assert(0 <= j && j < m);
        assert(dj == 0);
        assert(p1 == new_result.data() + i * n + k);
        assert(_p2 == data.data() + i * m * n + j * n + k);

        assert(new_result.data() <= p1);
        assert(p1 < new_result.data() + new_result.size());
        assert(data.data() <= _p2);
        assert(_p2 < data.data() + data.size());
#endif // !NDEBUG
        *const_cast<float *>(p1) += *_p2; // [i][k] [i][j][k]
      }
    }
  }
#endif // SIMPLE_AVERAGE

  for (auto &e : new_result)
    e /= m;
  return new_result;
}

static inline vector<float> sample(const vector<float> &data,
                                   const vector<u32> &sizes, size_t dimension,
                                   const vector<u32> &pos) {

  const auto m = sizes[dimension], n = array_stride(sizes, dimension + 1),
             _m = (u32)pos.size();

  vector<float> new_result;
  new_result.resize(data.size() / m * _m);

  const auto &_data = data;
  const u32 l = (u32)data.size() / m / n;
  assert(data.size() == l * m * n);

  for (auto p : pos)
    assert(p < m);

  // Will this be optimized well?
  for (size_t i = 0; i < l; ++i) {
    for (size_t k = 0; k < n; ++k) {
      for (size_t _j = 0; _j < _m; ++_j) {
        new_result[i * _m * n + _j * n + k] = data[i * m * n + pos[_j] * n + k];
      }
    }
  }

  return new_result;
}

static inline void select_patch(const vector<patch_info> &patches) {
  u32 p = 0;
  cout << "Select a patch to analyze. Input an invalid index to discard."
       << endl;
#if GRAPHICS_ENABLED
  if (current < patches.size())
    cout << "Current patch selected in graphics mode is: " << current << "."
         << endl;
#endif // GRAPHICS_ENABLED

  cout << "Patch index: ";
  cin >> p;
  if (p >= patches.size()) {
    return;
  }
}

static inline void analyze(const vector<patch_info> &patches,
                           const vector<frame> &frames,
                           const vector<float> &nodes,
                           const vector<u32> &element_sizes,
                           const vector<u32> &element_indices,
                           const vector<u32> &element_index_map) {
  vector<float> data;
  vector<u32> sizes;

  auto elem_avail = [&]() {
    return !nodes.empty() && !element_sizes.empty() && !element_indices.empty();
  };

  /// polygon sizes, polygon indices, element sizes, element number ,surface
  /// numbers
  tuple<vector<u32>, vector<u32>, vector<u8>, vector<u32>, vector<u8>> polygons;
  bool polygons_calculated = false;
  auto calculate_polygon = [
    &element_sizes, &element_indices, &polygons, &element_index_map,
    &calculated = polygons_calculated
  ]() -> const auto & {
    if (!calculated) {
      if (element_index_map.empty())
        polygons = get_polygon<true, false>(element_sizes, element_indices, {});
      else
        polygons = get_polygon<true, true>(element_sizes, element_indices,
                                           element_index_map);
      calculated = true;
    }
    return polygons;
  };

  // indices_of_nodes_on_boundary
  set<u32> nodes_on_boundary;
  bool nodes_on_boundary_calculated = false;
  auto calculate_node_on_boundary = [
    &nodes, &patches, &element_sizes, &element_indices, &nodes_on_boundary,
    &calculate_polygon, &calculated = nodes_on_boundary_calculated
  ]() -> const auto & {
    if (!calculated) {
      auto &[polygon_sizes, polygon_indices, _, element_numbers,
             surface_numbers] = calculate_polygon();
      nodes_on_boundary = node_on_boundary(patches, nodes);
      calculated = true;
    }
    return nodes_on_boundary;
  };

  tuple<vector<u32>, vector<float>, vector<float>> node_data;
  bool node_data_calculated = false;
  auto calculate_node = [
        &calculated = node_data_calculated, &node_data, &calculate_polygon,
        &calculate_node_on_boundary, &patches, &nodes, &frames
  ]() -> const auto & {
    auto &indices_of_node_on_boundary = calculate_node_on_boundary();
    if (!calculated) {
      ;
      calculated = true;
    }
    return node_data;
  };

  // sizes_of_polygon_vertices_on_boundary,
  // indices_of_polygon_vertices_on_boundary,
  // numbers_of_element_whose_polygon_on_boundary
  // numbers_of_surface_on_boundary
  tuple<vector<u32>, vector<u32>, vector<u32>, vector<u8>> polygons_on_boundary;
  bool polygons_on_boundary_calculated = false;
  auto calculate_polygon_on_boundary = [
    &nodes, &patches, &element_sizes, &element_indices, &polygons_on_boundary,
    &calculate_polygon, &calculated = polygons_on_boundary_calculated
  ]() -> const auto & {
    if (!calculated) {
      auto &[polygon_sizes, polygon_indices, element_sizes, element_numbers,
             surface_numbers] = calculate_polygon();
      if (element_numbers.empty())
        polygons_on_boundary = polygon_on_boundary<false>(
            patches, nodes, polygon_sizes, polygon_indices, element_numbers,
            surface_numbers);
      else
        polygons_on_boundary = polygon_on_boundary<true>(
            patches, nodes, polygon_sizes, polygon_indices, element_numbers,
            surface_numbers);
      calculated = true;
    }
    return polygons_on_boundary;
  };

  // Polygon surface number, polygon surface centroid positions and polygon
  // average
  tuple<vector<u8>, vector<float>, vector<float>> polygons_average;
  bool average_calculated = false;
  auto calculate_polygon_average = [
        &calculated = average_calculated, &polygons_average, &calculate_polygon,
        &calculate_polygon_on_boundary, &patches, &nodes, &frames
  ]() -> const auto & {
    auto &[polygon_sizes, polygon_indices, element_sizes, element_indices,
           surface_numbers] = calculate_polygon();

    auto &[sizes_of_polygon_vertices_on_boundary,
           indices_of_polygon_vertices_on_boundary,
           indices_of_element_whose_polygon_on_boundary,
           numbers_of_surfaces_on_boundary] = calculate_polygon_on_boundary();

    if (!calculated) {
      polygons_average =
          polygon_average(patches, nodes, numbers_of_surfaces_on_boundary,
                          sizes_of_polygon_vertices_on_boundary,
                          indices_of_polygon_vertices_on_boundary, frames);
    }
    return polygons_average;
  };

  while (true) {
    cout << R"(
Commands:
m - Show current dimensions.
P - Set default float-point number precision.
w - Input in console to overlap current result.
S - Sample.
c - Read CSV file as current result.)";

    if (!data.empty())
      cout <<
          R"(
f - Flatten current result to specified dimension.
a - Calculate average on specified dimension of current result.
p - Print current result to console.
B - Save current result as binary to file.
C - Save current result as CSV file.)";

    if (!patches.empty())
      cout << R"(
s - Select a patch.
M - Visualize merged patches.)";

#if GRAPHICS_ENABLED
    if (!patches.empty() && elem_avail())
      cout << R"(
y - Visualize polygons on current patch.
n - Visualize nodes on current patch.
r - Visualize elements primitives on current patch.
Y - Visualize polygons on boundary.
N - Visualize nodes on boundary.
R - Visualize elements primitives on boundary.)";
#endif // GRAPHICS_ENABLED

    if (!patches.empty() && !frames.empty() && elem_avail())
      cout << R"(
l - Write APDL output.
A - Review output by visualizing them.)";

    if (!frames.empty() && selected_patch < patches.size())
      cout <<
          R"(
t - Copy patch data in a frame as current result.
T - Copy data in all frames of this patch as current result.)";

#if GRAPHICS_ENABLED
    if (!data.empty())
      cout << R"(
v - Visualize current result.)";
#endif // GRAPHICS_ENABLED

    cout <<
        R"(
d - Discard.
)";
    char opt = 'd';
    cin >> opt;
    switch (opt) {
    case 'd':
      return;
    case 's':
      select_patch(patches);
      break;
    case 'w': {
      cout << "Size(0 to stop): ";
      u32 sz = 0;
      sizes.clear();
      do {
        sz = 0;
        cin >> sz;
        if (sz == 0)
          break;
        sizes.push_back(sz);
      } while (true);

      u32 data_size = array_stride(sizes, 0);

      float e = 0.f;
      data.clear();
      data.reserve(data_size);
      while (data.size() < data_size) {
        cin >> e;
        data.push_back(e);
      }
    } break;
    case 'm': {
      if (sizes.empty())
        break;
      auto begin = sizes.cbegin(), end = sizes.cend();
      cout << *begin;
      ++begin;
      for (; begin != end; ++begin)
        cout << '*' << *begin;
    } break;
#if GRAPHICS_ENABLED
      /*
    case 'y':
      if (selected_patch < patches.size() && elem_avail()) {
        auto [polygon_sizes, polygon_indices, _0, _1] =
            get_polygon<false, false>(element_sizes, element_indices, {});
        auto [sizes_of_polygon_vertices_on_boundary,
              indices_of_polygon_vertices_on_boundary, _2, _3] =
            polygon_on_boundary<false>({patches[selected_patch]}, nodes,
                                       polygon_sizes, polygon_indices, {},);

        visualize_polygons(nodes, sizes_of_polygon_vertices_on_boundary,
                           indices_of_polygon_vertices_on_boundary);
      }
      break;*/
    case 'n':
      if (selected_patch < patches.size() && elem_avail()) {
        auto indices_of_node_on_boundary =
            node_on_boundary({patches[selected_patch]}, nodes);

        visualize_nodes(nodes, vector<u32>(indices_of_node_on_boundary.begin(),
                                           indices_of_node_on_boundary.end()));
      }
      break;
    case 'r':
      if (selected_patch < patches.size() && elem_avail()) {
        auto primitive_indices =
            get_primitives(nodes, element_sizes, element_indices, wireframe);
        auto primitive_indices_on_boundary = primitive_on_boundary(
            {patches[selected_patch]}, nodes, primitive_indices, wireframe);
        visualize_primitives(nodes, primitive_indices_on_boundary);
      }
      break;
    case 'Y':
      if (!patches.empty() && elem_avail()) {
        auto &[polygon_sizes, polygon_indices, _0, _1, _2] =
            calculate_polygon();
        auto &[sizes_of_polygon_vertices_on_boundary,
               indices_of_polygon_vertices_on_boundary,
               indices_of_element_whose_polygon_on_boundary, surface_numbers] =
            calculate_polygon_on_boundary();

        visualize_polygons(nodes, sizes_of_polygon_vertices_on_boundary,
                           indices_of_polygon_vertices_on_boundary);
      }
      break;
    case 'N':
      if (!patches.empty() && elem_avail()) {
        auto indices_of_node_on_boundary = node_on_boundary(patches, nodes);

        visualize_nodes(nodes, vector<u32>(indices_of_node_on_boundary.begin(),
                                           indices_of_node_on_boundary.end()));
      }
      break;
    case 'R':
      if (!patches.empty() && elem_avail()) {
        auto primitive_indices =
            get_primitives(nodes, element_sizes, element_indices, wireframe);
        auto primitive_indices_on_boundary =
            primitive_on_boundary(patches, nodes, primitive_indices, wireframe);
        visualize_primitives(nodes, primitive_indices_on_boundary);
      }
      break;
    case 'p': {
      u16 precision = input_precision();
      cout << setprecision(precision);

      auto [m, n] = get_matrix_size(sizes);
      if (n != 0)
        print_patch(data, m, n);
    } break;
    case 'A':
      if (!patches.empty() && elem_avail()) {
        auto &[_, centroid, boundary_data] = calculate_polygon_average();
        if (centroid.empty() || boundary_data.empty()) {
          cerr << "No surfaces on boundary found\n";
          break;
        }
        assert(frames.size() <= numeric_limits<u32>::max());
        visualize_nodes(centroid, boundary_data, (u32)frames.size());
      }
      break;
    case 'M':
      if (!patches.empty()) {
        visualize_regions(merge(patches));
      }
      break;
#endif // GRAPHICS_ENABLED
    case 'P': {
      u16 precision = input_precision();
      set_default_precision(precision);
    } break;
    case 'l':
      if (!patches.empty() && !frames.empty() && elem_avail()) {
        auto opt1 = request_file_by_name([](const path &p) { return true; },
                                         "APDL output");
        auto opt2 = request_file_by_name([](const path &p) { return true; },
                                         "APDL output directory");
        if (!opt1)
          break;
        auto out = ofstream(opt1.value());
        if (!out) {
          cout << "Open Failed." << endl;
        }
        if (!opt2)
          break;
        auto &dir = opt2.value();
        if (!exists(dir)) {
          auto suc = create_directory(dir);
          if (!suc) {
            cout << "Create Failed." << endl;
          }
        }
        out << "/PREP7" << endl;

        const auto N = frames.size();

        auto &[on_boundary_vertex_sizes, on_boundary_vertex_indices,
               element_numbers, surface_numbers] =
            calculate_polygon_on_boundary();
        auto &[on_boundary_surface_numbers, _2, boundary_data] =
            calculate_polygon_average();

        assert(on_boundary_surface_numbers.size() ==
               on_boundary_vertex_sizes.size());

        vector<vector<u32>::const_iterator> ps;
        ps.reserve(on_boundary_surface_numbers.size() + 1);
        vector<vector<float>::const_iterator> Ps;
        Ps.reserve(on_boundary_surface_numbers.size() + 1);
        {
          auto p = on_boundary_vertex_indices.begin();
          auto e = on_boundary_vertex_indices.end();
          auto P = boundary_data.begin();
          auto E = boundary_data.end();

          for (auto sz : on_boundary_vertex_sizes) {
            assert(p < e);
            assert(P < E);
            ps.push_back(p);
            Ps.push_back(P);
            p += sz;
            P += N;
          }
          ps.push_back(p);
          Ps.push_back(P);
          assert(p == e);
          assert(P == E);
        }

        assert(element_numbers.size() == on_boundary_surface_numbers.size());
#pragma omp parallel for schedule(dynamic, 5000)
        for (long i = 0; i < on_boundary_vertex_sizes.size(); ++i) {
          write_table(out, opt2.value(), "HFLUX", element_numbers[i],
                      on_boundary_surface_numbers[i], ps[i], ps[i + 1], frames,
                      Ps[i], Ps[i + 1]);
        }
        out << "FINISH" << endl
            << "/SOLU" << endl
            << "ALLSEL,ALL" << endl
            << "/NERR,,99999999" << endl;
      }
      break;
    case 'S': {
      size_t dimension = -1;
      cin >> dimension;
      if (dimension >= sizes.size() || sizes.size() == 1) {
        cout << "Dimension invalid." << endl;
        break;
      }
      const auto sz = sizes[dimension];
      cout << "Size of this dimension: " << sz << endl;
      vector<u32> pos;
      while (true) {
        u32 i = -1;
        cin >> i;
        if (i >= sz)
          break;
        pos.push_back(i);
      }

      data = sample(data, sizes, dimension, pos);
      sizes[dimension] = (u32)pos.size();

    } break;
#if GRAPHICS_ENABLED
    case 'v': {
      if (!data.empty()) {
        auto [m, n] = get_matrix_size(sizes);
        if (!n)
          break;
        visualize_data(data, m, n);
      }
    } break;
#endif // GRAPHICS_ENABLED
    case 'B': {
      cout << "File name: ";
      string path;
      getline(cin, path);
      ofstream fout(path, ios::binary);
      if (!fout) {
        cout << "Failed to open file." << endl;
        break;
      }
      save_patch_as_binary(data, sizes, fout);
    } break;
    case 'c': {
      auto opt = request_file_by_name(
          [](const path &p) { return is_directory(p); }, "csv file");
      if (opt.has_value()) {
        auto &path = opt.value();
        ifstream fin = ifstream(path);
        if (!fin) {
          cout << "Failed to open file." << endl;
          break;
        }

        auto [_data, _m, _n] = from_csv(fin);
        data = std::move(_data);
        sizes = {_m, _n};
      }
    } break;
    case 'C': {
      optional<path> opt = request_file_by_name(
          [](const path &p) { return exists(p); }, "existed");
      if (!opt)
        break;
      auto &path = opt.value();
      ofstream fout(path);
      if (!fout) {
        cout << "Failed to open file." << endl;
        break;
      }

      u16 precision = input_precision();
      fout << setprecision(precision);

      auto [m, n] = get_matrix_size(sizes);
      save_patch_as_csv_text(data, m, n, fout);
    } break;
    case 'f': {
      size_t d = 2;
      cin >> d;
      reduce_dimension(sizes, d);
      if (sizes.size() > d)
        cout << "Flattening failed. Current dimension is " << sizes.size()
             << "." << endl;
    } break;
    case 'a': {
      // Consider [l][m][n]
      size_t dimension = -1;
      cin >> dimension;
      if (dimension >= sizes.size() || sizes.size() == 1) {
        cout << "Dimension invalid." << endl;
        break;
      }
      data = average(data, sizes, dimension);
      sizes.erase(sizes.cbegin() + dimension);

    } break;
    case 't':
      if (!frames.empty() && selected_patch < patches.size()) {
        const auto &patch = patches[selected_patch];
        cout << "Frame index: ";
        auto f = frames.size();
        cin >> f;
        if (f >= frames.size()) {
          cout << "Frame index invalid." << endl;
          break;
        }
        data = frames[f].data[selected_patch].data;
        sizes = {patch.K(), patch.J(), patch.I()};
      }
      break;
    case 'T':
      if (!frames.empty() && selected_patch < patches.size()) {
        const auto &patch = patches[selected_patch];
        sizes = {(u32)frames.size(), patch.K(), patch.J(), patch.I()};
        data.clear();
        data.reserve(array_stride(sizes, 0));
        for (const auto &frame : frames) {
          const auto &patch = frame.data[selected_patch];
          data.insert(data.cend(), patch.data.cbegin(), patch.data.cend());
        }
      }
      break;
    default:
      break;
    }
  }
}