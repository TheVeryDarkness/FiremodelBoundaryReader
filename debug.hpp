#pragma once

#include "fds_basic.hpp"
#include "shared.hpp"
#include <iostream>
#include <vector>
using std::clog;
using std::endl;
using std::initializer_list;
using std::ostream;
using std::vector;

static inline void print_point(ostream &o, const vector<float> &nodes, u32 i) {
  o << (nodes[3 * i] - mesh.x0) / mesh.cell_size << " "
    << (nodes[3 * i + 1] - mesh.y0) / mesh.cell_size << " "
    << (nodes[3 * i + 2] - mesh.z0) / mesh.cell_size;
}

static inline void print_patch(ostream &o, const patch_info &patch) {
  o << "[" << patch.I1 << ", " << patch.I2 << "]"
    << "[" << patch.J1 << ", " << patch.J2 << "]"
    << "[" << patch.K1 << ", " << patch.K2 << "]";
}

static inline void print_debug(ostream &o, const vector<float> &nodes,
                               initializer_list<u32> is,
                               const patch_info &patch, const char *msg) {
  for (auto i : is) {
    o << "(";
    print_point(clog, nodes, i);
    o << ")";
  }
  o << " " << msg << " ";
  print_patch(clog, patch);
  o << endl;
}

static inline void print_points(ostream &o, const vector<float> &nodes, u32 i,
                                patch_info &patch) {}

static inline vector<u32> debug_primitive_on_boundary(
    const vector<patch_info> &patches, const vector<float> &nodes,
    const vector<u32> &primitive_indices, const bool wireframe) {
  assert(nodes.size() % 3 == 0);
  assert(!nodes.empty());
  assert(!patches.empty());
  vector<u32> res;

  auto _nodes = node_on_each_patch_boundary(patches, nodes);

  constexpr float r = .01f;

  // Loop for nodes on each boundary
  size_t i_patch = 0;
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

          print_debug(clog, nodes, {i0, i1}, patches[i_patch], "inside");
          break;
        } else {
          print_debug(clog, nodes, {i0, i1}, patches[i_patch], "outside");
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

          print_debug(clog, nodes, {i0, i1, i2}, patches[i_patch], "inside");
          break;
        } else {
          print_debug(clog, nodes, {i0, i1, i2}, patches[i_patch], "outside");
        }
      }
    ++i_patch;
  }
  return res;
}