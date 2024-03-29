#pragma once
#include "apdl_io.hpp"
#include "fs.hpp"
#include "proxy.hpp"
#include "types.hpp"
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
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
using std::map;
using std::move;
using std::ofstream;
using std::string;
using std::tuple;
using std::vector;
using std::ws;
using std::filesystem::path;

/// @brief Read nodes data from stream associated with standard file.
/// @param in: Stream associated with standard file.
/// @return Nodes coordinates.
static inline vector<float> read_standard_nodes(istream &in) {
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
static inline tuple<vector<u32>, vector<u32>>
read_standard_elements(istream &in) {
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

/// @brief Read node numbers from stream associated with standard file.
/// @param in: Stream associated with standard file.
/// @return Nodes numbers.
static inline vector<u32> read_standard_node_numbers(istream &in) {
  vector<u32> res;
  while (in.peek() != char_traits<char>::eof() && in) {
    u32 node_number = 0;
    in >> node_number;
    res.push_back(node_number);
    in >> ws;
  }
  return res;
}

/// @brief Read element numbers from stream associated with standard file.
/// @param in: Stream associated with standard file.
/// @return Element numbers.
static inline vector<u32> read_standard_element_numbers(istream &in) {
  vector<u32> res;
  while (in.peek() != char_traits<char>::eof() && in) {
    u32 element_number = 0;
    in >> element_number;
    res.push_back(element_number);
    in >> ws;
  }
  return res;
}

static inline void write_standard_nodes(ostream &out,
                                        const vector<float> &nodes) {
  assert(nodes.size() % 3 == 0);
  out << nodes.size() / 3 << endl;
  size_t i = 0;
  for (auto n : nodes) {
    out << n;

    ++i;
    if (i == 3) {
      out << endl;
      i = 0;
    } else
      out << ' ';
  }
}

static inline void write_standard_elements(ostream &out,
                                           const vector<u32> &size,
                                           const vector<u32> &indices) {
  map<u32, vector<u32>> m;
  {
    size_t i = 0;
    for (auto sz : size) {
      auto iter = m.find(sz);
      if (iter == m.cend()) {
        m.emplace(sz, vector<u32>{});
      }

      auto &vec = m.at(sz);
      vec.insert(vec.cend(), indices.cbegin() + i, indices.cbegin() + i + sz);
      i += sz;
    }
  }

  for (auto &[sz, vec] : m) {
    assert(vec.size() % sz == 0);
    out << vec.size() / sz << ' ' << sz << endl;

    size_t i = 0;
    for (auto index : vec) {
      out << index;
      ++i;
      if (i == sz) {
        i = 0;
        out << endl;
      } else
        out << ' ';
    }
  }
}

static inline void write_standard_element_numbers(ostream &out,
                                                  const vector<u32> &numbers) {
  for (auto &number : numbers) {
    out << number << '\n';
  }
}

static inline void write_standard_node_numbers(ostream &out,
                                               const vector<u32> &numbers) {
  for (auto &number : numbers) {
    out << number << '\n';
  }
}

/// @brief Read nodes and elements data from user-selected files.
/// @return Tuple that contains nodes coordinates, vertices counts, vertex
/// indices.
static inline vector<float> read_nodes() {
  char opt = 'd';
  cout << R"(
s - Standard format.
a - ANSYS mapdl file.
)";
  cin >> opt;
  switch (opt) {
  case 's': {
    auto opt_nodes =
        request_file_by_name([](const path &p) { return exists(p); }, "nodes");
    if (!opt_nodes)
      return {};
    auto in_nodes = ifstream(opt_nodes.value());
    auto nodes = read_standard_nodes(in_nodes);
    return std::move(nodes);
  } break;
  case 'a': {
    auto opt_apdl =
        request_file_by_name([](const path &p) { return exists(p); }, "APDL");
    if (!opt_apdl)
      return {};
    auto in = ifstream(opt_apdl.value());
    auto [nodes, _0, _1, _2] = read_mapdl<true, false, false, false, true>(in);
    return std::move(nodes);
  } break;
  default:
    COMMAND_NOT_FOUND;
  case 'd':
    return {};
  }
}

/// @brief Read nodes and elements data from user-selected files.
/// @return Nodes coordinates, element node counts, vertex indices, node
/// numbers, element numbers.
static inline tuple<vector<float>, vector<u32>, vector<u32>, vector<u32>,
                    vector<u32>>
read_nodes_and_elements() {
  char opt = 'd';
  cout << R"(
s - Standard format.
a - ANSYS mapdl file.
)";
  cin >> opt;
  switch (opt) {
  case 's': {
    auto opt_nodes =
        request_file_by_name([](const path &p) { return exists(p); }, "nodes");
    auto opt_elems = request_file_by_name(
        [](const path &p) { return exists(p); }, "elements");
    auto opt_node_nums = request_file_by_name(
        [](const path &p) { return exists(p); }, "node numbers");
    auto opt_elem_nums = request_file_by_name(
        [](const path &p) { return exists(p); }, "element numbers");
    if (!opt_nodes || !opt_elems)
      return {};
    if (!opt_node_nums || !opt_elem_nums)
      clog << "Element numbers may not read.\n";
    auto in_nodes = ifstream(opt_nodes.value());
    auto nodes = read_standard_nodes(in_nodes);

    auto in_elems = ifstream(opt_elems.value());
    auto [elems, sizes] = read_standard_elements(in_elems);

    vector<u32> node_numbers;
    if (opt_node_nums) {
      auto in_nums = ifstream(opt_node_nums.value());
      node_numbers = read_standard_node_numbers(in_nums);
    }

    vector<u32> elem_numbers;
    if (opt_elem_nums) {
      auto in_nums = ifstream(opt_elem_nums.value());
      elem_numbers = read_standard_element_numbers(in_nums);
    }
    return {move(nodes), move(sizes), move(elems), move(node_numbers),
            move(elem_numbers)};
  }
  case 'a': {
    auto opt_apdl =
        request_file_by_name([](const path &p) { return exists(p); }, "APDL");
    if (!opt_apdl)
      return {};
    auto in = ifstream(opt_apdl.value());
    auto &&[nodes, size_and_elements, node_numbers, element_numbers] =
        read_mapdl<true, true, true, true, true>(in);
    auto &&[size, elemements] = std::move(size_and_elements);

    cout
        << "Whether to save nodes and elements data in standard format?(y or n)"
        << endl;
    char opt = 'n';
    cin >> opt;
    if (opt == 'y') {
      auto opt_nodes =
          request_file_by_name([](const path &p) { return true; }, "nodes");
      auto opt_elems =
          request_file_by_name([](const path &p) { return true; }, "elements");
      auto opt_node_nums = request_file_by_name(
          [](const path &p) { return true; }, "node numbers");
      auto opt_elem_nums = request_file_by_name(
          [](const path &p) { return true; }, "element numbers");

      if (!opt_nodes || !opt_elems || !opt_node_nums || !opt_elem_nums)
        return {};

      {
        auto out_nodes = ofstream(opt_nodes.value());
        write_standard_nodes(out_nodes, nodes);
        out_nodes.close();
      }

      {
        auto out_elems = ofstream(opt_elems.value());
        write_standard_elements(out_elems, size, elemements);
        out_elems.close();
      }

      {
        auto out_node_nums = ofstream(opt_node_nums.value());
        write_standard_node_numbers(out_node_nums, node_numbers);
        out_node_nums.close();
      }

      {
        auto out_elem_nums = ofstream(opt_elem_nums.value());
        write_standard_element_numbers(out_elem_nums, element_numbers);
        out_elem_nums.close();
      }
      clog << "Done.\n";
    }
    return {move(nodes), move(size), move(elemements), move(node_numbers),
            move(element_numbers)};
  }
  default:
    COMMAND_NOT_FOUND;
  case 'd':
    return {};
  }
}