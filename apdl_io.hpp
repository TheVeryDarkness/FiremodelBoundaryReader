#pragma once

#include "fds_basic.hpp"
#include "io.hpp"
#include <cassert>
#include <filesystem>
#include <fstream>
#include <string>

using std::clog;
using std::conditional_t;
using std::getline;
using std::ofstream;
using std::string;
using std::string_view;
using std::to_string;
using std::filesystem::absolute;
using std::filesystem::path;

template <size_t ParameterCount, bool withSize, size_t SizeAt, size_t NumberAt>
static inline tuple<vector<u32>, u32, u32> read_element_in_apdl(istream &in,
                                                                u32 ei) {
  static_assert(SizeAt < ParameterCount);
  u32 unknown[ParameterCount] = {};
  for (auto &u : unknown)
    in >> u;

  if constexpr (withSize) {
    auto size = unknown[SizeAt];

    vector<u32> indices = {};
    indices.resize(size);

    for (size_t i = 0; i < size; ++i) {
      auto &n = indices.at(i);
      in >> n;
      if (n == 0) {
        cerr << "Element " << ei << " has a node 0." << endl;
      }
      n--;
    }
    return {std::move(indices), size, unknown[NumberAt]};
  } else {
    vector<u32> indices = {};
    while (true) {
      while (in.peek() == ' ') {
        in.get();
      }
      if (in.peek() == '\n') {
        break;
      }
      u32 n = 0;
      in >> n;
      if (n == 0) {
        cerr << "Element " << ei << " has a zero node." << endl;
      }
      indices.push_back(n - 1);
    }
    assert(indices.size() <= numeric_limits<u32>::max());
    const u32 size = (u32)indices.size();
    return {std::move(indices), size, unknown[NumberAt]};
  }
}

/// @brief Read from mapdl file.
/// @param in: Stream that reads from mapdl file.
/// @return Nodes coordinates, vertex count, element vertex indices.
template <bool readNodes, bool readElements, bool readElementNumbers>
static inline tuple<conditional_t<readNodes, vector<float>, tuple<>>,
                    tuple<conditional_t<readElements, vector<u32>, tuple<>>,
                          conditional_t<readElements, vector<u32>, tuple<>>>,
                    conditional_t<readElementNumbers, vector<u32>, tuple<>>>
read_mapdl(istream &in) {
  tuple<conditional_t<readNodes, vector<float>, tuple<>>,
        tuple<conditional_t<readElements, vector<u32>, tuple<>>,
              conditional_t<readElements, vector<u32>, tuple<>>>,
        conditional_t<readElementNumbers, vector<u32>, tuple<>>>
      res;
  auto &[nodes, size_and_element, numbers] = res;
  auto &[vertex_count, elements] = size_and_element;

  while (in.peek() != char_traits<char>::eof() && in) {
    string trash;
    string line;
    getline(in, line);
    string_view l = line;
    auto comma1 = l.find_first_of(",");
    if (comma1 == string_view::npos)
      continue;

    auto first = l.substr(0, comma1);
    // clog << first << endl;
    if constexpr (readNodes)
      if (first == "nblock") {
        getline(in, trash);

        while (in) {
          in >> ws;
          if (in.peek() == '-')
            break;

          u32 index;
          float coordinates[3];
          in >> index;

          // Assume continuous
          assert(index == nodes.size() / 3 + 1);

          for (auto &c : coordinates)
            in >> c;
          coordinates[2] = -coordinates[2];
          for (auto &c : coordinates)
            nodes.push_back(c);
        }
      }
    if constexpr (readElements)
      if (first == "eblock") {
        getline(in, trash);

        auto comma2 = line.find_first_of(',', comma1 + 1);
        auto second = line.substr(comma1 + 1, comma2 - comma1 - 1);
        auto comma3 = line.find_first_of(',', comma2 + 1);
        auto third = line.substr(comma2 + 1, comma3 - comma2 - 1);

        while (in) {
          in >> ws;
          if (in.peek() == '-')
            break;

          assert(vertex_count.size() <= numeric_limits<u32>::max());
          if (third == "solid") {
            auto [indices, size, number] =
                read_element_in_apdl<11, true, 8, 10>(in,
                                                      (u32)vertex_count.size());
            for (auto n : indices)
              elements.push_back(n);
            vertex_count.push_back(size);
            if constexpr (readElementNumbers)
              numbers.push_back(number);
          } else {
            auto [indices, size, number] = read_element_in_apdl<5, false, 0, 0>(
                in, (u32)vertex_count.size());
            for (auto n : indices)
              elements.push_back(n);
            vertex_count.push_back(size);
            if constexpr (readElementNumbers)
              numbers.push_back(number);
          }
        }
      }

    in >> ws;
  }
  if constexpr (readNodes)
    if (nodes.empty())
      cerr << "Nodes definition not found" << endl;
  if constexpr (readElements)
    if (elements.empty())
      cerr << "Elements definition not found" << endl;
  return res;
}

static inline void write_table(
    ostream &o, const path &directory, const char *name, u32 element_number,
    u32 surface_number, vector<u32>::const_iterator surface_indices_begin,
    vector<u32>::const_iterator surface_indices_end,
    const vector<frame> &frames, vector<float>::const_iterator vec_begin,
    vector<float>::const_iterator vec_end) {

  static const path ext = path("txt");
  assert(vec_begin + frames.size() == vec_end);

  if (all_of(vec_begin, vec_end, [](float data) { return data == 0.0f; }))
    return;

  auto p = directory /
           path(to_string(element_number) + '_' + to_string(surface_number))
               .replace_extension(ext);
  ofstream tout;
  tout.open(p);
  assert(tout);
  tout << "TIME " << name << '\n';

  for (const frame &frame : frames) {
    assert(vec_begin != vec_end);
    tout << frame.time << ' ' << *vec_begin * 1000 << '\n';
    ++vec_begin;
  }
  tout.close();

#pragma omp critical
  {
    o << "*DIM," << name << element_number << '_' << surface_number << ",TABLE,"
      << frames.size() << ",1,1,TIME,\n";
    o << "*TREAD," << name << element_number << '_' << surface_number << ','
      << element_number << '_' << surface_number << ",txt,"
      << absolute(directory).generic_string() << ",1\n";
    o << "ESEL,S,,," << element_number << "\n";
    o << "SFE,ALL," << surface_number << ',' << name << ",%" << name
      << element_number << '_' << surface_number << "%\n";
  }
}