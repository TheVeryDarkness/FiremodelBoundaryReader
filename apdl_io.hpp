#pragma once

#include "io.hpp"
#include <cassert>
#include <string>

using std::clog;
using std::getline;
using std::string;
using std::string_view;

template <size_t ParameterCount, bool withSize, size_t SizeAt>
static inline tuple<vector<u32>, u32> read_element(istream &in, u32 ei) {
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
    return {std::move(indices), size};
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
    return {std::move(indices), size};
  }
}

/// @brief Read from mapdl file.
/// @param in: Stream that reads from mapdl file.
/// @return Nodes coordinates, vertex count, element vertex indices.
static inline tuple<vector<float>, vector<u32>, vector<u32>>
read_mapdl(istream &in) {
  tuple<vector<float>, vector<u32>, vector<u32>> res;
  auto &[nodes, vertex_count, elements] = res;

  while (in.peek() != char_traits<char>::eof() && in) {
    string trash;
    string line;
    getline(in, line);
    string_view l = line;
    auto comma1 = l.find_first_of(",");
    if (comma1 == string_view::npos)
      continue;

    auto first = l.substr(0, comma1);
    clog << first << endl;
    if (first == "nblock") {
      getline(in, trash);

      // Be cautious not to use first below.
      auto &[nodes, vertex_count, elements] = res;
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
        for (auto &c : coordinates)
          nodes.push_back(c);
      }
    } else if (first == "eblock") {
      getline(in, trash);

      auto comma2 = line.find_first_of(',', comma1 + 1);
      auto second = line.substr(comma1 + 1, comma2 - comma1 - 1);
      auto comma3 = line.find_first_of(',', comma2 + 1);
      auto third = line.substr(comma2 + 1, comma3 - comma2 - 1);

      // Be cautious not to use first below.
      auto &[nodes, vertex_count, elements] = res;
      while (in) {
        in >> ws;
        if (in.peek() == '-')
          break;

        assert(vertex_count.size() <= numeric_limits<u32>::max());
        if (third == "solid") {
          auto [indices, size] =
              read_element<11, true, 8>(in, (u32)vertex_count.size());
          for (auto n : indices)
            elements.push_back(n);
          vertex_count.push_back(size);
        } else {
          auto [indices, size] =
              read_element<5, false, 0>(in, (u32)vertex_count.size());
          for (auto n : indices)
            elements.push_back(n);
          vertex_count.push_back(size);
        }
      }
    }

    in >> ws;
  }
  if (nodes.empty())
    cerr << "Nodes definition not found" << endl;
  if (elements.empty())
    cerr << "Elements definition not found" << endl;
  return res;
}