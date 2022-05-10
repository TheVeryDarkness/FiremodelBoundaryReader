#pragma once
#include "proxy.hpp"
#include <filesystem>
#include <iostream>
#include <optional>
using std::endl;
using std::optional;
using std::vector;
using std::filesystem::directory_entry;
using std::filesystem::directory_iterator;
using std::filesystem::path;

template <typename Pred>
static inline optional<directory_entry>
request_file_by_id(directory_iterator di, Pred &&pred, const char *type) {
  size_t i = 0;
  vector<directory_entry> entries;
  for (const auto &entry : di)
    if (pred(entry)) {
      cout << i << " - " << entry.path().generic_string() << endl;
      entries.push_back(entry);
      ++i;
    }
  if (entries.empty()) {
    clog << "No " << type << " found." << endl;
    return {};
  }
  cout << "There are " << entries.size() << " " << type << " in total." << endl;
  cout << "Index(Invalid index to discard): ";
  cin >> i;
  if (i < entries.size())
    return entries[i];
  clog << "Not a valid index." << endl;
  return {};
}

template <typename Pred>
static inline optional<path> request_file_by_name(Pred &&pred,
                                                  const char *type) {
  cout << "Path to " << type << " file: ";
  string s;
  getline(cin, s);
  if (s.empty())
    getline(cin, s);
  path p = s;
  if (pred(p))
    return p;
  clog << "Not " << type << "." << endl;
  return {};
}