#pragma once
#include <iostream>
#include <sstream>

using std::string;

class stdin_proxy {
  std::stringstream sin;
  bool string_empty() {
    auto pos = sin.tellg();
    auto size = sin.str().size();
    bool res = pos == std::ios::pos_type(-1) || pos == size;
    return res;
  }

public:
  template <typename Any> stdin_proxy &operator>>(Any &&any) {
    sin >> std::ws;
    if (string_empty()) {
      std::cin >> std::forward<Any>(any);
    } else {
      sin >> std::forward<Any>(any);
    }
    return *this;
  }
  void endl() { sin << std::endl; }

  template <typename Any> stdin_proxy &operator<<(Any &&any) {
    sin.clear();
    sin << std::forward<Any>(any);
    return *this;
  }
  friend stdin_proxy &getline(stdin_proxy &in, string &s) {
    if (in.string_empty())
      std::getline(std::cin, s);
    else
      std::getline(in.sin, s);
    return in;
  }
  friend stdin_proxy &getline(stdin_proxy &in, string &s, char delim) {
    if (in.string_empty())
      std::getline(std::cin, s, delim);
    else
      std::getline(in.sin, s, delim);
    return in;
  }
} cin;
