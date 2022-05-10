#pragma once
#include <iostream>
#include <sstream>

using std::cerr;
using std::clog;
using std::enable_if_t;
using std::is_function_v;
using std::ostream;
using std::string;

#define COMMAND_NOT_FOUND cerr << "Command not found.\n"
#define FILE_OPEN_FAILED cerr << "File open failed.\n"
#define DIRECTORY_CREATE_FAILED cerr << "Directory create failed.\n"

static class stdin_proxy {
  std::stringstream sin;
  bool mode = false; // true for file in, false for stdin
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
    mode = !string_empty();
    return *this;
  }
  void endl() { sin << std::endl; }

  template <typename Any> stdin_proxy &operator<<(Any &&any) {
    sin.clear();
    sin << std::forward<Any>(any);
    mode = true;
    return *this;
  }
  friend stdin_proxy &getline(stdin_proxy &in, string &s) {
    if (in.string_empty())
      std::getline(std::cin, s);
    else {
      std::getline(in.sin, s);
    }
    in.mode = !in.string_empty();
    return in;
  }
  friend stdin_proxy &getline(stdin_proxy &in, string &s, char delim) {
    if (in.string_empty())
      std::getline(std::cin, s, delim);
    else
      std::getline(in.sin, s, delim);
    return in;
  }
  bool get_mode() { return mode; }
} cin;

static bool display_when_script_mode = false;
static class stdout_proxy {
  bool displayable() { return display_when_script_mode || !cin.get_mode(); }

public:
  template <typename Any /*, enable_if_t<!is_function_v<Any>, int> = 0*/>
  stdout_proxy &operator<<(Any &&any) {
    if (displayable())
      std::cout << std::forward<Any>(any);
    return *this;
  }
  stdout_proxy &operator<<(ostream &(*_Pfn)(ostream &)) {
    if (displayable())
      std::cout << _Pfn;
    return *this;
  }

  ostream &original() { return std::cout; }

} cout;
