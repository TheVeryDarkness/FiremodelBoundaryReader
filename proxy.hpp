#pragma once
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

using std::cerr;
using std::clog;
using std::enable_if_t;
using std::is_function_v;
using std::istream;
using std::ostream;
using std::string;
using std::vector;

#define COMMAND_NOT_FOUND_OPT(opt) cerr << "Option '" << opt << "' not found.\n"
#define COMMAND_NOT_FOUND cerr << "Option '" << opt << "' not found.\n"
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
  void check() {
    /*
    if (sin.fail())
      cerr << "Error detected during executing script.\n";
    mode = !string_empty();
  */
  }

public:
  template <typename Any> stdin_proxy &operator>>(Any &&any) {
    check();

    sin >> std::ws;
    if (string_empty()) {
      std::cin >> std::forward<Any>(any);
    } else {
      sin >> std::forward<Any>(any);
      sin >> std::ws;
    }
    mode = !string_empty();
    return *this;
  }
  void endl() {
    check();
    sin << std::endl;
  }
  void ws() {
    check();
    sin >> std::ws;
    if (string_empty())
      std::cin >> std::ws;
  }
  auto peek() { return string_empty() ? std::cin.peek() : sin.peek(); }
  auto &ignore() { return string_empty() ? std::cin.ignore() : sin.ignore(); }

  operator bool() { return !string_empty() || std::cin; }

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
      in.sin >> std::ws;
    }
    in.mode = !in.string_empty();
    return in;
  }
  friend stdin_proxy &getline(stdin_proxy &in, string &s, char delim) {
    if (in.string_empty())
      std::getline(std::cin, s, delim);
    else {
      std::getline(in.sin, s, delim);
      in.sin >> std::ws;
    }
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

template <typename Ty>
static inline vector<Ty> read_until(stdin_proxy &in, char delim = '\\') {
  vector<Ty> res;
  while (in) {
    in.ws();
    if (in.peek() == delim) {
      in.ignore();
      break;
    } else {
      res.push_back({});
      in >> res.back();
    }
  }
  return res;
}