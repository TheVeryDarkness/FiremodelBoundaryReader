#pragma once
#include "types.hpp"
#include <array>
#include <iomanip>
#include <iostream>

using std::array;
using std::cerr;
using std::char_traits;
using std::decay_t;
using std::endl;
using std::istream;
using std::numeric_limits;
using std::ostream;
using std::setw;
using std::tuple;
using std::ws;

template <size_t sz> static inline void read(istream &in, array<char, sz> &s) {
  CHECK_STREAM(in);
  static_assert(sz > 0);
  in.get(s.data(), sz);
  for (auto rp = s.rbegin() + 1, rend = s.rend(); rp != rend; ++rp)
    if (*rp == ' ')
      *rp = 0;
    else
      break;
}

// Actually I'm not sure about their meanings.
constexpr static inline const char string_separator[] = "\x1E\x00\x00\x00";
constexpr static inline const char integer_separator[] = "\x04\x00\x00\x00";
constexpr static inline const char line_separator[] = "\x24\x00\x00\x00";

template <size_t sz>
static inline void write_line(ostream &out, const array<char, sz> &s) {
  out.write(s.data(), sz) << endl;
}

template <typename Ty>
static inline ostream &write_number(ostream &out, const Ty &val,
                                    const char *sep = " ") {
  return out << setw(numeric_limits<Ty>::digits10 + 1) << val << sep;
}

static inline ostream &write_bytes(ostream &out, char c,
                                   const char *sep = " ") {
  return out << setw(sizeof(c)) << (unsigned)(unsigned char)(c) << sep;
}

template <size_t sz>
static inline void check(istream &in, const char (&s)[sz]) {
  CHECK_STREAM(in);
  static_assert(sz > 0);
  char buf[sz];
  buf[sz - 1] = 0;
  in.read(buf, sz - 1);
  if (strcmp(buf, s)) {
    cerr << "At pos " << in.tellg() << ", expect ";
    for (const char *p = s; p != s + sz; ++p)
      write_bytes(cerr, *p);
    cerr << ", but get ";
    for (const char *p = buf; p != buf + sz; ++p)
      write_bytes(cerr, *p);
    cerr << endl;
    std::terminate();
  }
}

template <typename Ty, size_t n>
static inline void write_separated(ostream &o, const char *sep,
                                   std::array<Ty, n> list) {
  auto begin = list.begin();
  auto end = list.end();
  if (begin == end)
    return;
  o << setw(sizeof(Ty) * 2) << *begin;
  ++begin;
  for (; begin != end; ++begin)
    o << sep << setw(sizeof(Ty) * 2) << *begin;
}

template <typename Ty> static inline Ty read_integer(istream &in) {
  CHECK_STREAM(in);
  constexpr auto size = sizeof(Ty);
  char buf[size];
  in.read(buf, size);
  // May cause error if endian not matched.
  Ty res = *reinterpret_cast<Ty *>(buf);

  // static_assert(std::is_unsigned_v<Ty>);
  // Ty res = 0;
  // for (const char* p = buf; p < buf + size; ++p)  // Big endian
  // for (const char *p = buf + size; p >= buf; --p) // Little endian
  //  (res <<= 8) += static_cast<unsigned char>(*p);
  return res;
}
static inline std::uint32_t read_uint32(istream &in) {
  return read_integer<std::uint32_t>(in);
}
static inline std::int32_t read_int32(istream &in) {
  return static_cast<std::int32_t>(read_integer<std::uint32_t>(in));
}

static inline std::float_t read_float(istream &in) {
  return read_integer<std::float_t>(in);
}

template <typename Ty> static inline void write_binary(ostream &o, Ty &&data) {
  o.write(reinterpret_cast<const char *>(&data), sizeof(data));
}
template <typename Ty>
static inline void write_vector_binary(ostream &o, const vector<Ty> &data) {
  o.write(reinterpret_cast<const char *>(data.data()),
          sizeof(Ty) * data.size());
}
