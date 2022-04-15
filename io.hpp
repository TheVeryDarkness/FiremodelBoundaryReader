#pragma once
#include "types.hpp"
#include <array>
#include <iomanip>
#include <iostream>

using std::array;
using std::cerr;
using std::decay_t;
using std::endl;
using std::istream;
using std::numeric_limits;
using std::ostream;
using std::setw;
using std::tuple;

template <size_t sz> static inline void read(istream &in, array<char, sz> &s) {
  CHECK_STREAM(in);
  static_assert(sz > 0);
  in.get(s.data(), sz);
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

std::float_t read_float(istream &in) { return read_integer<std::float_t>(in); }

static inline tuple<array<char, 30 + 1>, array<char, 30 + 1>,
                    array<char, 30 + 1>>
read_file_header(istream &fin) {
  tuple<array<char, 30 + 1>, array<char, 30 + 1>, array<char, 30 + 1>> res;
  auto &[label, bar_label, units] = res;

  check(fin, string_separator);
  read(fin, label);
  check(fin, string_separator);
  check(fin, string_separator);
  read(fin, bar_label);
  check(fin, string_separator);
  check(fin, string_separator);
  read(fin, units);
  check(fin, string_separator);

  return res;
}

static inline vector<patch_info> read_patches(istream &fin) {
  u32 n_patch;

  check(fin, integer_separator);
  n_patch = read_uint32(fin);
  check(fin, integer_separator);

  vector<patch_info> patches;
  patches.reserve(n_patch);
  for (u32 i = 0; i < n_patch; ++i) {
    check(fin, line_separator);
    u32 I1 = read_uint32(fin);
    u32 I2 = read_uint32(fin);
    u32 J1 = read_uint32(fin);
    u32 J2 = read_uint32(fin);
    u32 K1 = read_uint32(fin);
    u32 K2 = read_uint32(fin);
    i32 IOR = read_int32(fin);
    u32 OBST_INDEX = read_uint32(fin);
    u32 NM = read_uint32(fin);
    patches.push_back({I1, I2, J1, J2, K1, K2, IOR, OBST_INDEX, NM});
    check(fin, line_separator);
  }
  return patches;
}

static inline vector<frame> read_frames(istream &fin,
                                        const vector<patch_info> patches) {
  vector<frame> frames;
  while (fin.peek() != decay_t<decltype(fin)>::traits_type::eof() &&
         !fin.eof()) {
    vector<patch_data> current(patches.size(), patch_data{});
    check(fin, integer_separator);
    float stime = read_float(fin);
    check(fin, integer_separator);

    for (size_t ip = 0; ip < patches.size(); ++ip) {
      u32 patch_size = read_uint32(fin);
      const auto &info = patches[ip];
      auto _size = info.size();
      CHECK_FORMAT(_size * sizeof(float) == patch_size);
      auto &data = current[ip].data;
      data.reserve(_size);
      for (size_t i = 0; i < info.size(); ++i) {
        float val = read_float(fin);
        data.push_back(val);
      }
      u32 patch_end = read_uint32(fin);
      CHECK_FORMAT(patch_end == patch_size);
    }

    frames.push_back(frame{stime, std::move(current)});
  }
  return frames;
}

/*
 * @retval Label, Bar Label, Units, Patches
 */
static inline tuple<array<char, 30 + 1>, array<char, 30 + 1>,
                    array<char, 30 + 1>, vector<patch_info>>
read_file_header_and_patches(istream &fin) {
  auto &&header = read_file_header(fin);
  auto &&patches = read_patches(fin);
  return std::tuple_cat(std::move(header), std::make_tuple(std::move(patches)));
}

/*
 * @retval Label, Bar Label, Units, Patches, Frames
 */
static inline tuple<array<char, 30 + 1>, array<char, 30 + 1>,
                    array<char, 30 + 1>, vector<patch_info>, vector<frame>>
read_file(istream &fin) {
  auto &&header = read_file_header(fin);
  auto &&patches = read_patches(fin);
  auto &&frames = read_frames(fin, patches);
  return std::tuple_cat(std::move(header),
                        std::make_tuple(std::move(patches), std::move(frames)));
}

ostream &operator<<(ostream &o, const patch_info &patch) {
  write_number(o, patch.I1);
  write_number(o, patch.I2);
  write_number(o, patch.J1);
  write_number(o, patch.J2);
  write_number(o, patch.K1);
  write_number(o, patch.K2);
  o << " " << patch.IOR_repr() << " ";
  write_number(o, patch.OBST_INDEX);
  write_number(o, patch.NM);
  write_number(o, patch.size());
  return o;
}

void print_header(ostream &o, const array<char, 30 + 1> &label,
                  const array<char, 30 + 1> &bar_label,
                  const array<char, 30 + 1> &units) {
  write_line(o << "Label:     ", label);
  write_line(o << "Bar Label: ", bar_label);
  write_line(o << "Units:     ", units);
}

void print_patches(ostream &o, const vector<patch_info> patches) {
  constexpr auto int_len = numeric_limits<i32>::digits10 + 1;
  constexpr auto uint_len = numeric_limits<u32>::digits10 + 1;
  o << setw(uint_len) << "I1"
    << " " << setw(uint_len) << "I2"
    << " " << setw(uint_len) << "J1"
    << " " << setw(uint_len) << "J2"
    << " " << setw(uint_len) << "K1"
    << " " << setw(uint_len) << "K2"
    << " " << setw(3) << "IOR"
    << " " << setw(uint_len) << "OBST_INDEX"
    << " " << setw(uint_len) << "MESH_INDEX"
    << " " << setw(uint_len) << "SIZE"
    << " " << endl;
  for (const auto &patch : patches)
    o << patch << endl;
}

void print_frames(ostream &o, const vector<frame> frames) {
  size_t i = 0;
  for (const auto &f : frames) {
    o << "Frame " << i << " at " << f.time << "s." << endl;
    ++i;
  }
}

template <typename Ty> void write_binary(ostream &o, Ty &&data) {
  o.write(reinterpret_cast<const char *>(&data), sizeof(data));
}
template <typename Ty>
void write_vector_binary(ostream &o, const vector<Ty> &data) {
  o.write(reinterpret_cast<const char *>(data.data()),
          sizeof(Ty) * data.size());
}