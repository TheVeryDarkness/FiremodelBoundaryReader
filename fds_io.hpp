#pragma once

#include "fds_basic.hpp"
#include "io.hpp"
#include "proxy.hpp"

using std::ws;

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

static inline fds_boundary_file
read_frames(istream &fin, const vector<patch_info> patches, u32 frames_count) {
  fds_boundary_file frames;

  auto &data_map = frames.data;
  for (u32 i_patch = 0; i_patch < patches.size(); ++i_patch) {
    auto size = patches[i_patch].size();
    auto [iter, rep] = data_map.emplace(i_patch, patch_datas{size, {}});
    iter->second.data.reserve(frames_count * size);
  }
  while (fin.peek() != decay_t<decltype(fin)>::traits_type::eof() &&
         !fin.eof()) {
    check(fin, integer_separator);
    float stime = read_float(fin);
    check(fin, integer_separator);

    for (u32 ip = 0; ip < patches.size(); ++ip) {
      u32 patch_size = read_uint32(fin);
      const auto &info = patches[ip];
      auto _size = info.size();
      CHECK_FORMAT(_size * sizeof(float) == patch_size);
      auto &data = data_map[ip].data;
      for (size_t i = 0; i < info.size(); ++i) {
        float val = read_float(fin);
        data.push_back(val);
      }
      u32 patch_end = read_uint32(fin);
      CHECK_FORMAT(patch_end == patch_size);
    }

    frames.times.push_back(stime);
  }
  return frames;
}

static inline fds_boundary_file
read_frames_of_specified_patches(istream &fin, const vector<patch_info> patches,
                                 const vector<u32> &patch_indices,
                                 u32 frames_count) {
  fds_boundary_file frames;

  auto &data_map = frames.data;
  for (u32 i_patch : patch_indices) {
    assert(i_patch < patches.size());
    auto sz = patches[i_patch].size();
    auto [iter, rep] = data_map.emplace(i_patch, patch_datas{sz, {}});
    iter->second.data.reserve(sz * frames_count);
  }
  // for (u32 i_patch = 0; i_patch < patches.size(); ++i_patch)
  while (fin.peek() != decay_t<decltype(fin)>::traits_type::eof() &&
         !fin.eof()) {
    check(fin, integer_separator);
    const float stime = read_float(fin);
    check(fin, integer_separator);

    for (u32 ip = 0; ip < patches.size(); ++ip) {
      u32 patch_size = read_uint32(fin);
      const auto &info = patches[ip];
      const auto _size = info.size();
      CHECK_FORMAT(_size * sizeof(float) == patch_size);
      if (data_map.find(ip) != data_map.cend()) {
        auto &data = data_map[ip].data;
        for (size_t i = 0; i < _size; ++i) {
          float val = read_float(fin);
          data.push_back(val);
        }
      } else {
        fin.ignore(patch_size);
      }
      u32 patch_end = read_uint32(fin);
      CHECK_FORMAT(patch_end == patch_size);
    }

    frames.times.push_back(stime);
  }
  return frames;
}

static inline u32 calculate_frames_count(const vector<patch_info> &patches,
                                         u32 file_size) {
  u32 frame_bytes = 3;
  for (auto &p : patches)
    frame_bytes += p.size() + 2;
  frame_bytes *= 4;
  u32 header_bytes = (4U + 30U + 4U) * 3U + (3U + 11U * patches.size()) * 4U;
  u32 remained_bytes = (file_size - header_bytes);
  assert(remained_bytes % frame_bytes == 0);
  u32 frames_count = remained_bytes / frame_bytes;
  return frames_count;
}

/*
 * @return Label, Bar Label, Units, Patches
 */
static inline tuple<array<char, 30 + 1>, array<char, 30 + 1>,
                    array<char, 30 + 1>, vector<patch_info>>
read_file_header_and_patches(istream &fin) {
  auto &&header = read_file_header(fin);
  auto &&patches = read_patches(fin);
  return std::tuple_cat(std::move(header), std::make_tuple(std::move(patches)));
}

/*
 * @return Label, Bar Label, Units, Patches, Frames
 */
static inline tuple<array<char, 30 + 1>, array<char, 30 + 1>,
                    array<char, 30 + 1>, vector<patch_info>, fds_boundary_file>
read_file(istream &fin, u32 file_size) {
  auto &&header = read_file_header(fin);
  auto &&patches = read_patches(fin);

  const auto frames_count = calculate_frames_count(patches, file_size);

  auto &&frames = read_frames(fin, patches, frames_count);
  return std::tuple_cat(std::move(header),
                        std::make_tuple(std::move(patches), std::move(frames)));
}

/*
 * @return Label, Bar Label, Units, Patches, Frames
 */
static inline tuple<array<char, 30 + 1>, array<char, 30 + 1>,
                    array<char, 30 + 1>, vector<patch_info>, fds_boundary_file>
read_file_of_specified_patches(istream &fin, u32 file_size) {
  auto &&header = read_file_header(fin);
  auto &&patches = read_patches(fin);

  vector<u32> indices = read_until<u32>(cin);

  for (auto index : indices)
    if (index >= patches.size())
      clog << index << " is out of bound.\n";

  const auto frames_count = calculate_frames_count(patches, file_size);

  auto &&frames =
      read_frames_of_specified_patches(fin, patches, indices, frames_count);
  return std::tuple_cat(std::move(header),
                        std::make_tuple(std::move(patches), std::move(frames)));
}

static inline ostream &operator<<(ostream &o, const patch_info &patch) {
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

static inline void print_header(ostream &o, const array<char, 30 + 1> &label,
                                const array<char, 30 + 1> &bar_label,
                                const array<char, 30 + 1> &units) {
  write_line(o << "Label:     ", label);
  write_line(o << "Bar Label: ", bar_label);
  write_line(o << "Units:     ", units);
}

static inline void print_patches(ostream &o, const vector<patch_info> patches) {
  constexpr auto int_len = numeric_limits<i32>::digits10 + 1;
  constexpr auto uint_len = numeric_limits<u32>::digits10 + 1;
  o << setw(uint_len) << "INDEX" << setw(uint_len) << "I1"
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
  u32 i_patch = 0;
  for (const auto &patch : patches) {
    o << setw(uint_len) << i_patch << patch << endl;
    ++i_patch;
  }
}

static inline void print_frames(ostream &o, const fds_boundary_file frames) {
  size_t i = 0;
  for (const auto &f : frames.times) {
    o << "Frame " << i << " at " << f << "s." << endl;
    ++i;
  }
}