#pragma once
#include "types.hpp"
#include <tuple>

using std::conditional_t;
using std::lexicographical_compare;
using std::tuple;

struct patch_info {
  const u32 I1;
  const u32 I2;
  const u32 J1;
  const u32 J2;
  const u32 K1;
  const u32 K2;
  const i32 IOR;        // Orientation. 123 stand for XYZ.
  const u32 OBST_INDEX; // Obstruction index
  const u32 NM;         // Mesh index
  patch_info(u32 I1, u32 I2, u32 J1, u32 J2, u32 K1, u32 K2, i32 IOR,
             u32 OBST_INDEX, u32 NM)
      : I1(I1), I2(I2), J1(J1), J2(J2), K1(K1), K2(K2), IOR(IOR),
        OBST_INDEX(OBST_INDEX), NM(NM) {
    CHECK_FORMAT(I2 >= I1);
    CHECK_FORMAT(J2 >= J1);
    CHECK_FORMAT(K2 >= K1);
    CHECK_FORMAT(IOR != 0);
    CHECK_FORMAT(IOR <= 3);
    CHECK_FORMAT(IOR >= -3);
  }
  patch_info(const patch_info &) noexcept = default;
  template <typename Callable> void for_each(Callable &&callable) const {
    callable(I1);
    callable(I2);
    callable(J1);
    callable(J2);
    callable(K1);
    callable(K2);
    callable(IOR);
    callable(OBST_INDEX);
    callable(NM);
  }

  const char (&IOR_repr() const noexcept)[3] {
    switch (IOR) {
    case 1:
      return "+X";
    case -1:
      return "-X";
    case 2:
      return "+Y";
    case -2:
      return "-Y";
    case 3:
      return "+Z";
    case -3:
      return "-Z";
    default:
      return "??";
    }
  }

  u32 I() const noexcept { return I2 - I1 + 1; }
  u32 J() const noexcept { return J2 - J1 + 1; }
  u32 K() const noexcept { return K2 - K1 + 1; }

  u32 size() const noexcept { return I() * J() * K(); }
};
struct patch_data {
  vector<float> data;
};
struct frame {
  float time;
  vector<patch_data> data;
};

enum class data_category { temperature, other };

template <size_t len>
constexpr static inline bool compare(const array<char, 30 + 1> &a,
                                     const char (&b)[len]) {
  static_assert(len < 30);
  return lexicographical_compare(a.data(), a.data() + len + 1, b,
                                 b + len + 1) == 0;
}

constexpr static inline data_category
get_data_category(const array<char, 30 + 1> &units) {
  if (compare(units, "C")) {
    return data_category::temperature;
  } else
    return data_category::other;
}

static struct {
  float cell_size = 1;
  float x0 = 0;
  float y0 = 0;
  float z0 = 0;
} mesh;

template <bool withPatchIndices>
static inline tuple<
    tuple<vector<u32>, conditional_t<withPatchIndices, vector<u32>, tuple<>>>,
    vector<u32>>
from_patches(const vector<patch_info> &patches, bool wireframe) {
  tuple<
      tuple<vector<u32>, conditional_t<withPatchIndices, vector<u32>, tuple<>>>,
      vector<u32>>
      res;
  auto &[vertices, indices] = res;
  auto &[position, patchIndices] = vertices;

  position.reserve(patches.size() * 12);
  if constexpr (withPatchIndices)
    patchIndices.reserve(patches.size() * 4);
  indices.reserve(wireframe ? indices.size() * 8 : indices.size() * 6);

  u32 i_patch = 0;
  for (const auto &patch : patches) {
    const u32 N = (u32)position.size() / 3;
    if (wireframe) {
      switch (patch.IOR) {
      case 1:
      case -1:
        indices.push_back(N + 0);
        indices.push_back(N + 1);
        indices.push_back(N + 1);
        indices.push_back(N + 3);
        indices.push_back(N + 3);
        indices.push_back(N + 2);
        indices.push_back(N + 2);
        indices.push_back(N + 0);
        break;
      case 2:
      case -2:
        indices.push_back(N + 0);
        indices.push_back(N + 2);
        indices.push_back(N + 2);
        indices.push_back(N + 1);
        indices.push_back(N + 1);
        indices.push_back(N + 3);
        indices.push_back(N + 3);
        indices.push_back(N + 0);
        break;
      case 3:
      case -3:
        indices.push_back(N + 0);
        indices.push_back(N + 1);
        indices.push_back(N + 1);
        indices.push_back(N + 2);
        indices.push_back(N + 2);
        indices.push_back(N + 3);
        indices.push_back(N + 3);
        indices.push_back(N + 0);
        break;
      default:
        break;
      }
    } else {
      indices.push_back(N + 0);
      indices.push_back(N + 1);
      indices.push_back(N + 2);
      switch (patch.IOR) {
      case 1:
      case -1:
        indices.push_back(N + 1);
        indices.push_back(N + 2);
        indices.push_back(N + 3);
        break;
      case 2:
      case -2:
        indices.push_back(N + 0);
        indices.push_back(N + 1);
        indices.push_back(N + 3);
        break;
      case 3:
      case -3:
        indices.push_back(N + 0);
        indices.push_back(N + 2);
        indices.push_back(N + 3);
        break;
      default:
        break;
      }
    }

    position.push_back(patch.I1);
    position.push_back(patch.J1);
    position.push_back(patch.K1);

    position.push_back(patch.I2);
    position.push_back(patch.J1);
    position.push_back(patch.K2);

    position.push_back(patch.I2);
    position.push_back(patch.J2);
    position.push_back(patch.K1);

    position.push_back(patch.I1);
    position.push_back(patch.J2);
    position.push_back(patch.K2);

    if constexpr (withPatchIndices)
      for (size_t i = 0; i < 4; ++i)
        patchIndices.push_back(i_patch);

    ++i_patch;
  }

  return res;
}

static inline tuple<tuple<vector<u32>, vector<float>>, vector<u32>>
from_data(const vector<patch_info> &patches, const frame &frame) {
  tuple<tuple<vector<u32>, vector<float>>, vector<u32>> res;
  auto &[vertices, indices] = res;
  auto &[position, _data] = vertices;

  u32 points = 0;
  for (auto &p : patches) {
    points += p.size();
  }

  _data.reserve(points);
  position.reserve(points * 3);
  indices.reserve(points);

  for (auto &p : frame.data)
    _data.insert(_data.cend(), p.data.cbegin(), p.data.cend());

  for (const auto &patch : patches) {
    for (auto k = patch.K1; k <= patch.K2; ++k) {
      for (auto j = patch.J1; j <= patch.J2; ++j) {
        for (auto i = patch.I1; i <= patch.I2; ++i) {
          position.push_back(i);
          position.push_back(j);
          position.push_back(k);
          indices.push_back(indices.size());
        }
      }
    }
  }
  assert(_data.size() == points);
  assert(position.size() == points * 3);
  assert(indices.size() == points);

  return res;
}