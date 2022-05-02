#pragma once
#include "types.hpp"
#include <algorithm>
#include <array>
#include <bitset>
#include <cassert>
#include <iostream>
#include <map>
#include <set>
#include <tuple>

using std::abs;
using std::all_of;
using std::any_of;
using std::array;
using std::bitset;
using std::clog;
using std::conditional_t;
using std::invoke_result_t;
using std::lexicographical_compare;
using std::lroundf;
using std::map;
using std::max;
using std::max_element;
using std::min;
using std::min_element;
using std::minmax_element;
using std::numeric_limits;
using std::pair;
using std::set;
using std::sort;
using std::tuple;

template <typename Ty> static inline bool all_same(const vector<Ty> &vec) {
  return all_of(vec.begin(), vec.end(),
                [&vec](auto v) { return v == vec.front(); });
}

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

  template <size_t sz> u32 length() const noexcept {
    switch (sz) {
    case 0:
      return I();
    case 1:
      return J();
    case 2:
      return K();
    default:
      static_assert(sz == 0 || sz == 1 || sz == 2);
      return numeric_limits<u32>::max();
    }
  }

  template <size_t sz> pair<u32, u32> border() const noexcept {
    switch (sz) {
    case 0:
      return {I1, I2};
    case 1:
      return {J1, J2};
    case 2:
      return {K1, K2};
    default:
      static_assert(sz == 0 || sz == 1 || sz == 2);
      return {u32(0), numeric_limits<u32>::max()};
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

static inline size_t selected_patch = 0;

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

static inline vector<float> from_frame(const frame &frame, u32 size) {
  vector<float> res;
  res.reserve(size);

  for (auto &p : frame.data)
    res.insert(res.cend(), p.data.cbegin(), p.data.cend());
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

  position.reserve(points * 3);
  indices.reserve(points);

  _data = from_frame(frame, points);

  for (const auto &patch : patches) {
    for (auto k = patch.K1; k <= patch.K2; ++k) {
      for (auto j = patch.J1; j <= patch.J2; ++j) {
        for (auto i = patch.I1; i <= patch.I2; ++i) {
          position.push_back(i);
          position.push_back(j);
          position.push_back(k);
          indices.push_back((u32)indices.size());
        }
      }
    }
  }
  assert(_data.size() == points);
  assert(position.size() == points * 3);
  assert(indices.size() == points);

  return res;
}

struct patch_region {
  u32 P1;
  u32 P2;
  u32 Q1;
  u32 Q2;

public:
  constexpr patch_region(u32 P1, u32 P2, u32 Q1, u32 Q2) noexcept
      : P1(P1), P2(P2), Q1(Q1), Q2(Q2) {}
  explicit constexpr patch_region() noexcept : patch_region(0, 0, 0, 0) {}
  constexpr patch_region(const patch_region &) noexcept = default;
  constexpr static patch_region from_patch(const patch_info &patch) {
    switch (patch.IOR) {
    case 1:
    case -1:
      return {patch.J1, patch.J2, patch.K1, patch.K2};
    case 2:
    case -2:
      return {patch.I1, patch.I2, patch.K1, patch.K2};
    case 3:
    case -3:
      return {patch.I1, patch.I2, patch.J1, patch.J2};
    default:
      assert(false);
      return patch_region{};
    }
  }

  template <u32 i> constexpr tuple<u32, u32> border() const {
    switch (i) {
    case 0:
      return {P1, P2};
    case 1:
      return {Q1, Q2};
    default:
      static_assert(i == 0 || i == 1);
      assert(false);
      return {0, numeric_limits<u32>::max()};
    }
  }
  static bool in_range(u32 x, u32 x1, u32 x2) { return x1 <= x && x <= x2; }
  static bool xor_range(u32 x, u32 x1, u32 x2, u32 X1, u32 X2) {
    return in_range(x, x1, x2) ^ in_range(x, X1, X2);
  }
  bitset<9> cross(const array<u32, 4> &P, const array<u32, 4> &Q) const {
    bitset<9> res;
    for (size_t i = 0; i < 3; ++i)
      if (P1 * 2 <= P[i] + P[i + 1] && P[i] + P[i + 1] <= P2 * 2)
        for (size_t j = 0; j < 3; ++j)
          if (Q1 * 2 <= Q[j] + Q[j + 1] && Q[j] + Q[j + 1] <= Q2 * 2)
            res[i * 3 + j] = true;

    return res;
  }
  [[nodiscard]] array<patch_region, 4> cross(const patch_region &p) const {
    array<patch_region, 4> res = {patch_region{}, patch_region{},
                                  patch_region{}, patch_region{}};
    array<u32, 4> x = {P1, P2, p.P1, p.P2};
    sort(x.begin(), x.end());
    array<u32, 4> y = {Q1, Q2, p.Q1, p.Q2};
    sort(y.begin(), y.end());
    bitset<9> a = cross(x, y);
    bitset<9> b = p.cross(x, y);
    bitset<9> c = a | b;

    size_t count = 0;
    const auto _push = [&count, &res](u32 x1, u32 x2, u32 y1, u32 y2) {
      assert(count < 4);
      res[count] = {x1, x2, y1, y2};
      assert(res[count].legal());
      ++count;
    };
    for (size_t i = 0; i < 3; ++i)
      if (x[i] != x[i + 1]) {
        if (c[3 * i] && c[3 * i + 1] && c[3 * i + 2])
          _push(x[i], x[i + 1], y[0], y[3]);
        else if (c[3 * i] && c[3 * i + 1])
          _push(x[i], x[i + 1], y[0], y[2]);
        else if (c[3 * i + 1] && c[3 * i + 2])
          _push(x[i], x[i + 1], y[1], y[3]);
        else
          for (size_t j = 0; j < 3; ++j)
            if (c[3 * i + j]) {
              _push(x[i], x[i + 1], y[j], y[j + 1]);
            }
      }
    return res;
  }
  /// @brief Merge another region with this. Changes are made only if true is
  /// returned.
  /// @param p
  /// @return Whether successfully merged.
  [[nodiscard]] constexpr bool merge(const patch_region &p) {
    auto [p1, p2] = p.border<0>();
    auto [q1, q2] = p.border<1>();
    if (P1 == p1 && P2 == p2) {
      if (Q1 == q2) {
        Q1 = q1;
        return true;
      }
      if (Q2 == q1) {
        Q2 = q2;
        return true;
      }
    }
    if (Q1 == q1 && Q2 == q2) {
      if (P1 == p2) {
        P1 = p1;
        return true;
      }
      if (P2 == p1) {
        P2 = p2;
        return true;
      }
    }
    return false;
  }

  /// @brief Check if one region is near another
  constexpr bool near(const patch_region &p) const {
    auto [x1, x2] = border<0>();
    auto [y1, y2] = border<1>();
    auto [X1, X2] = p.border<0>();
    auto [Y1, Y2] = p.border<1>();
    if (X1 < x2 && x1 < X2)
      if (y1 == Y2 || y2 == Y1)
        return true;
    if (Y1 < y2 && y1 <= Y2)
      if (x1 == X2 || x2 == X1)
        return true;
    return false;
  }

  constexpr bool contains(u32 p, u32 q) const noexcept {
    return P1 <= p && p <= P2 && Q1 <= q && q <= Q2;
  }
  /// @brief Check if one region intersects with another
  /// @param that
  /// @return
  constexpr bool intersect(const patch_region &that) const noexcept {
    return P1 < that.P2 && that.P1 < P2 && Q1 < that.Q2 && that.Q1 < Q2;
  }
  constexpr bool empty() const noexcept {
    auto res = P1 == P2 || Q1 == Q2;
    assert(!res || (P1 == 0 && P2 == 0));
    return res;
  }
  constexpr bool legal() const noexcept { return !empty(); }

  constexpr tuple<u32, u32, u32, u32> minmax() const noexcept {
    return {P1, P2, Q1, Q2};
  }
  /// @brief
  /// @param p
  /// @return Near; Changes made; fragments
  [[nodiscard]] constexpr tuple<bool, bool, array<patch_region, 4>>
  try_push(const patch_region &p) {
    constexpr array<patch_region, 4> null = {patch_region{}, patch_region{},
                                             patch_region{}, patch_region{}};

    if (intersect(p)) {
      auto &&c = cross(p);

      if (!c[0].empty()) {
        for (auto rp = c.rbegin(), re = c.rend(); rp != re; ++rp) {
          if (!rp->empty()) {
            *this = *rp;
            *rp = patch_region{};
            break;
          }
        }
        return {false, true, c};
      }
      assert(c[0].empty());
      assert(c[1].empty());
      assert(c[2].empty());
      assert(c[3].empty());
    }
    if (merge(p))
      return {false, true, null};
    if (near(p)) {
      return {true, false, null};
    }
    return {false, false, null};
  }
};

constexpr static array<patch_region, 4> null_region = {
    patch_region{}, patch_region{}, patch_region{}, patch_region{}};

struct connected_region {
  vector<patch_region> regions;

  /// @brief
  /// @param p Patch region to merge
  /// @return true if changes are made, region fragments are returned if only
  /// changes are made
  [[nodiscard]] tuple<bool, array<patch_region, 4>>
  push(const patch_region &p) {
    for (auto &region : regions) {
      auto [near, suc, frags] = region.try_push(p);
      if (suc)
        return {suc, frags};
      if (!suc)
        assert(all_of(frags.begin(), frags.end(),
                      [](const patch_region &pr) { return pr.empty(); }));
      if (near) {
        regions.push_back(p);
        return {true, null_region};
      }
    }
    return {false, null_region};
  }

  bool contains(u32 P, u32 Q) const noexcept {
    for (auto &reg : regions) {
      if (reg.contains(P, Q))
        return true;
    }
    return false;
  }
  bool contains_all(const vector<u32> &P, const vector<u32> &Q) const noexcept {
    assert(P.size() == Q.size());
    auto p1 = P.begin(), e1 = P.end();
    auto p2 = Q.begin(), e2 = Q.end();

    for (; p1 != e1; ++p1, ++p2) {
      assert(p2 != e2);
      for (auto &reg : regions) {
        if (!reg.contains(*p1, *p2))
          return false;
      }
    }
    assert(p1 == e1 && p2 == e2);
    return true;
  }

  bool legal() const {
    return !empty() &&
           all_of(regions.begin(), regions.end(),
                  [](const patch_region &pr) { return !pr.empty(); });
  }
  bool empty() const { return regions.empty(); }
  tuple<u32, u32, u32, u32> minmax() const {
    auto [P1, P2, Q1, Q2] = regions.front().minmax();
    for (auto &cr : regions) {
      const auto [p1, p2, q1, q2] = regions.front().minmax();
      P1 = min(P1, p1);
      P2 = min(P2, p2);
      Q1 = min(Q1, q1);
      Q2 = min(Q2, q2);
    }
    return {P1, P2, Q1, Q2};
  }

  static bool insert(vector<patch_region> &fragments,
                     const array<patch_region, 4> &arr) {
#ifndef NDEBUG
    bool empty_found = false;
    for (auto &f : arr)
      if (!f.empty())
        assert(!empty_found);
      else
        empty_found = true;
#endif // !NDEBUG

    bool not_empty = false;
    for (auto &f : arr)
      if (!f.empty()) {
        fragments.push_back(f);
        not_empty = true;
      } else
        break;
    return not_empty;
  }
  void merge_fragments(vector<patch_region> &&fragments) {
  RESTART:
    for (auto p = fragments.cbegin(); p != fragments.cend(); ++p) {
      auto &frag = *p;
      auto &&[suc, frags] = push(frag);
      if (suc) {
        insert(fragments, frags);
        fragments.erase(p);
        goto RESTART;
      }
    }
    assert(fragments.empty());
  }
  void merge() {
    vector<patch_region> fragments;
  RESTART:
    for (auto p1 = regions.begin(), e = regions.end(); p1 != e; ++p1) {
      for (auto p2 = p1 + 1; p2 != e; ++p2) {
        auto [near, suc, frags] = p1->try_push(*p2);
        insert(fragments, frags);
        if (suc) {
          regions.erase(p2);
          goto RESTART;
        }
      }
    }
    merge_fragments(std::move(fragments));
  }
  [[nodiscard]] bool merge(const connected_region &that) {
    assert(that.legal());

    bool changed = false;
    vector<patch_region> fragments;

    for (const auto &r : that.regions) {
      auto &&[suc, frag] = push(r);
      if (suc)
        changed = true;
      else {
        assert(all_of(frag.cbegin(), frag.cend(),
                      [](const patch_region &pr) { return pr.empty(); }));
      }
      insert(fragments, frag);
    }
    merge_fragments(std::move(fragments));

    return changed;
  }
};

struct connected_regions {
  vector<connected_region> regions;

  void push(const patch_region &p) {
    for (auto &region : regions) {
      auto &&[suc, frag] = region.push(p);

      if (!suc) {
        assert(frag[0].empty());
        assert(frag[1].empty());
        assert(frag[2].empty());
        assert(frag[3].empty());
      } else {
        for (auto &f : frag)
          if (!f.empty())
            push(f);
        return;
      }
    }
    regions.push_back(connected_region{{p}});
  }
  void merge() {
  RESTART:
    for (auto p1 = regions.begin(), e = regions.end(); p1 != e; ++p1) {
      for (auto p2 = p1 + 1; p2 != e; ++p2) {
        if (p1->merge(*p2)) {
          regions.erase(p2);
          goto RESTART;
        }
      }
    }
    for (auto &cr : regions)
      cr.merge();
  }
  bool in_single_connected_region(u32 P, u32 Q) {
    return any_of(
        regions.cbegin(), regions.cend(),
        [&P, &Q](const connected_region &cr) { return cr.contains(P, Q); });
  }
  bool in_single_connected_region(const vector<u32> &P, const vector<u32> &Q) {
    return any_of(
        regions.cbegin(), regions.cend(),
        [&P, &Q](const connected_region &cr) { return cr.contains_all(P, Q); });
  }
  bool legal() const {
    return !empty() &&
           all_of(regions.begin(), regions.end(),
                  [](const connected_region &pr) { return pr.legal(); });
  }
  bool empty() const { return regions.empty(); }
  tuple<u32, u32, u32, u32> minmax() const {
    auto [P1, P2, Q1, Q2] = regions.front().minmax();
    for (auto &cr : regions) {
      const auto [p1, p2, q1, q2] = regions.front().minmax();
      P1 = min(P1, p1);
      P2 = min(P2, p2);
      Q1 = min(Q1, q1);
      Q2 = min(Q2, q2);
    }
    return {P1, P2, Q1, Q2};
  }
};

static const float tolerance = .01f;

static inline bool out_of_tolerance(float n, long i) {
  return abs(n - i) > tolerance;
}

class regions {
  // X Y Z
  array<map<u32, connected_regions>, 3> r;

  template <size_t dim> void push(const patch_info &patch) {
    static_assert(dim < 3);
    auto &m = r[dim];
    auto [R, _] = patch.border<dim>();
    assert(R == _);
    auto p = m.find(R);
    auto region = patch_region::from_patch(patch);
    if (p == m.cend()) {
      m.emplace(R, connected_regions{{connected_region{{region}}}});
    } else {
      connected_regions &cr = p->second;
      assert(cr.legal());
#ifndef NDEBUG
      auto copy = cr;
      copy.push(region);
      if (!copy.legal()) {
        copy = cr;
        copy.push(region);
      }
#endif // !NDEBUG

      cr.push(region);
      assert(cr.legal());
    }
  }

  template <size_t dim>
  vector<set<u32>>
  filter_for_each_connected_region(const vector<float> &P,
                                   const vector<float> &Q,
                                   const vector<float> &R) const {
    vector<set<u32>> res;
    u32 i_node = 0;
    auto &m = r[dim];
    for (const auto &[_r, crs] : m) {
      for (auto &cr : crs.regions) {
        set<u32> current;
        auto p1 = P.begin();
        const auto p2 = P.end();
        auto q1 = Q.begin();
        const auto q2 = Q.end();
        auto r1 = R.begin();
        const auto r2 = R.end();
        u32 i_node = 0;
        for (; p1 != p2; ++p1, ++q1, ++r1, ++i_node) {
          assert(q1 != q2);
          assert(r1 != r2);
          long r = lroundf(*r1);
          if (out_of_tolerance(*r1, r) || r != _r)
            continue;
          long p = lroundf(*p1);
          long q = lroundf(*q1);
          if (cr.contains(p, q))
            current.insert(i_node);
        }
        assert(p1 == p2);
        assert(q1 == q2);
        assert(r1 == r2);
        res.push_back(std::move(current));
      }
    }
    return res;
  }

public:
  regions(regions &&) = default;
  regions() = default;

  void push(const patch_info &p) {
    switch (p.IOR) {
    case 1:
    case -1: {
      push<0>(p);
    } break;
    case 2:
    case -2: {
      push<1>(p);
    } break;
    case 3:
    case -3: {
      push<2>(p);
    } break;
    default:
      assert(false);
      break;
    }
  }
  void merge() {
    for (auto &m : r)
      for (auto &[R, crs] : m)
        crs.merge();
  }
  template <size_t dim>
  bool in_single_connected_region(u32 R, vector<u32> P, vector<u32> Q) {
    auto &m = r[dim];
    auto iter = m.find(R);
    if (iter == m.cend())
      return false;
    return iter->second.in_single_connected_region(P, Q);
  }
  bool in_single_connected_region(vector<u32> I, vector<u32> J, vector<u32> K) {
    assert(I.size() == J.size());
    assert(I.size() == K.size());
    bool i = all_same(I);
    bool j = all_same(J);
    bool k = all_same(K);
    assert(!i || !j || !k);
    assert((i ^ j ^ k) || (!i && !j && !k));
    if (i)
      return in_single_connected_region<0>(I.front(), J, K);
    if (j)
      return in_single_connected_region<1>(J.front(), I, K);
    if (k)
      return in_single_connected_region<2>(K.front(), I, J);
    return false;
  }

  vector<set<u32>>
  filter_for_each_connected_region(const vector<float> &nodes) {
    vector<float> x, y, z;
    x.reserve(nodes.size() / 3);
    y.reserve(nodes.size() / 3);
    z.reserve(nodes.size() / 3);
    for (auto iter = nodes.begin(), end = nodes.end(); iter != end;) {
      x.push_back(((*iter++ - mesh.x0) / mesh.cell_size));
      y.push_back(((*iter++ - mesh.y0) / mesh.cell_size));
      z.push_back(((*iter++ - mesh.z0) / mesh.cell_size));
    }

    vector<set<u32>> res = filter_for_each_connected_region<0>(y, z, x);
    {
      auto _0 = filter_for_each_connected_region<1>(x, z, y);
      res.insert(res.cend(), _0.cbegin(), _0.cend());
    }
    {
      auto _1 = filter_for_each_connected_region<2>(x, y, z);
      res.insert(res.cend(), _1.cbegin(), _1.cend());
    }
    return res;
  }
  auto begin() const { return r.begin(); }
  auto end() const { return r.end(); }
  template <size_t dim> tuple<u32, u32, u32, u32, u32, u32> minmax() const {
    auto &m = r[dim];
    auto &first = *m.begin();
    auto R1 = first.first;
    auto R2 = first.first;
    auto [P1, P2, Q1, Q2] = first.second.minmax();
    for (auto &[R, crs] : m) {
      auto [p1, p2, q1, q2] = crs.minmax();
      P1 = min(P1, p1);
      P2 = max(P2, p2);
      Q1 = min(Q1, q1);
      Q2 = max(Q2, q2);
      R1 = min(R1, R);
      R2 = max(R2, R);
    }
    return {P1, P2, Q1, Q2, R1, R2};
  }
  tuple<u32, u32, u32, u32, u32, u32> minmax() const {
    auto [y01, y02, z01, z02, x01, x02] = minmax<0>();
    auto [x11, x12, z11, z12, y11, y12] = minmax<1>();
    auto [x21, x22, y21, y22, z21, z22] = minmax<2>();
    return {min({x01, x11, x21}), max({x02, x12, x22}), min({y01, y11, y21}),
            max({y02, y12, y22}), min({z01, z11, z21}), max({z02, z12, z22})};
  }
};

template <size_t i>
static inline void push_patch_region(vector<u32> &vec, const patch_region &reg,
                                     u32 R) {
  if constexpr (i == 0)
    vec.push_back(R);
  vec.push_back(reg.P1);
  if constexpr (i == 1)
    vec.push_back(R);
  vec.push_back(reg.Q1);
  if constexpr (i == 2)
    vec.push_back(R);

  if constexpr (i == 0)
    vec.push_back(R);
  vec.push_back(reg.P1);
  if constexpr (i == 1)
    vec.push_back(R);
  vec.push_back(reg.Q2);
  if constexpr (i == 2)
    vec.push_back(R);

  if constexpr (i == 0)
    vec.push_back(R);
  vec.push_back(reg.P2);
  if constexpr (i == 1)
    vec.push_back(R);
  vec.push_back(reg.Q1);
  if constexpr (i == 2)
    vec.push_back(R);

  if constexpr (i == 0)
    vec.push_back(R);
  vec.push_back(reg.P2);
  if constexpr (i == 1)
    vec.push_back(R);
  vec.push_back(reg.Q2);
  if constexpr (i == 2)
    vec.push_back(R);
}

static inline tuple<tuple<vector<u32>>, vector<u32>>
from_patches(const regions &rgs, bool wireframe) {
  tuple<tuple<vector<u32>>, vector<u32>> res;
  auto &[_, indices] = res;
  auto &[position] = _;

  indices.reserve(wireframe ? indices.size() * 8 : indices.size() * 6);

  u32 IOR = 0;
  for (const auto &map : rgs) {
    for (const auto &[R, crs] : map) {
      for (const auto &cr : crs.regions) {
        for (const auto &reg : cr.regions) {
          const u32 N = (u32)position.size() / 3;
          if (wireframe) {
            indices.push_back(N + 0);
            indices.push_back(N + 1);
            indices.push_back(N + 1);
            indices.push_back(N + 3);
            indices.push_back(N + 3);
            indices.push_back(N + 2);
            indices.push_back(N + 2);
            indices.push_back(N + 0);
          } else {
            indices.push_back(N + 0);
            indices.push_back(N + 1);
            indices.push_back(N + 3);
            indices.push_back(N + 0);
            indices.push_back(N + 3);
            indices.push_back(N + 2);
          }
          switch (IOR) {
          case 0:
            push_patch_region<0>(position, reg, R);
            break;
          case 1:
            push_patch_region<1>(position, reg, R);
            break;
          case 2:
            push_patch_region<2>(position, reg, R);
            break;
          default:
            assert(false);
            break;
          }
        }
      }
    }

    ++IOR;
  }

  return res;
}

static inline tuple<float, float> frame_minmax(const frame &frame) {
  float _min = frame.data.front().data.front();
  float _max = frame.data.front().data.front();
  for (auto &patch_frame : frame.data) {
    auto &data = patch_frame.data;
    auto [__min, __max] = minmax_element(data.begin(), data.end());
    _min = min(_min, *__min);
    _max = max(_max, *__max);
  }
  return {_min, _max};
}

static inline tuple<float, float> data_minmax(const vector<float> &data) {
  auto [_min, _max] = minmax_element(data.cbegin(), data.cend());
  clog << '[' << *_min << ',' << *_max << "]\n";

  return {*_min, *_max};
}

static inline regions merge(const vector<patch_info> &patches) {
  regions res;
  for (auto &patch : patches) {
    res.push(patch);
  }
  res.merge();
  return res;
}