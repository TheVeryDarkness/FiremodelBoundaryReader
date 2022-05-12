#pragma once
#include "fs.hpp"
#include "graphics.hpp"
#include "io.hpp"
#include "proxy.hpp"
#include "types.hpp"
#include <cassert>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <set>
#include <stack>
#include <vector>

using std::conditional_t;
using std::decay_t;
using std::endl;
using std::ifstream;
using std::ios;
using std::isinf;
using std::isnan;
using std::make_pair;
using std::map;
using std::ofstream;
using std::ostream;
using std::set;
using std::setprecision;
using std::stack;
using std::stringstream;
using std::tuple;
using std::filesystem::create_directory;
using std::filesystem::exists;

static inline u16 default_precision = 0;
template <bool enableDefaultPrecision = true>
static inline u16 input_precision() {
  if (enableDefaultPrecision && default_precision)
    return default_precision;
  u16 precision = 0;
  cout << "Precision:";
  cin >> precision;
  return precision;
}
static inline void set_default_precision(u16 precision) {
  default_precision = precision;
}

static inline void as_csv(ostream &o, const vector<float> &data, u32 m, u32 n,
                          const char *sep = ",") {
  u32 i = 0;
  for (auto v : data) {
    o << v;
    ++i;
    if (i % n == 0) {
      o << endl;
    } else
      o << sep;
  }
  o << endl;
}

template <size_t sz = 2>
static inline tuple<vector<float>, u32, u32>
from_csv(istream &in, const char (&sep)[sz] = ",") {
  tuple<vector<float>, u32, u32> res;
  auto &[vec, m, n] = res;

  while (in) {
    string line;
    getline(in, line);
    stringstream ss(line);

    size_t i = 0;
    float e;
    ss >> e;
    vec.push_back(e);
    while (ss) {
      check(in, sep);
      ss >> e;
      vec.push_back(e);
      ++i;
    }
    if (i != n)
      throw file_error();
  }

  return res;
}

static inline void print_patch(const vector<float> &data, u32 m, u32 n) {
  cout << endl;
  as_csv(cout.original(), data, m, n, " ");
}

static inline void save_patch_as_binary(const vector<float> &data,
                                        const vector<u32> &sizes,
                                        ofstream &fout) {
  write_vector_binary(fout, sizes);
  write_binary(fout, u32(0));
  write_vector_binary(fout, data);
}

struct function_info {
  using func = float (&)(float);
  const char *desc;
  float (&fun)(float);
  bool (*pred)(float) = nullptr;
};

bool positive(float num) noexcept { return num > 0; }
bool more_than_minus_one(float num) noexcept { return num > -1; }

static inline const map<string, function_info> &get_function_map() {
  using func = function_info::func;
  static map<string, function_info> map;
  if (map.empty()) {
    map.emplace("log", function_info{"Natural logarithm", std::log, positive});
    map.emplace("log1p", function_info{"Natural logarithm of plus 1",
                                       std::log1p, more_than_minus_one});
  }
  return map;
}

static bool disableDataCheck = false;

static inline void apply_function(vector<float> &data) {
  string func;
  cin >> func;
  auto &f_map = get_function_map();
  for (auto &[name, info] : f_map) {
    auto &[desc, _0, _1] = info;
    cout << setw(8) << name << " - " << desc << '\n';
  }
  auto iter = f_map.find(func);
  if (iter == f_map.end())
    return;
  auto &[name, tup] = *iter;
  auto &[desc, f, pred] = tup;
  if (pred && !disableDataCheck) {
    for (auto d : data)
      if (!pred(d)) {
        cout << "Function may not be able to be applied to these "
                "data.\nContinue? (y or n)";
        char opt = 'n';
        cin >> opt;
        if (opt != 'y')
          return;
        break;
      }
  }
  cout << desc << '\n';
  for (auto &d : data)
    d = f(d);
}

static inline void analysis_settings() {
  cout << "C - Disable data check." << disableDataCheck << endl;
  char opt = 'd';
  cin >> opt;
  switch (opt) {
  case 'C':
    cin >> disableDataCheck;
    break;
  default:
    COMMAND_NOT_FOUND;
    break;
  }
}

static inline tuple<u32, u32> get_matrix_size(const vector<u32> &sizes) {
  if (sizes.size() > 2) {
    cerr << "Data dimension " << sizes.size() << " > 2." << endl;
    return {0, 0};
  } else if (sizes.size() == 0) {
    return {0, 1};
  } else if (sizes.size() == 2) {
    auto m = sizes[0], n = sizes[1];
    return {m, n};
  } else {
    auto m = (u32)1, n = sizes[0];
    return {m, n};
  }
}

static inline void save_patch_as_csv_text(const vector<float> &data, u32 m,
                                          u32 n, ofstream &fout) {
  as_csv(fout, data, m, n);
}

static inline void reduce_dimension(vector<u32> &sizes, size_t target) {
  for (size_t i = 0; i < sizes.size() && sizes.size() > target; ++i) {
    if (sizes[i] == 1) {
      sizes.erase(sizes.cbegin() + i);
      --i;
    }
  }
}

/// @brief Get data offset in the flattened multi-dimension array.
static inline u32 array_offset(const vector<u32> &sizes,
                               const vector<u32> &index) {
  u32 offset = 0;
  auto s1 = sizes.begin(), s2 = sizes.end();
  auto i1 = index.begin(), i2 = index.end();
  for (; s1 != s2 && i1 != i2; ++s1, ++i1)
    offset *= *s1, offset += *i1;
  return offset;
}

/// @return Stride on the dimension
static inline u32 array_stride(const vector<u32> &sizes, size_t d) {
  u32 res = 1;
  auto p = sizes.begin() + d, end = sizes.end();
  for (; p != end; ++p)
    res *= *p;
  return res;
}

static inline vector<float> average(const vector<float> &data,
                                    vector<u32> &sizes, size_t dimension) {

  const auto m = sizes[dimension], n = array_stride(sizes, dimension + 1);

  vector<float> new_result;
  new_result.resize(data.size() / m);

  const auto &_data = data;
  const u32 l = (u32)data.size() / m / n;
  assert(data.size() == l * m * n);
#define SIMPLE_AVERAGE
#ifdef SIMPLE_AVERAGE
  // Will this be optimized well?
  for (size_t i = 0; i < l; ++i) {
    for (size_t k = 0; k < n; ++k) {
      for (size_t j = 0; j < m; ++j) {
        new_result[i * n + k] += data[i * m * n + j * n + k];
      }
    }
  }
#else
  for (const float *p0 = new_result.data(), // [i][]
           *const end0 = new_result.data() + new_result.size(),
                   *_p0 = data.data(); // [i][][]
       p0 != end0; p0 += n, _p0 += m * n) {
#ifndef NDEBUG
    assert(p0 >= new_result.data());
    const u32 i = (u32)(p0 - new_result.data()) / n;
    const u32 di = (u32)(p0 - new_result.data()) % n;
    assert(di == 0);
    assert(_p0 == data.data() + i * m * n);
#endif // !NDEBUG
    for (auto *p1 = p0, // [i][k]
             *const end1 = p0 + n,
              *_p1 = _p0; // [i][][k]
         p1 != end1; ++p1, ++_p1) {
#ifndef NDEBUG
      const u32 k = (u32)(p1 - p0);
      assert(0 <= k && k < n);
      assert(p1 == new_result.data() + i * n + k);
      assert(_p1 == data.data() + i * m * n + k);
#endif // !NDEBUG
      for (auto *_p2 = _p1, *const _end2 = _p1 + m * n; _p2 != _end2;
           _p2 += n) {
#ifndef NDEBUG
        const u32 j = (u32)(_p2 - _p1) / n;
        const u32 dj = (u32)(_p2 - _p1) % n;
        assert(0 <= j && j < m);
        assert(dj == 0);
        assert(p1 == new_result.data() + i * n + k);
        assert(_p2 == data.data() + i * m * n + j * n + k);

        assert(new_result.data() <= p1);
        assert(p1 < new_result.data() + new_result.size());
        assert(data.data() <= _p2);
        assert(_p2 < data.data() + data.size());
#endif // !NDEBUG
        *const_cast<float *>(p1) += *_p2; // [i][k] [i][j][k]
      }
    }
  }
#endif // SIMPLE_AVERAGE

  for (auto &e : new_result)
    e /= m;
  return new_result;
}

static inline vector<float> sample(const vector<float> &data,
                                   const vector<u32> &sizes, size_t dimension,
                                   const vector<u32> &pos) {

  const auto m = sizes[dimension], n = array_stride(sizes, dimension + 1),
             _m = (u32)pos.size();

  vector<float> new_result;
  new_result.resize(data.size() / m * _m);

  const auto &_data = data;
  const u32 l = (u32)data.size() / m / n;
  assert(data.size() == l * m * n);

  for (auto p : pos)
    assert(p < m);

  // Will this be optimized well?
  for (size_t i = 0; i < l; ++i) {
    for (size_t k = 0; k < n; ++k) {
      for (size_t _j = 0; _j < _m; ++_j) {
        new_result[i * _m * n + _j * n + k] = data[i * m * n + pos[_j] * n + k];
      }
    }
  }

  return new_result;
}

static inline void select_patch(const vector<patch_info> &patches) {
  u32 p = 0;
  cout << "Select a patch to analyze. Input an invalid index to discard."
       << endl;

  cout << "Patch index: ";
  cin >> p;
  if (p >= patches.size()) {
    return;
  }
  selected_patch = p;
}

static inline void interlop(vector<float> &a, const vector<float> &b, float wa,
                            float wb) {
  assert(a.size() == b.size());
  auto pa = a.begin(), ea = a.end();
  auto pb = b.begin(), eb = b.end();
  for (; pa != ea; ++pa, ++pb) {
    assert(pb != eb);
    *pa = (*pa * wa + *pb * wb) / (wa + wb);
  }
  assert(pb == eb);
}

static inline void analyze(const vector<patch_info> &patches,
                           const fds_boundary_file &frames) {
  vector<pair<vector<float>, vector<u32>>> data_and_size;
  // vector<float> data;
  // vector<u32> sizes;

  const auto frames_empty = frames.data.empty();

  while (true) {
    cout << R"(
Commands:
_ - Analysis settings.
m - Show current dimensions.
P - Set default float-point number precision.
w - Input in console to overlap current result.
+ - Push.)";

    if (!data_and_size.empty())
      cout <<
          R"(
y - Apply math functions to all cells in current result.
l - Flatten current result to specified dimension.
a - Calculate average on specified dimension of current result.
b - Calculate border of current result.
p - Print current result to console.
B - Save current result as binary to file.
C - Save current result as CSV file.
n - Replace nan with specified value.
N - Replace inf with specified value.
S - Sample.
c - Read CSV file as current result.
t - Get frame at a time. May interlop between two frame.
- - Pop.)";

    if (data_and_size.size() >= 2) {
      auto rp = data_and_size.rbegin();
      if ((rp + 0)->second == (rp + 1)->second)
        cout << R"(
L - Interlop.)";
    }

    if (!patches.empty())
      cout << R"(
s - Select a patch.
M - Visualize merged patches.)";

    if (!frames_empty && selected_patch < patches.size())
      cout <<
          R"(
f - Copy patch data in a frame as current result.
F - Copy data in all frames of this patch as current result.)";

#if GRAPHICS_ENABLED
    if (!data_and_size.empty())
      cout << R"(
v - Visualize current result.)";
#endif // GRAPHICS_ENABLED

    cout <<
        R"(
d - Discard.
)";
    char opt = 'd';
    cin >> opt;
    switch (opt) {
    case 'd':
      return;
    case '_':
      analysis_settings();
      break;
    case '+':
      data_and_size.push_back({});
      break;
    case '-':
      data_and_size.pop_back();
      break;
    case 'n':
      if (!data_and_size.empty()) {
        float val = 0.f;
        cin >> val;
        size_t i = 0;
        for (auto &cell : data_and_size.back().first)
          if (isnan(cell)) {
            cell = val;
            ++i;
          }
        if (i > 0)
          clog << "Replaced " << i << " nan with " << val << '\n';
      }
      break;
    case 'N':
      if (!data_and_size.empty()) {
        float val = 0.f;
        cin >> val;
        size_t i = 0;
        for (auto &cell : data_and_size.back().first)
          if (isinf(cell)) {
            cell = val;
            ++i;
          }
        if (i > 0)
          clog << "Replaced " << i << " inf with " << val << '\n';
      }
      break;
    case 'z':
      if (!data_and_size.empty()) {
        size_t i = 0;
        for (auto &cell : data_and_size.back().first)
          if (fpclassify(cell) == FP_SUBNORMAL) {
            ++i;
          }
        clog << "Found " << i
             << " subnormal float-point numbers. This typically means some "
                "errors.\n";
      }
      break;
    case 's':
      select_patch(patches);
      break;
    case 'w':
      if (!data_and_size.empty()) {
        cout << "Size(0 to stop): ";
        u32 sz = 0;

        auto &[data, sizes] = data_and_size.back();
        sizes.clear();
        do {
          sz = 0;
          cin >> sz;
          if (sz == 0)
            break;
          sizes.push_back(sz);
        } while (true);

        u32 data_size = array_stride(sizes, 0);

        float e = 0.f;
        data.clear();
        data.reserve(data_size);
        while (data.size() < data_size) {
          cin >> e;
          data.push_back(e);
        }
      }
      break;
    case 'm':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        if (sizes.empty())
          break;
        auto &out = cout.original();
        auto begin = sizes.cbegin(), end = sizes.cend();
        out << *begin;
        ++begin;
        for (; begin != end; ++begin)
          out << '*' << *begin;
      }
      break;
    case 'p':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        u16 precision = input_precision();
        cout << setprecision(precision);

        auto [m, n] = get_matrix_size(sizes);
        if (n != 0)
          print_patch(data, m, n);
      }
      break;
    case 'P': {
      u16 precision = input_precision<false>();
      set_default_precision(precision);
    } break;
    case 'S':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        size_t dimension = -1;
        cin >> dimension;
        if (dimension >= sizes.size()) {
          cerr << "Dimension invalid." << endl;
          break;
        }
        const auto sz = sizes[dimension];
        cout << "Size of this dimension: " << sz << endl;
        cout << "'\\' to stop\n";
        vector<u32> pos = read_until<u32>(cin);
        for (auto p : pos)
          if (p >= sz) {
            cerr << "Out of border.\n";
            break;
          }

        data = sample(data, sizes, dimension, pos);
        sizes[dimension] = (u32)pos.size();
      }
      break;
#if GRAPHICS_ENABLED
    case 'v':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        auto [m, n] = get_matrix_size(sizes);
        if (!n)
          break;
        visualize_data(data, m, n);
      }
      break;
#endif // GRAPHICS_ENABLED
    case 'B':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        cout << "File name: ";
        string path;
        getline(cin, path);
        ofstream fout(path, ios::binary);
        if (!fout) {
          FILE_OPEN_FAILED;
          break;
        }
        save_patch_as_binary(data, sizes, fout);
      }
      break;
    case 'c':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        auto opt = request_file_by_name(
            [](const path &p) { return is_directory(p); }, "csv file");
        if (opt.has_value()) {
          auto &path = opt.value();
          ifstream fin = ifstream(path);
          if (!fin) {
            FILE_OPEN_FAILED;
            break;
          }

          auto [_data, _m, _n] = from_csv(fin);
          data = std::move(_data);
          sizes = {_m, _n};
        }
      }
      break;
    case 'C':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        optional<path> opt =
            request_file_by_name([](const path &p) { return true; }, "output");
        if (!opt)
          break;
        auto &path = opt.value();
        ofstream fout(path);
        if (!fout) {
          FILE_OPEN_FAILED;
          break;
        }

        u16 precision = input_precision();
        fout << setprecision(precision);

        auto [m, n] = get_matrix_size(sizes);
        save_patch_as_csv_text(data, m, n, fout);
      }
      break;
    case 'y':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        apply_function(data);
      }
      break;
    case 'l':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        size_t d = 2;
        cin >> d;
        reduce_dimension(sizes, d);
        if (sizes.size() > d)
          cerr << "Flattening failed. Current dimension is " << sizes.size()
               << ".\n";
      }
      break;
    case 'a':
      if (!data_and_size.empty()) {
        auto &[data, sizes] = data_and_size.back();
        // Consider [l][m][n]
        size_t dimension = -1;
        cin >> dimension;
        if (dimension >= sizes.size() || sizes.size() == 1) {
          cerr << "Dimension invalid.\n";
          break;
        }
        data = average(data, sizes, dimension);
        sizes.erase(sizes.cbegin() + dimension);
      }
      break;
    case 'b':
      if (!data_and_size.empty()) {
        const auto &[data, size] = data_and_size.back();
        auto [pmin, pmax] = minmax_element(data.begin(), data.end());
        clog << "Range: " << *pmin << " - " << *pmax << '\n';
      }
      break;
    case 'L':
      if (data_and_size.size() >= 2) {
        auto rp = data_and_size.rbegin();
        if ((rp + 0)->second == (rp + 1)->second) {
          auto latter = std::move(data_and_size.back());
          float a = 1.f;
          float b = 0.f;
          cin >> a >> b;
          if (b != 0.0f)
            interlop(data_and_size.back().first, latter.first, a, b);
        }
      }
      break;
    case 't':
      if (!frames_empty && selected_patch < patches.size()) {
        float t = 0;
        cin >> t;

        auto &_data = frames.data.at(selected_patch);
        auto &data = _data.data;
        auto sz = _data.size;
        if (t < frames.times.front() || t > frames.times.back())
          cerr << "Out of border. Terminated.\n";
        else
          for (size_t i = 0; i < frames.times.size(); ++i)
            if (t < frames.times[i]) {
              assert(i > 0);
              vector<float> last = {data.begin() + (i - 1) * sz,
                                    data.begin() + i * sz};
              const vector<float> next = {data.begin() + i * sz,
                                          data.begin() + (i + 1) * sz};
              interlop(last, next, frames.times[i] - t,
                       t - frames.times[i - 1]);
              clog << "Time: " << frames.times[i - 1] << ", " << frames.times[i]
                   << "\n";
              auto &patch = patches[selected_patch];
              data_and_size.push_back(
                  make_pair(std::move(last),
                            vector<u32>{patch.K(), patch.J(), patch.I()}));
              break;
            } else if (t == frames.times[i]) {
              auto &patch = patches[selected_patch];
              clog << "Time: " << frames.times[i] << "\n";
              data_and_size.push_back(
                  make_pair(vector<float>{data.begin() + i * sz,
                                          data.begin() + (i + 1) * sz},
                            vector<u32>{patch.K(), patch.J(), patch.I()}));
              break;
            }
      }
      break;
    case 'f':
      if (!frames_empty && selected_patch < patches.size()) {
        auto item = make_pair(vector<float>{}, vector<u32>{});
        auto &[data, sizes] = item;
        const auto &patch = patches[selected_patch];
        cout << "Frame index: ";
        auto f = frames.times.size();
        cin >> f;
        if (f >= frames.times.size()) {
          cerr << "Frame index invalid." << endl;
          break;
        }
        data.clear();
        auto &_patch = frames.data.at(selected_patch);
        auto &patch_data = _patch.data;
        auto &patch_size = _patch.size;
        data.insert(data.end(), patch_data.begin() + patch_size * f,
                    patch_data.begin() + patch_size * (f + 1));
        sizes = {patch.K(), patch.J(), patch.I()};
        data_and_size.push_back(std::move(item));
      }
      break;
    case 'F':
      if (!frames_empty && selected_patch < patches.size()) {
        auto item = make_pair(vector<float>{}, vector<u32>{});
        auto &[data, sizes] = item;
        const auto &patch = patches[selected_patch];
        sizes = {(u32)frames.times.size(), patch.K(), patch.J(), patch.I()};
        data.reserve(array_stride(sizes, 0));
        data = frames.data.at(selected_patch).data;
        data_and_size.push_back(std::move(item));
      }
      break;
    default:
      COMMAND_NOT_FOUND;
      break;
    }
  }
}