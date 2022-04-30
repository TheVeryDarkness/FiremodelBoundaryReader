/// @file main.cpp
/// @author Boyi Huang
/// @mainpage Written by Boyi Huang, who was a student in Tongji University when
/// he wrote this. C++17 is required for filesystem operation. GLAD and GLFW3 is
/// required for boundary patch display. For reading .bf files generated by fds.
/// See Source/dump.f90 in fds reposity for more information, especailly codes
/// that are writing data to LU_BNDF(NF,NM). See also section Boundary Files in
/// FDS_User_Guide.pdf. As for opengl, most of the codes about opengl are
/// adapted from https://github.com/JoeyDeVries/LearnOpenGL.
#include "analysis.hpp"
#include "apdl_io.hpp"
#include "elements_io.hpp"
#include "fds_io.hpp"
#include "fs.hpp"
#include "graphics.hpp"
#include "io.hpp"
#include "types.hpp"
#include <array>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>

using std::char_traits;
using std::conditional_t;
using std::cout;
using std::endl;
using std::getline;
using std::hex;
using std::ifstream;
using std::ios;
using std::optional;
using std::setprecision;
using std::string;
using std::filesystem::current_path;
using std::filesystem::directory_entry;
using std::filesystem::directory_iterator;
using std::filesystem::is_regular_file;
using std::filesystem::path;

static inline void search_frame_by_time(const vector<frame> &frames) {
  using std::cin;
  using std::cout;
  constexpr auto help = R"(
Choose search mode.
> - Earliest frame later than given time.
< - Latest frame earlier than given time.
+ - Earliest frame not earlier than given time.
- - Latest frame not later than given time.
d - Discard.
)";
  while (true) {
    cout << help;

    char op4 = 'd';
    cin >> op4;
    const auto time_search =
        [&frames = frames](bool (*pred)(float, float, float), bool right) {
          float t = NAN;
          cin >> t;
          for (size_t i = 1; i < frames.size(); ++i)
            if (pred(frames[i - 1].time, t, frames[i].time)) {
              cout << "Frame " << i << " at " << frames[right ? i : i - 1].time
                   << "s." << endl;
              return;
            }
          cout << "Not found." << endl;
        };
    switch (op4) {
    case '>': {
      time_search(
          [](float prev, float t, float cur) { return prev <= t && t < cur; },
          true);
    } break;
    case '+': {
      time_search(
          [](float prev, float t, float cur) { return prev < t && t <= cur; },
          true);
    } break;
    case '<': {
      time_search(
          [](float cur, float t, float next) { return cur <= t && t < next; },
          false);
    } break;
    case '-': {
      time_search(
          [](float cur, float t, float next) { return cur < t && t <= next; },
          false);
    } break;
    case 'd':
      return;
    default:
      cout << "Option not found." << endl;
      break;
    }
  }
}

static inline void search_frame(const vector<frame> &frames) {
  while (true) {
    cout << R"(
Select item to search by.
i - Index.
t - Time.
d - discard.
)";
    char op3 = 'd';
    cin >> op3;
    switch (op3) {
    case 'i': {
      cout << "Input the index. There are " << frames.size()
           << " frames in total." << endl;
      u32 f = 0;
      cin >> f;
      if (f >= frames.size()) {
        cout << "Not a valid frame." << endl;
        break;
      }
      cout << "Frame " << f << " is at " << frames[f].time << "s." << endl;
    } break;
    case 't': {
      search_frame_by_time(frames);
    } break;
    case 'd':
      return;
    default:
      cout << "Option not found." << endl;
      break;
    }
  }
}

static inline void search_patch(const vector<patch_info> &patches) {
  while (true) {
    cout << R"(
Select item to search by.
i - Index.
d - discard.
)";
    char op3 = 'd';
    cin >> op3;
    switch (op3) {
    case 'i': {
      cout << "Input the index. There are " << patches.size()
           << " patches in total." << endl;
      u32 p = 0;
      cin >> p;
      if (p >= patches.size()) {
        cout << "Not a valid patch." << endl;
        break;
      }
      cout << "Patch " << p << ": " << patches[p] << endl;
    } break;
    case 'd':
      return;
    default:
      cout << "Option not found." << endl;
      break;
    }
  }
}

/*
 * @retval File stream that reads the corresponding file.
 */
static inline optional<ifstream> select_file() {
  constexpr auto help = R"(
Choose what to do next.
h - Show this help list.
d - Show all sub-directories in current directory and change current working directory.
s - Select a file to read.
r - Read file with specified name.
c - Change current working directory.
p - Change current directory to parent directory if there is one.
q - quit
)";
  cout << help;
  while (true) {
    auto cp = current_path();
    auto di = directory_iterator(cp);

    cout << endl << cp.string() << "> ";
    char opt0 = 'q';
    cin >> opt0;
    switch (opt0) {
    case 'd': {
      auto dir = request_file_by_id(
          di, [](const directory_entry &entry) { return entry.is_directory(); },
          "sub-directory");
      if (dir.has_value())
        current_path(dir.value());
    } break;
    case 's': {
      auto bf = request_file_by_id(
          di,
          [](const directory_entry &entry) {
            return entry.path().extension() == ".bf";
          },
          ".bf file");

      if (bf.has_value()) {
        auto entry = bf.value();
        auto fin = ifstream(entry.path(), ios::binary);
        if (fin.is_open())
          return fin;
        cout << "Failed to open." << endl;
      }
    } break;
    case 'c': {
      auto dir = request_file_by_name(
          [](const path &p) { return is_directory(p); }, "directory");
      if (dir.has_value())
        current_path(dir.value());
    } break;
    case 'r': {
      auto opt = request_file_by_name(
          [](const path &p) { return is_regular_file(p); }, "directory");

      if (opt.has_value()) {
        auto f = opt.value();
        auto fin = ifstream(f, ios::binary);
        if (fin.is_open())
          return fin;
        cout << "Failed to open." << endl;
      }
    } break;
    case 'p': {
      current_path(current_path().parent_path());
    } break;
    case 'q': {
      return {};
    }
    default:
      cout << "Command not found." << endl;
    case 'h': {
      cout << help;
    } break;
    }
  }
}

static inline tuple<array<char, 30 + 1>, array<char, 30 + 1>,
                    array<char, 30 + 1>, vector<patch_info>, vector<frame>>
read_file_with_mode(istream &fin) {

  cout << R"(
h - Header only.
p - Header and patches.
f - Header, patches and all frames.
d - Discard.
)";
  char c = 'd';
  cin >> c;
  switch (c) {
  case 'h': {
    cout << "Reading started." << endl;
    auto &&[label, bar_label, units] = read_file_header(fin);
    cout << "Reading finished." << endl;
    return {std::move(label), std::move(bar_label), std::move(units), {}, {}};
  } break;
  case 'p': {
    cout << "Reading started." << endl;
    auto [label, bar_label, units, patches] = read_file_header_and_patches(fin);
    cout << "Reading finished." << endl;
    return {std::move(label),
            std::move(bar_label),
            std::move(units),
            std::move(patches),
            {}};
  } break;
  case 'f': {
    cout << "Reading started." << endl;
    auto [label, bar_label, units, patches, frames] = read_file(fin);
    cout << "Reading finished." << endl;
    return {std::move(label), std::move(bar_label), std::move(units),
            std::move(patches), std::move(frames)};
  } break;
  default:
    return {};
  }
}

/// @retval Whether to quit.
/// @param label Boundary quantity label.
/// @param bar_label Boundary quantity label shown on the bar.
/// @param units Display units.
/// @param patches Information of patches stored in the .bf file.
/// @param frames Information and boundary quantity data of frames stored in the
/// .bf file.
static inline void process_file() {
  tuple<array<char, 30 + 1>, array<char, 30 + 1>, array<char, 30 + 1>,
        vector<patch_info>, vector<frame>>
      data;
  auto &[label, bar_label, units, patches, frames] = data;

  tuple<vector<float>, vector<u32>, vector<u32>, vector<u32>> elements_data;
  auto &[nodes, sizes, elems, nums] = elements_data;

  auto elem_available = [&]() {
    auto &[nodes, sizes, elems, nums] = elements_data;
    return !nodes.empty() && !sizes.empty() && !elems.empty();
  };

  auto help = [&]() {
    auto &[label, bar_label, units, patches, frames] = data;
    auto &[nodes, sizes, elems, nums] = elements_data;

    cout <<
        R"(
Commands: 
q - Quit.
h - Show this help.
r - Read file. Will override the file that was read before if there is one.
e - Read elements.
S - Execute a script.
a - Analyze.
u - Unload file.)";
    if (label.front() && bar_label.front() && units.front())
      cout << R"(
b - Show boundary quantity basic information.)";

    if (!frames.empty())
      cout << R"(
f - Show frames.
F - Search for frame.)"
#if GRAPHICS_ENABLED
              R"(
V - Visualize frame.)"
#endif // GRAPHICS_ENABLED
          ;

    if (elem_available())
      cout << R"(
R - Visualize elements.
Y - Visualize elements polygons.)";

    if (!patches.empty())
      cout << R"(
p - Show patches.
P - Search for Patch.)"
#if GRAPHICS_ENABLED
              R"(
v - Visualize patches geometry.)"
#endif // GRAPHICS_ENABLED
          ;

    if (!patches.empty() && elem_available())
      cout << R"(
N - Visualize nodes and patches geometry.)";

    if (!patches.empty() && !frames.empty())
      cout << R"(
a - Analyze patch data.)";

    cout << endl;
  };
  help();
  while (true) {
    cout << "Please input your command here: ";
    char op1 = 'q';
    cin >> op1;
    switch (op1) {
    case 'q':
      return;
    case 'u': {
      label.front() = bar_label.front() = units.front() = '\0';
      patches.clear();
      frames.clear();
    } break;
    case 'e': {
      elements_data = read_nodes_and_elements();
      if (!elem_available()) {
        cout << "Element data not valid." << endl;
      }
    }; break;
    case 'r': {
      auto in = select_file();
      if (!in.has_value())
        break;
      auto &fin = in.value();
      fin.peek();
      if (!fin) {
        cout << "Failed to open the file." << endl;
      }
      data = read_file_with_mode(fin);
    } break;
#if GRAPHICS_ENABLED
    case 'N':
      if (!patches.empty() && elem_available()) {
        visualize_patches_and_elements(nodes, sizes, elems, patches);
      } else
        goto CMDNF;
    case 'R':
      if (elem_available()) {
        visualize_3d_elements(nodes, sizes, elems);
      } else
        goto CMDNF;
      break;
    case 'Y':
      if (elem_available()) {
        auto [polygon_sizes, polygon_indices, _0, _1] =
            get_polygon<false, false>(sizes, elems, {});
        visualize_polygons(nodes, polygon_sizes, polygon_indices);
      } else
        goto CMDNF;
      break;
#endif // GRAPHICS_ENABLED
    case 'b': {
      print_header(cout, label, bar_label, units);
    } break;
    case 'a': {
      analyze(patches, frames, nodes, sizes, elems, nums);
    } break;
    case 'f': {
      if (!frames.empty())
        print_frames(cout, frames);
      else
        goto CMDNF;
    } break;
#if GRAPHICS_ENABLED
    case 'v': {
      if (!patches.empty())
        visualize_patch(patches);
      else
        goto CMDNF;
    } break;
#endif // GRAPHICS_ENABLED
#if GRAPHICS_ENABLED
    case 'V': {
      if (!frames.empty()) {
        visualize_frames(patches, frames, get_data_category(units));
      } else
        goto CMDNF;
    } break;
#endif // GRAPHICS_ENABLED
    case 'p': {
      if (!patches.empty())
        print_patches(cout, patches);
      else
        goto CMDNF;
    } break;
    case 'F': {
      if (!frames.empty())
        search_frame(frames);
      else
        goto CMDNF;
    } break;
    case 'P': {
      if (!patches.empty())
        search_patch(patches);
      else
        goto CMDNF;
    } break;
    case 'S': {
      auto opt = request_file_by_name([](const path &p) { return exists(p); },
                                      "script");
      if (!opt)
        break;
      ifstream in(opt.value());
      string line;
      while (in.peek() != char_traits<char>::eof() && in) {
        getline(in, line);
        cin << line;
        cin.endl();
      }
    }; break;
    case 'h':
      help();
      break;
    default:
    CMDNF:
      cout << "Command not found." << endl;
      help();
      break;
    }
  }
}

int main(int argc, char *argv[]) {
  try {
    if (argc >= 2) {
      if (strcmp(argv[1], "--") == 0) {
        clog << "Forwarding command line arguments." << endl;
        for (int i = 2; i < argc; ++i) {
          cin << argv[i];
          cin.endl();
        }
      } else if (strcmp(argv[1], "-s") == 0) {
        ifstream fin(argv[2]);
        if (!fin) {
          cerr << "Failed to open the script." << endl;
          return -1;
        }
        clog << "Loading script." << endl;
        string file;
        getline(fin, file, (char)char_traits<char>::eof());
        cin << file;
        cin.endl();
      } else {
        clog << "Command line arguments ignored." << endl;
      }
    }

    cerr << hex;
    process_file();

  } catch (const exception &e) {
    cerr << "Exception caught: " << e.what() << endl;
  }

  return 0;
}