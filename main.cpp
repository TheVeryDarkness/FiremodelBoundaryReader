// Written by Boyi Huang,
// who was a student in Tongji University when he wrote this.
// C++17 is required for filesystem operation.
// GLAD and GLFW3 is required for boundary patch display.
// For reading .bf files generated by fds.
// See Source/dump.f90 in fds reposity for more information,
// especailly codes that are writing data to LU_BNDF(NF,NM).
// See also section Boundary Files in FDS_User_Guide.pdf.
// As for opengl, most of the codes about opengl are adapted from
// https://github.com/JoeyDeVries/LearnOpenGL.
#include <array>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <optional>
#include <vector>

#if !defined(GRAPHICS_ENABLED)
#define GRAPHICS_ENABLED 1
#endif // !defined(GRAPHICS_ENABLED)

#if GRAPHICS_ENABLED
#include <glad/glad.h>
#include <glfw/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#endif // GRAPHICS_ENABLED

using namespace std;

class file_error : exception {
public:
  const char *what() const noexcept {
    return "File is not correct. Maybe something unexpected happened when "
           "running fds.";
  }
};
class stream_error : exception {
public:
  const char *what() const noexcept {
    return "File stream terminated unexpectedly. This is usually because of an "
           "unexpected EOF.";
  }
};

#define CHECK_STREAM(STREAM)                                                   \
  if (!(STREAM))                                                               \
    throw stream_error();

#define CHECK_FORMAT(EXPR)                                                     \
  if (!(EXPR))                                                                 \
    throw file_error();

template <size_t sz> static inline void read(istream &in, array<char, sz> &s) {
  CHECK_STREAM(in);
  static_assert(sz > 0);
  in.get(s.data(), sz);
}

// Actually I'm not sure about their meanings.
constexpr const char string_separator[] = "\x1E\x00\x00\x00";
constexpr const char integer_separator[] = "\x04\x00\x00\x00";
constexpr const char line_separator[] = "\x24\x00\x00\x00";

template <size_t sz>
static inline void write_line(ostream &out, array<char, sz> &s) {
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

using i64 = std::int64_t;
using u64 = std::uint64_t;
using i32 = std::int32_t;
using u32 = std::uint32_t;
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
  friend ostream &operator<<(ostream &o, const patch_info &patch) {
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
};
struct patch_data {
  vector<float> data;
};
struct frame {
  float time;
  vector<patch_data> data;
};

/*
 * @retval Label, Bar Label, Units, Patches, Frames
 */
static inline tuple<array<char, 30 + 1>, array<char, 30 + 1>,
                    array<char, 30 + 1>, vector<patch_info>, vector<frame>>
read_file(istream &fin) {
  array<char, 30 + 1> label, bar_label, units;
  u32 n_patch;

  check(fin, string_separator);
  read(fin, label);
  check(fin, string_separator);
  check(fin, string_separator);
  read(fin, bar_label);
  check(fin, string_separator);
  check(fin, string_separator);
  read(fin, units);
  check(fin, string_separator);

  check(fin, integer_separator);
  n_patch = read_uint32(fin);
  check(fin, integer_separator);

  vector<patch_info> patches;
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
    const auto &patch = patches.back();
    check(fin, line_separator);
  }

  vector<frame> vals;
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
      for (size_t i = 0; i < info.size(); ++i) {
        float val = read_float(fin);
        current[ip].data.push_back(val);
      }
      u32 patch_end = read_uint32(fin);
      CHECK_FORMAT(patch_end == patch_size);
    }

    vals.push_back(frame{stime, std::move(current)});
  }
  return {label, bar_label, units, patches, vals};
}

static inline void search_frame_by_time(const vector<frame> &frames) {
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
      break;
    }
  }
}

static inline void search_frame_or_patch(const vector<frame> &frames,
                                         const vector<patch_info> &patches) {
  while (true) {
    cout << R"(
Select item to search for.
f - Frame.
p - Patch.
d - Discard.
)" << endl;
    char op2 = 'd';
    cin >> op2;
    switch (op2) {
    case 'f': {
      search_frame(frames);
    } break;
    case 'p': {
      search_patch(patches);
    } break;
    case 'd':
      return;
    default:
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
h - Show this list.
d - Show all sub-directories in current directory and change current working directory.
r - Select a file to read.
c - Change current working directory.
w - Show current working directory.
p - Change current directory to parent directory if there is one.
q - quit
)";
  cout << help;
  while (true) {
    auto cp = filesystem::current_path();
    auto di = filesystem::directory_iterator(cp);

    char opt = 'q';
    cin >> opt;
    switch (opt) {
    case 'h': {
      cout << help;
    } break;
    case 'd': {
      size_t i = 0;
      vector<filesystem::directory_entry> entries;
      for (const auto &entry : di)
        if (entry.is_directory()) {
          cout << i << " - " << entry.path().generic_string() << endl;
          entries.push_back(entry);
          ++i;
        }
      if (entries.empty()) {
        cout << "No sub-directory found." << endl;
        break;
      }
      cout << "Directory index(Invalid index to discard): ";
      cin >> i;
      if (i < entries.size())
        filesystem::current_path(entries[i]);
    } break;
    case 'r': {
      vector<filesystem::path> files;
      for (const auto &entry : di) {
        const auto &path = entry.path();
        if (path.extension() == ".bf") {
          files.push_back(path);
        }
      }
      cout << "There are " << files.size()
           << " .bf files in total in current directory." << endl;
      if (files.empty()) {
        cout << ".bf files not found in current directory." << endl;
        continue;
      }
      size_t i = 0;
      for (const auto &file : files) {
        cout << i << " - " << file.string() << endl;
        ++i;
      }
      cout << "File index: ";
      cin >> i;
      if (i >= files.size()) {
        cout << "Not a valid file index." << endl;
        continue;
      }
      auto fin = ifstream(files[i], ios::binary);
      if (fin.is_open())
        return fin;
      cout << "Failed to open." << endl;
      continue;
    } break;
    case 'c': {
      cout << "Directory: ";
      string s;
      cin >> s;
      filesystem::path p = s;
      if (filesystem::is_directory(p))
        filesystem::current_path(p);
      else
        cout << "Not a directory." << endl;
    } break;
    case 'w': {
      cout << filesystem::current_path() << endl;
    } break;
    case 'p': {
      filesystem::current_path(filesystem::current_path().parent_path());
    } break;
    case 'q': {
      return {};
    }
    default:
      break;
    }
  }
}

static inline void show_patch(const vector<patch_info> &patches,
                              const vector<frame> &frames) {
  u32 f = 0, p = 0;
  cout << "Select a frame and patch." << endl;
  cout << "Frame index: ";
  cin >> f;
  cout << "Patch index: ";
  cin >> p;
  u32 m = 0, n = 0;
  if (f >= frames.size()) {
    cout << "Not a valid frame." << endl;
    return;
  }
  if (p >= frames[f].data.size()) {
    cout << "Not a valid patch." << endl;
    return;
  }
  const auto &patch = frames[f].data[p];
  switch (patches[p].IOR) {
  case 1:
  case -1:
    m = patches[p].J();
    n = patches[p].K();
    break;
  case 2:
  case -2:
    m = patches[p].I();
    n = patches[p].K();
    break;
  case 3:
  case -3:
    m = patches[p].I();
    n = patches[p].J();
    break;
  default:
    CHECK_FORMAT(false);
  }

  u32 precision = 0;
  cout << "Precision:";
  cin >> precision;
  cout << setprecision(precision);

  u32 i = 0;
  for (auto v : patch.data) {
    cout << v;
    ++i;
    if (i % m == 0) {
      cout << endl;
    } else
      cout << ' ';
  }
  cout << endl;
}

#if GRAPHICS_ENABLED

static inline glm::vec3 cameraPos = glm::vec3(-1.0f, 0.0f, 0.0f);
static inline glm::vec3 cameraFront = glm::vec3(1.0f, 0.0f, 0.0f);
static inline glm::vec3 cameraUp = glm::vec3(0.0f, 0.0f, 1.0f);
/// @var time between current frame and last frame
static inline float deltaTime = 0.0f;
static inline float lastFrame = 0.0f;

#define STR(NUM) #NUM
#define DETECT_ERROR                                                           \
  do {                                                                         \
    auto gl_err = glGetError();                                                \
    if (gl_err != GL_NO_ERROR)                                                 \
      clog << "OpenGL Error " << gl_err << " detected at " __FILE__ "("        \
           << __LINE__ << ")" << endl;                                         \
                                                                               \
    const char *error = nullptr;                                               \
    auto glfw_err = glfwGetError(&error);                                      \
    if (glfw_err != GLFW_NO_ERROR)                                             \
      clog << "GLFW Error " << glfw_err                                        \
           << " detected at " __FILE__ "(" STR(__LINE__) ")" << endl;          \
    if (error && *error)                                                       \
      clog << *error << endl;                                                  \
  } while (0);

static bool cursor_enabled = false;
void processInput(GLFWwindow *window) {
  static bool tab_pressed = false;
  if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
    glfwSetWindowShouldClose(window, true);
  if (glfwGetKey(window, GLFW_KEY_TAB) == GLFW_PRESS) {
    if (!tab_pressed) {
      cursor_enabled = !cursor_enabled;
      glfwSetInputMode(window, GLFW_CURSOR,
                       cursor_enabled ? GLFW_CURSOR_NORMAL
                                      : GLFW_CURSOR_DISABLED);
    }
    tab_pressed = true;
  } else
    tab_pressed = false;
}

void onFramebufferSizeChange(GLFWwindow *window, int width, int height) {
  glViewport(0, 0, width, height);
}

void onMouseMove(GLFWwindow *window, double xposIn, double yposIn) {
  if (cursor_enabled)
    return;

  float xpos = static_cast<float>(xposIn);
  float ypos = static_cast<float>(yposIn);
  static bool firstMouse = true;
  static GLfloat lastX;
  static GLfloat lastY;

  if (firstMouse) {
    lastX = xpos;
    lastY = ypos;
    firstMouse = false;
  }

  float xoffset = xpos - lastX;
  /// reversed since y-coordinates go from bottom to top
  float yoffset = lastY - ypos;
  lastX = xpos;
  lastY = ypos;

  float sensitivity = 0.01f; // change this value to your liking
  xoffset *= sensitivity;
  yoffset *= sensitivity;

  cameraFront += glm::cross(cameraFront, cameraUp) * xoffset;
  cameraFront /= glm::length(cameraFront);

  auto up = cameraUp;
  cameraUp -= cameraFront * yoffset;
  cameraFront += up * yoffset;

  cameraUp /= glm::length(cameraUp);
  cameraFront /= glm::length(cameraFront);
}

void onScroll(GLFWwindow *window, double xoffset, double yoffset) {
  static GLfloat rate = 1;
  if (cursor_enabled) {
    rate *= (GLfloat)exp2(yoffset);
  } else {
    GLfloat dx = (GLfloat)xoffset * rate;
    GLfloat dy = (GLfloat)yoffset * rate;
    cameraPos += cameraFront * dy;
  }
}

auto vertexShaderSource =
    R"(
#version 330 core
layout (location = 0) in vec3 aPos;

uniform mat4 projection;
uniform mat4 view;

void main() {
 gl_Position = projection * view * vec4(aPos.x, aPos.y, aPos.z, 1.0);
}
)";
auto fragmentShaderSource = R"(
#version 330 core
out vec4 FragColor;
void main() {
 FragColor = vec4(1.0f, 1.0f, 1.0f, 0.6f);
}
)";

static inline int show_patch_position(const vector<patch_info> &patches) {
  bool wireframe = true;
  GLuint screenWidth = 800;
  GLuint screenHeight = 600;

  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

  GLFWwindow *window =
      glfwCreateWindow(screenWidth, screenHeight, "BoundaryReader", NULL, NULL);
  if (window == NULL) {
    cout << "Failed to create GLFW window." << endl;
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window);
  glfwSetFramebufferSizeCallback(window, onFramebufferSizeChange);
  glfwSetCursorPosCallback(window, onMouseMove);
  glfwSetScrollCallback(window, onScroll);

  glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    cout << "Failed to initialize GLAD" << endl;
    return -1;
  }

  glEnable(GL_DEPTH_TEST);

  DETECT_ERROR;

  int success;
  char infoLog[512];

  auto vertexShader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
  glCompileShader(vertexShader);

  glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
  if (!success) {
    glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
    cerr << infoLog << endl;
  }

  auto fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
  glCompileShader(fragmentShader);

  glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
  if (!success) {
    glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
    cerr << infoLog << endl;
  }

  DETECT_ERROR;

  // link shaders
  auto shaderProgram = glCreateProgram();
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  glLinkProgram(shaderProgram);

  // check for linking errors
  glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
  if (!success) {
    glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
    cerr << infoLog << endl;
  }
  glDeleteShader(vertexShader);
  glDeleteShader(fragmentShader);

  DETECT_ERROR;

  vector<GLuint> vertices = {};
  vector<GLuint> indices = {};
  for (const auto &patch : patches) {
    const GLuint N = (GLuint)vertices.size() / 3;
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

    vertices.push_back(patch.I1);
    vertices.push_back(patch.J1);
    vertices.push_back(patch.K1);

    vertices.push_back(patch.I2);
    vertices.push_back(patch.J1);
    vertices.push_back(patch.K2);

    vertices.push_back(patch.I2);
    vertices.push_back(patch.J2);
    vertices.push_back(patch.K1);

    vertices.push_back(patch.I1);
    vertices.push_back(patch.J2);
    vertices.push_back(patch.K2);
  }

  GLuint VBO, VAO, EBO;
  glGenVertexArrays(1, &VAO);
  glGenBuffers(1, &VBO);
  glGenBuffers(1, &EBO);

  glBindVertexArray(VAO);

  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER,
               vertices.size() * sizeof(decltype(vertices)::value_type),
               vertices.data(), GL_STATIC_DRAW);

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               indices.size() * sizeof(decltype(indices)::value_type),
               indices.data(), GL_STATIC_DRAW);

  DETECT_ERROR;

  glVertexAttribPointer(0, 3, GL_UNSIGNED_INT, GL_FALSE, 3 * sizeof(GLuint),
                        (void *)0);
  glEnableVertexAttribArray(0);

  DETECT_ERROR;

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glBindVertexArray(0);

  DETECT_ERROR;

  glUseProgram(shaderProgram);

  DETECT_ERROR;

  constexpr auto u32_max = numeric_limits<u32>::max();
  u32 I1 = u32_max, I2 = 0, J1 = u32_max, J2 = 0, K1 = u32_max, K2 = 0;
  for (const auto &patch : patches) {
    I1 = min(patch.I1, I1);
    I2 = max(patch.I2, I2);
    J1 = min(patch.I1, J1);
    J2 = max(patch.I1, J2);
    K1 = min(patch.I1, K1);
    K2 = max(patch.I1, K2);
  }
  float far = 100.;
  if (patches.size() > 0) {
    u32 I, J, K;
    I = I2 - I1;
    J = J2 - J1;
    K = K2 - K1;
    far = float(sqrt(I * I + J * J + K * K) * 2);
  }

  glm::mat4 projection = glm::perspective(
      glm::radians(60.0f), float(screenWidth) / screenHeight, 0.1f, far);

  glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "projection"), 1,
                     false, glm::value_ptr(projection));

  DETECT_ERROR;

  // Render loop
  while (!glfwWindowShouldClose(window)) {
    // camera/view transformation
    glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "view"), 1, false,
                       glm::value_ptr(view));

    DETECT_ERROR;

    processInput(window);

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(shaderProgram);
    glBindVertexArray(VAO);
    if (wireframe)
      glDrawElements(GL_LINES, (GLuint)indices.size(), GL_UNSIGNED_INT, 0);
    else
      glDrawElements(GL_TRIANGLES, (GLuint)indices.size(), GL_UNSIGNED_INT, 0);

    glfwSwapBuffers(window);
    glfwPollEvents();

    DETECT_ERROR;
  }

  DETECT_ERROR;

  glDeleteVertexArrays(1, &VAO);
  glDeleteBuffers(1, &VBO);
  glDeleteProgram(shaderProgram);

  glfwTerminate();
  return 0;
}
#endif // GRAPHICS_ENABLED

/// @retval Whether to quit.
static inline bool process_file(array<char, 30 + 1> &label,
                                array<char, 30 + 1> &bar_label,
                                array<char, 30 + 1> &units,
                                const vector<patch_info> &patches,
                                const vector<frame> &frames) {
  constexpr auto help =
      R"(
Input to get information.
q - Quit.
h - Show this help information.
b - Show boundary quantity information.
f - Show frames information.
p - Show patches information.
o - Show patch data.)"
#if GRAPHICS_ENABLED
      R"(
g - Show patch position in graphics.)"
#endif // GRAPHICS_ENABLED
      R"(
s - Search for frame or patch.
u - Unload file.
)";
  cout << help << endl;
  while (true) {
    cout << "Please input your command here: ";
    char op1 = 'q';
    cin >> op1;
    switch (op1) {
    case 'q':
      return true;
      break;
    case 'u':
      return false;
    case 'b':
      write_line(cout << "Label:     ", label);
      write_line(cout << "Bar Label: ", bar_label);
      write_line(cout << "Units:     ", units);
      break;
    case 'o': {
      show_patch(patches, frames);
    } break;
    case 'f': {
      size_t i = 0;
      for (const auto &f : frames) {
        cout << "Frame " << i << " at " << f.time << "s." << endl;
        ++i;
      }
    } break;
#if GRAPHICS_ENABLED
    case 'g':
      show_patch_position(patches);
      break;
#endif // GRAPHICS_ENABLED
    case 'p': {
      constexpr auto int_len = numeric_limits<i32>::digits10 + 1;
      constexpr auto uint_len = numeric_limits<u32>::digits10 + 1;
      cout << setw(uint_len) << "I1"
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
        cout << patch << endl;
    } break;
    case 's': {
      search_frame_or_patch(frames, patches);
    } break;
    default:
      cout << "Command not found." << endl;
    case 'h':
      cout << help << endl;
      break;
    }
  }
}

int main(int argc, char *argv[]) {
  while (true) {
    try {
      cerr << hex;

      auto in = select_file();
      if (!in.has_value())
        return 1;
      auto &fin = in.value();
      fin.peek();
      if (!fin) {
        cout << "Failed to open the file." << endl;
      }

      cout << "Reading started." << endl;

      auto [label, bar_label, units, patches, frames] = read_file(fin);

      cout << "Reading finished." << endl;

      bool quitting = process_file(label, bar_label, units, patches, frames);
      if (quitting)
        return 0;

    } catch (const exception &e) {
      cerr << "Exception caught: " << e.what() << endl;
    }
  }

  return 0;
}