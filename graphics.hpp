#pragma once
#include "fds_basic.hpp"
#include "types.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>
using std::accumulate;
using std::cerr;
using std::clog;
using std::cout;
using std::endl;
using std::extent_v;
using std::get;
using std::index_sequence;
using std::initializer_list;
using std::integral_constant;
using std::is_invocable_r_v;
using std::is_same_v;
using std::make_index_sequence;
using std::make_tuple;
using std::map;
using std::max;
using std::max_element;
using std::min;
using std::numeric_limits;
using std::remove_all_extents_t;
using std::set;
using std::tuple_element_t;
using std::underlying_type_t;

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
#include "element_basic.hpp"
#include "shared.hpp"

#if GRAPHICS_ENABLED

template <typename> struct gl_type_enum;
template <>
struct gl_type_enum<GLuint> : integral_constant<GLenum, GL_UNSIGNED_INT> {};
template <> struct gl_type_enum<GLint> : integral_constant<GLenum, GL_INT> {};
template <>
struct gl_type_enum<GLfloat> : integral_constant<GLenum, GL_FLOAT> {};
template <>
struct gl_type_enum<GLdouble> : integral_constant<GLenum, GL_DOUBLE> {};

template <typename T>
constexpr static inline GLenum gl_type_enum_v = gl_type_enum<T>::value;

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

static inline bool cursor_enabled = false;
static inline bool patch_loop = false;
static inline size_t current = 0;
static inline size_t patches_count = 0;
static inline float key_move_sensity = 0.1f;

static inline void onKey(GLFWwindow *window, const int key, int scancode,
                         int action, const int mods) {
  if (key == GLFW_KEY_ESCAPE)
    if (action == GLFW_PRESS)
      glfwSetWindowShouldClose(window, true);

  if (key == GLFW_KEY_TAB)
    if (action == GLFW_PRESS) {
      cursor_enabled = !cursor_enabled;
      glfwSetInputMode(window, GLFW_CURSOR,
                       cursor_enabled ? GLFW_CURSOR_NORMAL
                                      : GLFW_CURSOR_DISABLED);
    }

  if (action == GLFW_REPEAT || action == GLFW_PRESS) {
    if (key == GLFW_KEY_LEFT) {
      cameraPos += glm::cross(cameraUp, cameraFront) * key_move_sensity;
    }
    if (key == GLFW_KEY_RIGHT) {
      cameraPos += glm::cross(cameraFront, cameraUp) * key_move_sensity;
    }
    if (key == GLFW_KEY_UP) {
      cameraPos += cameraFront * key_move_sensity;
    }
    if (key == GLFW_KEY_DOWN) {
      cameraPos -= cameraFront * key_move_sensity;
    }
  }

  const bool continuous = mods & GLFW_MOD_CAPS_LOCK;

  if (key == GLFW_KEY_LEFT_BRACKET)
    if (action == GLFW_PRESS || (continuous && action == GLFW_REPEAT)) {
      if (current > 0)
        --current;
      else if (patch_loop)
        current = patches_count ? patches_count - 1 : 0;
    }

  if (key == GLFW_KEY_RIGHT_BRACKET)
    if (action == GLFW_PRESS || (continuous && action == GLFW_REPEAT)) {
      if (current < patches_count) {
        ++current;
      } else if (patch_loop)
        current = 0;
    }

  string title = std::to_string(current);
  glfwSetWindowTitle(window, title.c_str());
}

static inline GLfloat sensitivity = 0.01f;
static inline void onMouseMove(GLFWwindow *window, double xposIn,
                               double yposIn) {

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
  /// Reversed since y-coordinates go from bottom to top
  float yoffset = lastY - ypos;
  lastX = xpos;
  lastY = ypos;

  if (cursor_enabled)
    return;

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

static inline GLfloat rate = 1;
static inline void onScroll(GLFWwindow *window, double xoffset,
                            double yoffset) {
  if (cursor_enabled) {
    rate *= (GLfloat)exp2(yoffset);
  } else {
    GLfloat dx = (GLfloat)xoffset * rate;
    GLfloat dy = (GLfloat)yoffset * rate;
    cameraPos += cameraFront * dy;
  }
}

static GLuint windowWidth = 1600;
static GLuint windowHeight = 1200;
static inline void onFramebufferSizeChange(GLFWwindow *window, int width,
                                           int height) {
  windowWidth = width;
  windowHeight = height;
  glViewport(0, 0, width, height);
}

static inline void onGLDebugMessage(GLenum source, GLenum type, unsigned int id,
                                    GLenum severity, GLsizei length,
                                    const char *message,
                                    const void *userParam) {
  clog << "Debug message (" << id << "): " << message << endl;

  switch (source) {
  case GL_DEBUG_SOURCE_API:
    clog << "Source: API";
    break;
  case GL_DEBUG_SOURCE_WINDOW_SYSTEM:
    clog << "Source: Window System";
    break;
  case GL_DEBUG_SOURCE_SHADER_COMPILER:
    clog << "Source: Shader Compiler";
    break;
  case GL_DEBUG_SOURCE_THIRD_PARTY:
    clog << "Source: Third Party";
    break;
  case GL_DEBUG_SOURCE_APPLICATION:
    clog << "Source: Application";
    break;
  case GL_DEBUG_SOURCE_OTHER:
    clog << "Source: Other";
    break;
  }
  std::cout << std::endl;

  switch (type) {
  case GL_DEBUG_TYPE_ERROR:
    clog << "Type: Error";
    break;
  case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR:
    clog << "Type: Deprecated Behaviour";
    break;
  case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR:
    clog << "Type: Undefined Behaviour";
    break;
  case GL_DEBUG_TYPE_PORTABILITY:
    clog << "Type: Portability";
    break;
  case GL_DEBUG_TYPE_PERFORMANCE:
    clog << "Type: Performance";
    break;
  case GL_DEBUG_TYPE_MARKER:
    clog << "Type: Marker";
    break;
  case GL_DEBUG_TYPE_PUSH_GROUP:
    clog << "Type: Push Group";
    break;
  case GL_DEBUG_TYPE_POP_GROUP:
    clog << "Type: Pop Group";
    break;
  case GL_DEBUG_TYPE_OTHER:
    clog << "Type: Other";
    break;
  }
  clog << endl;

  switch (severity) {
  case GL_DEBUG_SEVERITY_HIGH:
    clog << "Severity: high";
    break;
  case GL_DEBUG_SEVERITY_MEDIUM:
    clog << "Severity: medium";
    break;
  case GL_DEBUG_SEVERITY_LOW:
    clog << "Severity: low";
    break;
  case GL_DEBUG_SEVERITY_NOTIFICATION:
    clog << "Severity: notification";
    break;
  }
  clog << endl << endl;
}

static inline GLuint
create_shader_program(initializer_list<const char *> vertexShaderSource,
                      initializer_list<const char *> fragmentShaderSource) {
  int success;
  char infoLog[512];

  DETECT_ERROR;

  auto vertexShader = glCreateShader(GL_VERTEX_SHADER);
  glShaderSource(vertexShader, (GLsizei)vertexShaderSource.size(),
                 vertexShaderSource.begin(), NULL);
  glCompileShader(vertexShader);

  glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
  if (!success) {
    glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
    cerr << infoLog << endl;

#ifndef NDEBUG
    cerr << "Vertex Shader: " << endl;
    for (auto s : vertexShaderSource) {
      cerr << s;
    }
    cerr << "Fragment Shader: " << endl;
    for (auto s : fragmentShaderSource) {
      cerr << s;
    }
#endif // !NDEBUG
  }

  auto fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(fragmentShader, (GLsizei)fragmentShaderSource.size(),
                 fragmentShaderSource.begin(), NULL);
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

  return shaderProgram;
}

static bool fullScreen = false;
static bool wireframe = false;

/// @retval Start visualization or discard.
static inline bool visualization_settings() {
  while (true) {
    cout << "Settings:" << endl
         << "f - Display as wireframe: " << wireframe << endl
         << "c - Current patch:        " << current << endl
         << "m - Mouse move sensity:   " << rate << endl
         << "s - Scroll sensity:       " << sensitivity << endl
         << "F - Fullscreen:           " << fullScreen << endl
         << "W - Window width:         " << windowWidth << endl
         << "H - Window height:        " << windowHeight << endl
         << "C - Cell size:            " << mesh.cell_size << endl
         << "X - Origin x:             " << mesh.x0 << endl
         << "Y - Origin y:             " << mesh.y0 << endl
         << "Z - Origin z:             " << mesh.z0 << endl
         << "Options:" << endl
         << "h - Show graphics help." << endl
         << "g - Start patch visualiztion." << endl
         << "d - Discard." << endl;
    char opt = 'd';
    cin >> opt;
    switch (opt) {
    case 'f': {
      cout << "Display as wireframe: ";
      cin >> wireframe;
    } break;
    case 'F': {
      cout << "Display as full screen mode: ";
      cin >> fullScreen;
    } break;
    case 'W': {
      cout << "Window width: ";
      cin >> windowWidth;
    } break;
    case 'H': {
      cout << "Window height: ";
      cin >> windowHeight;
    } break;
    case 'c': {
      cout << "Current patch: ";
      cin >> current;
    } break;
    case 's': {
      cout << "Scroll sensity: ";
      cin >> rate;
    } break;
    case 'm': {
      cout << "Mouse move sensity: ";
      cin >> sensitivity;
    } break;
    case 'C': {
      cin >> mesh.cell_size;
    } break;
    case 'X': {
      cin >> mesh.x0;
    } break;
    case 'Y': {
      cin >> mesh.y0;
    } break;
    case 'Z': {
      cin >> mesh.z0;
    } break;
    case 'h': {
      cout << R"(
Press TAB to enable or disable cursor.
Press CAPS to enable or disable continuous key input.
Press [ or ] to change current highlighted patch.
Press ESC to terminate.
Move mouse to rotate the camera when cursor is disabled.
Direction key to move.
Scroll to move forward or backward when cursor is disabled.
Scroll to modify the scroll sensity when cursor is enabled.
)";
    } break;
    case 'g': {
      return true;
    } break;
    case 'd':
      return false;
    default:
      cout << "Option not found." << endl;
      break;
    }
  }
}

template <typename T, size_t i, typename U>
static inline void bind_attribute(const vector<U> &vec, const GLuint VBO) {
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, vec.size() * sizeof(U), vec.data(),
               GL_STATIC_DRAW);

  static_assert(is_same_v<remove_all_extents_t<T>, U>, "Not matched.");
  constexpr auto extent = extent_v<T, 0>;
  constexpr auto length = extent == 0 ? 1 : extent;

  glVertexAttribPointer(i, length, gl_type_enum_v<U>, GL_FALSE,
                        length * sizeof(U), 0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  DETECT_ERROR;
}

template <typename... types> struct type_list {};

template <typename... T, typename... U, size_t VBOCount, size_t... i>
static inline void bind_attributes_sequence(const tuple<vector<U>...> &vec,
                                            const GLuint (&VBOs)[VBOCount],
                                            type_list<T...>,
                                            index_sequence<i...>) {
  using t = tuple<T...>;
  static_assert(((i < VBOCount) && ...), "Out of bound.");
  (bind_attribute<tuple_element_t<i, t>, i>(get<i>(vec), VBOs[i]), ...);
}

template <typename... T, typename... U, size_t VBOCount>
static inline void bind_attributes(const tuple<vector<U>...> &vec,
                                   const GLuint (&VBOs)[VBOCount],
                                   type_list<T...> Ts) {
  static_assert(sizeof...(T) == VBOCount);
  static_assert(sizeof...(U) == VBOCount);
  bind_attributes_sequence(vec, VBOs, Ts, make_index_sequence<VBOCount>{});
}

template <typename Index, typename GetData, typename GetNearFar,
          typename GetDrawMode, typename... Vertex>
static inline int visualize(GetData &&getData, GetNearFar &&getNearFar,
                            GetDrawMode &&getDrawMode,
                            initializer_list<const char *> vertexShaderSource,
                            initializer_list<const char *> fragmentShaderSource,
                            type_list<Vertex...>) {

  static_assert(
      is_invocable_r_v<
          tuple<tuple<vector<remove_all_extents_t<Vertex>>...>, vector<Index>>,
          GetData>,
      "Can't get data.");

  constexpr static auto vertexAttributesCount = sizeof...(Vertex);

  glfwInit();
  glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
  glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
  glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
  glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

#ifndef NDEBUG
  glfwWindowHint(GLFW_OPENGL_DEBUG_CONTEXT, GL_TRUE);
#endif // !NDEBUG

  GLFWmonitor *monitor = fullScreen ? glfwGetPrimaryMonitor() : NULL;
  if (fullScreen && !monitor)
    clog << "Full screen failed." << endl;

  GLFWwindow *window = glfwCreateWindow(windowWidth, windowHeight,
                                        "BoundaryReader", monitor, NULL);
  if (window == NULL) {
    cout << "Failed to create GLFW window." << endl;
    glfwTerminate();
    return -1;
  }
  glfwMakeContextCurrent(window);
  glfwSetInputMode(window, GLFW_LOCK_KEY_MODS, GLFW_TRUE);
  glfwSetKeyCallback(window, onKey);
  glfwSetFramebufferSizeCallback(window, onFramebufferSizeChange);
  glfwSetCursorPosCallback(window, onMouseMove);
  glfwSetScrollCallback(window, onScroll);

  glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

  if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
    cout << "Failed to initialize GLAD." << endl;
    return -1;
  }

#ifndef NDEBUG
  int flags;
  glGetIntegerv(GL_CONTEXT_FLAGS, &flags);
  if (flags & GL_CONTEXT_FLAG_DEBUG_BIT) {
    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    glDebugMessageCallback(onGLDebugMessage, nullptr);
    glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DONT_CARE, 0, nullptr,
                          GL_TRUE);
  }
#endif // !NDEBUG

  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  DETECT_ERROR;

  GLuint shaderProgram =
      create_shader_program(vertexShaderSource, fragmentShaderSource);

  DETECT_ERROR;

  auto [vertices, indices] = getData();

  GLuint VBO[vertexAttributesCount], VAO, EBO;
  glGenVertexArrays(1, &VAO);
  glGenBuffers(vertexAttributesCount, VBO);
  glGenBuffers(1, &EBO);

  DETECT_ERROR;

  glBindVertexArray(VAO);

  DETECT_ERROR;

  bind_attributes(vertices, VBO, type_list<Vertex...>{});

  DETECT_ERROR;

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(Index),
               indices.data(), GL_STATIC_DRAW);

  DETECT_ERROR;

  for (size_t i = 0; i < vertexAttributesCount; ++i)
    glEnableVertexAttribArray(i);

  DETECT_ERROR;

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindVertexArray(0);

  DETECT_ERROR;

  glUseProgram(shaderProgram);

  DETECT_ERROR;

  auto [near, far] = getNearFar();

  glm::mat4 projection = glm::perspective(
      glm::radians(60.0f), float(windowWidth) / windowHeight, near, far);

  DETECT_ERROR;

  auto mode = getDrawMode();

  // Render loop
  while (!glfwWindowShouldClose(window)) {
    // camera/view transformation
    glm::mat4 view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);

    glm::mat4 pv = projection * view;
    // Projection matrix multiplied with view matrix
    glUniformMatrix4fv(glGetUniformLocation(shaderProgram, "pv"), 1, false,
                       glm::value_ptr(pv));

    DETECT_ERROR;

    glUniform1ui(glGetUniformLocation(shaderProgram, "highlighted"), current);

    DETECT_ERROR;

    glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glUseProgram(shaderProgram);
    glBindVertexArray(VAO);
    glDrawElements(mode, (GLuint)indices.size(), gl_type_enum_v<Index>, 0);
    glBindVertexArray(0);

    glfwSwapBuffers(window);
    // glfwPollEvents();
    glfwWaitEvents();

    DETECT_ERROR;
  }

  glBindBuffer(GL_ARRAY_BUFFER, 0);
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  DETECT_ERROR;

  glDeleteVertexArrays(1, &VAO);
  glDeleteBuffers(vertexAttributesCount, VBO);
  glDeleteProgram(shaderProgram);

  glfwTerminate();
  return 0;
}

constexpr static inline struct VertexShaderSource {
  const char *head_pos_pv_color = R"(
#version 330 core
layout (location = 0) in vec3 aPos;
uniform mat4 pv;
out vec4 color;
)";
  const char *index = R"(
layout(location = 1) in float index;
)";
  const char *temperature = R"(
layout(location = 1) in float index;
)";
  const char *highlight = R"(
uniform uint highlighted;
)";
  const char *main_begin = R"(
void main() {
 gl_Position = pv * vec4(aPos.x, aPos.y, aPos.z, 1.0);
)";
  const char *main_wave_length_to_rgb = R"(
 float rr = 2800;
 float rg = 2000;
 float rb = 1200;
 float ti = 2;
 float wave_length = 2897772.1 / temperature;
 color = vec4(
  0.41 * pow(2, -pow(abs(wave_length - 595), ti) / rr),
  0.82 * pow(2, -pow(abs(wave_length - 530), ti) / rg),
  0.40 * pow(2, -pow(abs(wave_length - 460), ti) / rb),
  .6f
 );
)";
  const char *main_highlight = R"(
 if (highlighted == index) {
  color = vec4(.8f, .1f, .0f, .8f);
 } else {
  color = vec4(.8f, .8f, .8f, .1f);
 }
)";
  const char *main_color = R"(
  color = vec4(.6f, .6f, .6f, .1f);
)";
  const char *main_end = R"(
}
)";
} vertexShaderSource;

constexpr static inline auto &fragmentShaderSource = R"(
#version 330 core

in vec4 color;

void main() {
 gl_FragColor = color;
}
)";

static inline GLenum defaultDrawMode() {
  return wireframe ? GL_LINES : GL_TRIANGLES;
}

static inline int visualize_patch(const vector<patch_info> &patches) {

  patches_count = patches.size();
  if (current > patches.size()) {
    current = 0;
  }

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint>(
        [&patches]() { return from_patches<true>(patches, wireframe); },
        [&patches]() -> tuple<float, float> {
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
            far = float(sqrt(I * I + J * J + K * K) * 5);
          }

          return {0.01f, far};
        },
        defaultDrawMode,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.index,
         vertexShaderSource.highlight, vertexShaderSource.main_begin,
         vertexShaderSource.main_highlight, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<GLuint[3], GLuint>{});
  }
}

static inline int visualize_data(const vector<float> &data, u32 m, u32 n) {
  assert(data.size() == m * n);
  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint>(
        [&data, n]() { return from_matrix_data(data, n, wireframe); },
        [&data, m, n]() -> tuple<float, float> {
          return {0.01f, 5 * max(initializer_list<float>{
                                 *max_element(data.cbegin(), data.cend()),
                                 (float)m, (float)n})};
        },
        defaultDrawMode,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_3d_elements(const vector<float> &nodes,
                                        const vector<u32> &vertex_count,
                                        const vector<u32> &elements) {

  u32 sum = accumulate(vertex_count.cbegin(), vertex_count.cend(), 0);
  assert(elements.size() == sum);
  float _max =
      max(initializer_list<float>{*max_element(nodes.cbegin(), nodes.cend())});
  float ratio = _max == 0 ? 1 : 100 / _max;

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint>(
        [&nodes, &vertex_count, &elements, sum, ratio]() {
          return from_elements(nodes, vertex_count, elements, wireframe, sum,
                               ratio);
        },
        [&nodes]() constexpr->tuple<float, float> {
          return {0.01f, 500.f};
        },
        defaultDrawMode,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_patches_and_elements(
    const vector<float> &nodes, const vector<u32> &vertex_count,
    const vector<u32> &elements, const vector<patch_info> &patches) {
  u32 sum = accumulate(vertex_count.cbegin(), vertex_count.cend(), 0);
  assert(elements.size() == sum);
  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint>(
        [&nodes, &vertex_count, &elements, &patches, sum]() {
          return from_patches_and_elements(nodes, vertex_count, elements,
                                           patches, wireframe, sum);
        },
        [&nodes]() constexpr->tuple<float, float> {
          return {0.01f, 500.f};
        },
        defaultDrawMode,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_nodes(const vector<float> &nodes,
                                  const vector<u32> &indices) {
  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint>(
        [&nodes, &indices ]() -> auto{
          glPointSize(10);
          return make_tuple(make_tuple(nodes), indices);
        },
        [&nodes]() constexpr->tuple<float, float> {
          return {0.01f, 500.f};
        },
        []() { return GL_POINTS; },
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_primitives_on_patch(const vector<float> &nodes,
                                                const vector<u32> &vertex_count,
                                                const vector<u32> &elements,
                                                set<u32> nodes_on_patch) {
  u32 sum = accumulate(vertex_count.cbegin(), vertex_count.cend(), 0);
  assert(elements.size() == sum);
  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint>(
        [&nodes, &nodes_on_patch, &vertex_count, &elements, sum ]() -> auto{
          auto t = from_elements_if_in_set(nodes, vertex_count, elements,
                                           wireframe, sum, nodes_on_patch);
          return t;
        },
        [&nodes]() constexpr->tuple<float, float> {
          return {0.01f, 500.f};
        },
        defaultDrawMode,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

#endif // GRAPHICS_ENABLED
