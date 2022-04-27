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

template <typename... types> struct type_list {};

template <size_t i, typename... types> struct subset_start_from;
template <size_t i, typename first, typename... res>
struct subset_start_from<i, first, res...> {
  using type = typename subset_start_from<i - 1, res...>::type;
};
template <typename... types> struct subset_start_from<0, types...> {
  using type = type_list<types...>;
};

template <size_t i, typename... types>
using subset_start_from_t = typename subset_start_from<i, types...>::type;

static inline glm::vec3 cameraPos = glm::vec3(-1.0f, 0.0f, 0.0f);
static inline glm::vec3 cameraFront = glm::vec3(1.0f, 0.0f, 0.0f);
static inline glm::vec3 cameraUp = glm::vec3(0.0f, 0.0f, 1.0f);
/// @var time between current frame and last frame
static inline double deltaTime = 0.0;
static inline double lastFrame = 0.0;
static inline double firstFrame = 0.0;

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
static inline bool index_loop = false;
static inline size_t current = 0;
static inline size_t index_max = 0;
static inline float key_move_sensity = 0.5f;

static inline bool keyX[2] = {};
static inline bool keyY[2] = {};

static inline void cameraMoveForward(float length) {
  cameraPos += cameraFront * length;
}
static inline void cameraMoveBackward(float length) {
  cameraPos -= cameraFront * length;
}
static inline void cameraMoveLeft(float length) {
  cameraPos += glm::cross(cameraUp, cameraFront) * length;
}
static inline void cameraMoveRight(float length) {
  cameraPos += glm::cross(cameraFront, cameraUp) * length;
}
static inline void cameraRotateX(float radian) {
  cameraFront += glm::cross(cameraFront, cameraUp) * radian;
}
static inline void cameraRotateY(float radian) {
  cameraFront += cameraUp * radian;
}

static inline void keyCameraMove(float deltaTime) {
  auto length = deltaTime * key_move_sensity;
  if (keyX[1] && !keyX[0])
    cameraMoveRight(length);
  if (keyX[0] && !keyX[1])
    cameraMoveLeft(length);
  if (keyY[1] && !keyY[0])
    cameraMoveForward(length);
  if (keyX[0] && !keyX[1])
    cameraMoveBackward(length);
}

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

  if (key == GLFW_KEY_LEFT)
    keyX[0] = action != GLFW_RELEASE;
  if (key == GLFW_KEY_RIGHT)
    keyX[1] = action != GLFW_RELEASE;
  if (key == GLFW_KEY_UP)
    keyY[1] = action != GLFW_RELEASE;
  if (key == GLFW_KEY_DOWN)
    keyY[0] = action != GLFW_RELEASE;

  const bool continuous = mods & GLFW_MOD_CAPS_LOCK;

  if (key == GLFW_KEY_LEFT_BRACKET)
    if (action == GLFW_PRESS || (continuous && action == GLFW_REPEAT)) {
      if (current > 0)
        --current;
      else if (index_loop)
        current = index_max ? index_max - 1 : 0;
    }

  if (key == GLFW_KEY_RIGHT_BRACKET)
    if (action == GLFW_PRESS || (continuous && action == GLFW_REPEAT)) {
      if (current < index_max) {
        ++current;
      } else if (index_loop)
        current = 0;
    }

  string title = std::to_string(current);
  glfwSetWindowTitle(window, title.c_str());
}

static inline GLfloat sensitivity = 0.02f;
static bool firstMouse = true;
static inline void onMouseMove(GLFWwindow *window, double xposIn,
                               double yposIn) {

  float xpos = static_cast<float>(xposIn);
  float ypos = static_cast<float>(yposIn);
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

static inline GLfloat rate = 4;
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
         << "m - Key move sensity:     " << key_move_sensity << endl
         << "m - Mouse move sensity:   " << sensitivity << endl
         << "s - Scroll sensity:       " << rate << endl
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
    case 'k': {
      cout << "Key move sensity: ";
      cin >> key_move_sensity;
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

template <typename Ty>
constexpr static inline size_t extent =
    extent_v<Ty, 0> == 0 ? 1 : extent_v<Ty, 0>;

template <typename T, size_t i, typename U>
static inline void bind_attribute(const vector<U> &vec, const GLuint VBO) {
  glBindBuffer(GL_ARRAY_BUFFER, VBO);
  glBufferData(GL_ARRAY_BUFFER, vec.size() * sizeof(U), vec.data(),
               GL_STATIC_DRAW);

  static_assert(is_same_v<remove_all_extents_t<T>, U>, "Not matched.");
  constexpr auto length = extent<T>;

  glVertexAttribPointer(i, length, gl_type_enum_v<U>, GL_FALSE,
                        length * sizeof(U), 0);

  glBindBuffer(GL_ARRAY_BUFFER, 0);

  DETECT_ERROR;
}

template <size_t offset, typename... T, typename... U, size_t VBOCount,
          size_t... i>
static inline void bind_attributes_sequence(const tuple<vector<U>...> &vec,
                                            const GLuint (&VBOs)[VBOCount],
                                            type_list<T...>,
                                            index_sequence<i...>) {
  using t = tuple<T...>;
  static_assert(((i < VBOCount) && ...), "Out of bound.");
  static_assert(((i < sizeof...(U)) && ...), "Out of bound.");
  (bind_attribute<tuple_element_t<i, t>, i>(get<i>(vec), VBOs[offset + i]),
   ...);
}

template <size_t offset, typename... T, typename... U, size_t VBOCount>
static inline void bind_attributes(const tuple<vector<U>...> &vec,
                                   const GLuint (&VBOs)[VBOCount],
                                   type_list<T...> Ts) {
  static_assert(offset + sizeof...(T) == VBOCount);
  static_assert(offset + sizeof...(U) == VBOCount);
  bind_attributes_sequence<offset>(vec, VBOs, Ts,
                                   make_index_sequence<VBOCount - offset>{});
}

template <typename Index, size_t refreshData, typename GetData,
          typename GetRefreshedData, typename GetNearFar, typename GetDrawMode,
          typename GetMinMax, typename GetDomain, typename... Vertex>
static inline int visualize(GetData &&getData,
                            GetRefreshedData &&getRefreshedData,
                            GetNearFar &&getNearFar, GetDrawMode &&getDrawMode,
                            GetMinMax &&getMinMax, GetDomain &&getDomain,
                            initializer_list<const char *> vertexShaderSource,
                            initializer_list<const char *> fragmentShaderSource,
                            type_list<Vertex...>) {

  static_assert(
      is_invocable_r_v<
          tuple<tuple<vector<remove_all_extents_t<Vertex>>...>, vector<Index>>,
          GetData>,
      "Can't get data.");

  constexpr static auto vertexAttributesCount = sizeof...(Vertex);
  static_assert(refreshData <= vertexAttributesCount);

  firstMouse = true;

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

  glDisable(GL_CULL_FACE);
  glDisable(GL_DEPTH_TEST);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);
  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

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

  bind_attributes<0>(vertices, VBO, type_list<Vertex...>{});

  DETECT_ERROR;

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(Index),
               indices.data(), GL_STATIC_DRAW);

#ifndef NDEBUG
  size_t out = 0;
  size_t vc = std::get<0>(vertices).size() /
              extent<std::tuple_element_t<0, tuple<Vertex...>>>;
  for (auto i : indices)
    if (i > vc)
      ++out;
  if (out > 0)
    cerr << out << " vertices out of nodes count detected." << endl;
#endif // !NDEBUG

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

  if constexpr (!is_same_v<GetMinMax, nullptr_t>) {
    const auto [min, max] = getMinMax();

    static_assert(is_same_v<decay_t<decltype(min)>, float>);
    static_assert(is_same_v<decay_t<decltype(max)>, float>);

    glUniform1f(glGetUniformLocation(shaderProgram, "min"), min);

    DETECT_ERROR;

    glUniform1f(glGetUniformLocation(shaderProgram, "max"), max);

    DETECT_ERROR;
  }

  if constexpr (!is_same_v<GetDomain, nullptr_t>) {
    const auto [x, y, z] = getDomain();

    static_assert(is_same_v<decay_t<decltype(x)>, float>);
    static_assert(is_same_v<decay_t<decltype(y)>, float>);
    static_assert(is_same_v<decay_t<decltype(z)>, float>);

    glUniform1f(glGetUniformLocation(shaderProgram, "X"), x);

    DETECT_ERROR;

    glUniform1f(glGetUniformLocation(shaderProgram, "Y"), y);

    DETECT_ERROR;

    glUniform1f(glGetUniformLocation(shaderProgram, "Z"), z);

    DETECT_ERROR;
  }

  firstFrame = lastFrame = glfwGetTime();

  // Render loop
  while (!glfwWindowShouldClose(window)) {
    double currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    keyCameraMove((double)deltaTime);

    if constexpr (refreshData < vertexAttributesCount) {
      auto refreshedData = getRefreshedData();
      bind_attributes<refreshData>(
          make_tuple(refreshedData), VBO,
          subset_start_from_t<refreshData, Vertex...>{});
    }

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
  const char *domain = R"(
uniform float X;
uniform float Y;
uniform float Z;
)";
  const char *data = R"(
layout(location = 1) in float data;
uniform float min;
uniform float max;
)";
  const char *highlight = R"(
uniform uint highlighted;
)";
  const char *main_begin = R"(
void main() {
 gl_Position = pv * vec4(aPos.x, aPos.y, aPos.z, 1.0);
)";
  // 2897772.1
  const char *main_data_to_rgb = R"(
 float rr = 2800;
 float rg = 2000;
 float rb = 1200;
 float ti = 2;
 float wave_length = 660 - 260 * (data - min) / (max - min);
 if (min == max) {
  color = vec4(1, 1, 1, .6f);
 } else {
  color = vec4(
   0.41 * pow(2, -pow(abs(wave_length - 595), ti) / rr),
   0.82 * pow(2, -pow(abs(wave_length - 530), ti) / rg),
   0.40 * pow(2, -pow(abs(wave_length - 460), ti) / rb),
   .6f
  );
 }
)";
  const char *main_position_as_color = R"(
 color = vec4(aPos.x / X, aPos.y / Y, aPos.z / Z, .8f);
)";
  const char *main_highlight = R"(
 if (highlighted == data) {
  color = vec4(.8f, .1f, .0f, .6f);
 } else {
  color = vec4(.6f, .6f, .6f, .1f);
 }
)";
  const char *main_default_color = R"(
  color = vec4(.6f, .6f, .6f, .1f);
)";
  const char *main_end = R"(
}
)";

} vertexShaderSource;

const char *get_color(data_category category) {
  switch (category) {
  case data_category::temperature:
    return vertexShaderSource.main_data_to_rgb;
  case data_category::other:
  default:
    return vertexShaderSource.main_default_color;
  }
}

constexpr static inline auto &fragmentShaderSource = R"(
#version 330 core

in vec4 color;

void main() {
 gl_FragColor = color;
 // gl_FragColor.a = gl_FragColor.a / (1 + gl_FragCoord.z / gl_FragCoord.w);
}
)";

static inline GLenum defaultDrawMode() {
  return wireframe ? GL_LINES : GL_TRIANGLES;
}

static inline float patch_far(const vector<patch_info> &patches) {
  constexpr auto u32_max = numeric_limits<u32>::max();
  u32 I1 = u32_max, I2 = 0, J1 = u32_max, J2 = 0, K1 = u32_max, K2 = 0;
  for (const auto &patch : patches) {
    I1 = min(patch.I1, I1);
    I2 = max(patch.I2, I2);
    J1 = min(patch.J1, J1);
    J2 = max(patch.J2, J2);
    K1 = min(patch.K1, K1);
    K2 = max(patch.K2, K2);
  }
  float far = 100.;
  if (patches.size() > 0) {
    u32 I, J, K;
    I = I2 - I1;
    J = J2 - J1;
    K = K2 - K1;
    far = max(far, float(sqrt(I * I + J * J + K * K) * 5));
  }
  return far;
}

constexpr static inline float defaultNear = .01f;
constexpr static inline float defaultFar = 1000.f;

static inline int visualize_patch(const vector<patch_info> &patches) {
  index_max = patches.size();

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint, 2>(
        [&patches]() { return from_patches<true>(patches, wireframe); },
        nullptr,
        [&patches]() -> tuple<float, float> {
          return {defaultNear, patch_far(patches)};
        },
        defaultDrawMode, nullptr, nullptr,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.data,
         vertexShaderSource.highlight, vertexShaderSource.main_begin,
         vertexShaderSource.main_highlight, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<GLuint[3], GLuint>{});
  }
}
static inline int visualize_frames(const vector<patch_info> &patches,
                                   const vector<frame> &frames,
                                   data_category category) {
  assert(!frames.empty());
  index_max = frames.size();
  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint, 1>(
        [&patches, &frames ]() -> auto{
          glPointSize(10);

          return from_data(patches, frames.front());
        },
        [&frames]() -> vector<float> {
          auto i = current % frames.size();
          auto &frame = frames[i];
          return from_frame(frame, 0);
        },
        [&patches]() constexpr->tuple<float, float> {
          return {defaultNear, patch_far(patches)};
        },
        []() { return GL_POINTS; },
        [&frames]() { return frame_minmax(frames.back()); }, nullptr,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.data,
         vertexShaderSource.main_begin, vertexShaderSource.main_data_to_rgb,
         vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<GLuint[3], float>{});
  }
}

static inline int visualize_data(const vector<float> &data, u32 m, u32 n) {
  assert(data.size() == m * n);
  const auto Z = *max_element(data.cbegin(), data.cend());
  const auto Y = (float)n;
  const auto X = (float)m;
  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint, 1>(
        [&data, n]() { return from_matrix_data(data, n, wireframe); }, nullptr,
        [&data, X, Y, Z]() -> tuple<float, float> {
          return {defaultNear, 5 * max(initializer_list<float>{X, Y, Z})};
        },
        defaultDrawMode, nullptr, [X, Y, Z]() { return make_tuple(X, Y, Z); },
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.domain,
         vertexShaderSource.main_begin,
         vertexShaderSource.main_position_as_color,
         vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_3d_elements(const vector<float> &nodes,
                                        const vector<u32> &vertex_count,
                                        const vector<u32> &elements) {
  u32 sum = accumulate(vertex_count.cbegin(), vertex_count.cend(), 0);
  assert(elements.size() == sum);
  float _max = *max_element(nodes.cbegin(), nodes.cend());
  float ratio = _max == 0 ? 1 : 100 / _max;

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint, 1>(
        [&nodes, &vertex_count, &elements, ratio]() {
          return from_elements(nodes, vertex_count, elements, wireframe, ratio);
        },
        nullptr,
        [&nodes]() constexpr->tuple<float, float> {
          return {defaultNear, defaultFar};
        },
        defaultDrawMode, nullptr, nullptr,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_default_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_patches_and_elements(
    const vector<float> &nodes, const vector<u32> &vertex_count,
    const vector<u32> &elements, const vector<patch_info> &patches) {
  u32 sum = accumulate(vertex_count.cbegin(), vertex_count.cend(), 0);
  assert(elements.size() == sum);
  const float _max1 = patch_far(patches);
  const float _max2 =
      max(initializer_list<float>{*max_element(nodes.cbegin(), nodes.cend())});
  const float _max = max(_max1, _max2);
  const float ratio = _max == 0 ? 1 : 100 / _max;

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint, 1>(
        [&nodes, &vertex_count, &elements, &patches, sum, ratio]() {
          return from_patches_and_elements(nodes, vertex_count, elements,
                                           patches, wireframe, sum, ratio);
        },
        nullptr,
        [&nodes]() constexpr->tuple<float, float> {
          return {defaultNear, defaultFar};
        },
        defaultDrawMode, nullptr, nullptr,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_default_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_nodes(const vector<float> &nodes,
                                  const vector<u32> &indices) {
  float _max = *max_element(nodes.cbegin(), nodes.cend());
  float ratio = _max == 0 ? 1 : 100 / _max;

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint, 1>(
        [&nodes, &indices, ratio ]() -> auto{
          glPointSize(20);
          auto _nodes = nodes;
          for (auto &n : _nodes)
            n *= ratio;
          return make_tuple(make_tuple(_nodes), indices);
        },
        nullptr,
        [&nodes]() constexpr->tuple<float, float> {
          return {defaultNear, defaultFar};
        },
        []() { return GL_POINTS; }, nullptr, nullptr,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_default_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_primitives(const vector<float> &nodes,
                                       const vector<u32> &primitive) {
  float _max = *max_element(nodes.cbegin(), nodes.cend());
  float ratio = _max == 0 ? 1 : 100 / _max;

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint, 1>(
        [&nodes, &primitive, ratio ]() -> auto{
          auto _nodes = nodes;
          for (auto &n : _nodes)
            n *= ratio;
          return make_tuple(make_tuple(_nodes), primitive);
        },
        nullptr,
        [&nodes]() constexpr->tuple<float, float> {
          return {defaultNear, defaultFar};
        },
        defaultDrawMode, nullptr, nullptr,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_default_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}
[[deprecated]] static inline int visualize_primitives_on_boundary(
    const vector<patch_info> &patches, const vector<float> &nodes,
    const vector<u32> &vertex_count, const vector<u32> &primitive) {
  float _max = *max_element(nodes.cbegin(), nodes.cend());
  float ratio = _max == 0 ? 1 : 100 / _max;

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint, 1>(
        [&nodes, &patches, &vertex_count, &primitive, ratio ]() -> auto{
          auto indices =
              primitive_on_boundary(patches, nodes, primitive, wireframe);
          auto _nodes = nodes;
          for (auto &n : _nodes)
            n *= ratio;
          return make_tuple(make_tuple(_nodes), indices);
        },
        nullptr,
        [&nodes]() constexpr->tuple<float, float> {
          return {defaultNear, defaultFar};
        },
        defaultDrawMode, nullptr, nullptr,
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.main_begin,
         vertexShaderSource.main_default_color, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<float[3]>{});
  }
}

static inline int visualize_polygons(const vector<float> &nodes,
                                     const vector<u32> &polygon_sizes,
                                     const vector<u32> &polygon_indices) {
  return visualize_primitives(
      nodes, from_polygons(polygon_sizes, polygon_indices, wireframe));
}

#endif // GRAPHICS_ENABLED
