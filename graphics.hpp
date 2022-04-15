#pragma once
#include "types.hpp"
#include <algorithm>
#include <functional>
#include <iostream>
#include <limits>
#include <vector>
using std::cerr;
using std::cin;
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
using std::max;
using std::min;
using std::numeric_limits;
using std::remove_all_extents_t;
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

#if GRAPHICS_ENABLED

template <typename> struct gl_type_enum;
template <>
struct gl_type_enum<GLuint> : integral_constant<GLenum, GL_UNSIGNED_INT> {};
template <>
struct gl_type_enum<GLfloat> : integral_constant<GLenum, GL_FLOAT> {};

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

static inline void onKey(GLFWwindow *window, int key, int scancode, int action,
                         int mods) {
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

  const bool continuous = mods & GLFW_MOD_CAPS_LOCK;

  if (key == GLFW_KEY_LEFT)
    if (action == GLFW_PRESS || (continuous && action == GLFW_REPEAT)) {
      if (current > 0)
        --current;
      else if (patch_loop)
        current = patches_count;
    }

  if (key == GLFW_KEY_RIGHT)
    if (action == GLFW_PRESS || (continuous && action == GLFW_REPEAT)) {
      if (current < patches_count) {
        ++current;
      } else if (patch_loop)
        current = 0;
    }
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

static GLuint windowWidth = 800;
static GLuint windowHeight = 600;
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
  glShaderSource(vertexShader, vertexShaderSource.size(),
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
  glShaderSource(fragmentShader, fragmentShaderSource.size(),
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

/// @retval Whether to start visualization.
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
    case 'h': {
      cout << R"(
Press TAB to enable or disable cursor.
Press CAPS to enable or disable continuous key input.
Press LEFT or RIGHT to change current highlighted patch.
Press ESC to terminate.
Move mouse to rotate the camera when cursor is disabled.
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

template <typename Index, typename GetData, typename GetFar, typename... Vertex>
static inline int visualize(GetData &&getData, GetFar &&getFar,
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
  glBufferData(GL_ELEMENT_ARRAY_BUFFER,
               indices.size() * sizeof(decltype(indices)::value_type),
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

  float far = getFar();

  glm::mat4 projection = glm::perspective(
      glm::radians(60.0f), float(windowWidth) / windowHeight, 0.01f, far);

  DETECT_ERROR;

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
    glDrawElements(wireframe ? GL_LINES : GL_TRIANGLES, (GLuint)indices.size(),
                   gl_type_enum_v<Index>, 0);
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
layout(location = 1) in uint index;
)";
  const char *highlight = R"(
uniform uint highlighted;
)";
  const char *main_begin = R"(
void main() {
 gl_Position = pv * vec4(aPos.x, aPos.y, aPos.z, 1.0);
)";
  // DEBUG:
  const char *main_highlight = R"(
 if (highlighted <= index) {
  color = vec4(.8f, .1f, .0f, .8f);
 } else {
  color = vec4(.8f, .8f, .8f, .2f);
 }
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

static inline int visualize_patch(const vector<patch_info> &patches) {

  patches_count = patches.size();
  if (current > patches.size()) {
    current = 0;
  }

  while (true) {
    if (!visualization_settings())
      return 0;

    visualize<GLuint>(
        [&patches]() {
          tuple<tuple<vector<GLuint>, vector<GLuint>>, vector<GLuint>> res;
          auto &[vertices, indices] = res;
          auto &[position, patchIndices] = vertices;

          position.reserve(patches.size() * 12);
          patchIndices.reserve(patches.size() * 4);
          indices.reserve(wireframe ? indices.size() * 8 : indices.size() * 6);

          size_t i_patch = 0;
          for (const auto &patch : patches) {
            const GLuint N = (GLuint)position.size() / 3;
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

            for (size_t i = 0; i < 4; ++i)
              patchIndices.push_back(i_patch);

            ++i_patch;
          }

          return res;
        },
        [&patches]() {
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

          return far;
        },
        {vertexShaderSource.head_pos_pv_color, vertexShaderSource.index,
         vertexShaderSource.highlight, vertexShaderSource.main_begin,
         vertexShaderSource.main_highlight, vertexShaderSource.main_end},
        {fragmentShaderSource}, type_list<GLuint[3], GLuint>{});
  }
}

#endif // GRAPHICS_ENABLED
