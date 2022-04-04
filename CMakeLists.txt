cmake_minimum_required(VERSION 3.0)

set(CMAKE_CXX_STANDARD 17)

project(main LANGUAGES CXX)

if(GRAPHICS_ENABLED)
    add_definitions(-DGRAPHICS_ENABLED=1)
else()
    add_definitions(-DGRAPHICS_ENABLED=0)
endif()

add_executable(reader main.cpp)

find_package(glad CONFIG REQUIRED)
target_link_libraries(reader PRIVATE glad::glad)

find_package(glfw3 CONFIG REQUIRED)
target_link_libraries(reader PRIVATE glfw)
