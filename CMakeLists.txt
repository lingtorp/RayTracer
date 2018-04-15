cmake_minimum_required(VERSION 3.9)
project(RayTracer)

set(CMAKE_CXX_STANDARD 17)

add_executable(RayTracer main.cpp vector.h hitable.h sphere.h world.h camera.h material.h random.h texture.h)

find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})
target_link_libraries(RayTracer ${SDL2_LIBRARIES})