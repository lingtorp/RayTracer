cmake_minimum_required(VERSION 3.9)
project(RayTracer)

set(CMAKE_CXX_STANDARD 17)

# Add pthread when compiling for Linux
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") # pthread for Linux
endif(UNIX)

add_executable(RayTracer main.cpp vector.h hitable.h sphere.h world.h camera.h material.h random.h texture.h)

find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})
target_link_libraries(RayTracer ${SDL2_LIBRARIES})