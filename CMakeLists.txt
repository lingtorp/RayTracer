cmake_minimum_required(VERSION 3.9)
project(RayTracer)

set(CMAKE_CXX_STANDARD 17)

add_executable(RayTracer main.cpp vector.h hitable.h sphere.h world.h camera.h material.h)