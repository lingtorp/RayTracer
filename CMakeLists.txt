cmake_minimum_required(VERSION 3.9)
project(RayTracer)

set(CMAKE_CXX_STANDARD 17)

# Add pthread when compiling for Linux
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
if (UNIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread") # pthread for Linux
endif(UNIX)

set(PROFILING_ENABLED TRUE)
if (UNIX AND PROFILING_ENABLED)
 #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pg") # Profiling with gprof
elseif(APPLE AND PROFILING_ENABLED)

endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-rtti")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-exceptions")

add_executable(RayTracer main.cpp vector.h hitable.h sphere.h world.h camera.h material.h random.h texture.h include/stb_image.h include/noise.h)

find_package(SDL2 REQUIRED)
include_directories(${SDL2_INCLUDE_DIRS})
target_link_libraries(RayTracer ${SDL2_LIBRARIES})