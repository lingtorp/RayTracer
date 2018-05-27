#include <iostream>
#include <fstream>
#include <random>
#include "random.h"
#include <SDL2/SDL.h>
#include <SDL2/SDL_video.h>
#include <thread>
#include "vector.h"
#include "world.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"

/// dot((p(t) - c, p(t) - c)) = R*R sphere equation in vector form
Vec3<> color(const Ray& r, const World& world, int depth) {
  Hit hit;
  if (world.hit(r, 0.001, MAXFLOAT, hit)) {
    Ray scattered;
    Vec3<> attenuation;
    if (depth < 50 && hit.mat->scatter(r, hit, attenuation, scattered)) {
      return attenuation*color(scattered, world, depth + 1);
    } else {
      return Vec3<>{0.0};
    }
  } else {
    Vec3<> dir = r.direction().normalized();
    double t = 0.5*(dir.y + 1.0);
    return (1.0 - t) * Vec3<>{1.0, 1.0, 1.0} + t * Vec3<>{0.5, 0.7, 1.0};
  }
}

struct Region {
  const size_t nx, ny;
  size_t x0, x1;
  size_t y0, y1;
};

void thread_work(const World& world, const Camera& cam, uint32_t* pixels, int ns, Region reg) {
  for (size_t j = reg.y0; j <= reg.y1; j++) {
    for (size_t i = reg.x0; i <= reg.x1; i++) {
      Vec3<> t_color{};
      for (int s = 0; s < ns; s++) {
        double u = double(i + drand48()) / double(reg.nx);
        double v = double(j + drand48()) / double(reg.ny);
        Ray r = cam.get_ray(u, v);
        t_color += color(r, world, 0);
      }
      t_color /= double(ns);
      t_color = {std::sqrt(t_color.x), std::sqrt(t_color.y), std::sqrt(t_color.z)}; // Gamma-2 correction
      auto ir = uint32_t(t_color.x * 255);
      auto ig = uint32_t(t_color.y * 255);
      auto ib = uint32_t(t_color.z * 255);
      auto ia = uint32_t(1);
      uint32_t pixel = 0;
      pixel += (ia << (8 * 3));
      pixel += (ir << (8 * 2));
      pixel += (ig << (8 * 1));
      pixel += (ib << (8 * 0));
      pixels[((reg.ny - j) * reg.nx) + i] = pixel;
    }
  }
}

int main() {
  SDL_Init(SDL_INIT_EVERYTHING);

  // Setup random generator state for Linux
#ifdef __LINUX__
  init_random_generator();
#endif

  const size_t nx = 720;
  const size_t ny = 400;
  const size_t ns = 10; // Number of samples per px
  Vec3<> lookfrom = {3, 3, 2};
  Vec3<> lookat = {0, 0, -1};
  double dist_to_focus = (lookfrom - lookat).length();
  Camera cam{lookfrom, lookat, 20, double(nx) / double(ny), 0.1, dist_to_focus};

  SDL_Window* window = SDL_CreateWindow("RayTracer", 0, 0, nx, ny, 0);
  SDL_Surface* scr = SDL_GetWindowSurface(window);
  uint32_t* pixels = (uint32_t*) scr->pixels;

  World world;
  auto checker_texture = new CheckerTexture{new ConstantTexture{0.2, 0.3, 0.1}, new ConstantTexture{0.9, 0.9, 0.9}};
  world.hitables.push_back(new Sphere{Vec3<>{0, -100.5, -1}, 100, new Lambertian{checker_texture}});
  world.hitables.push_back(new Sphere{Vec3<>{0, 0, -1}, 0.5, new Lambertian{new ConstantTexture{0.8, 0.3, 0.3}}});
  world.hitables.push_back(new Sphere{Vec3<>{1, 0, -1}, 0.5, new Metal{0.8, 0.8, 0.0, 0.3}});
  world.hitables.push_back(new Sphere{Vec3<>{-1, 0, -1}, 0.5, new Dielectric{1.5}});
  
  const size_t num_threads = std::thread::hardware_concurrency() == 0 ? 4 : std::thread::hardware_concurrency();
  std::cout << "Starting " << num_threads << " number of threads." << std::endl;
  std::vector<std::thread> threads{};
  auto start = std::chrono::high_resolution_clock::now();
  for (size_t i = 0; i < num_threads; i++) {
    const size_t num_rows = ny / num_threads;
    Region reg{nx, ny, 0, nx, i * num_rows, (i + 1) * num_rows};
    threads.emplace_back(std::thread{thread_work, world, cam, pixels, ns, reg});
  }
  for (size_t i = 0; i < num_threads; i++) {
    threads[i].join();
  }
  auto end = std::chrono::high_resolution_clock::now();
  auto diff = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  std::cout << diff << " ms, " << (ns*nx*ny)/(diff/1000.0) << " rays/s" << std::endl;
  
  SDL_UpdateWindowSurface(window);
  
  bool quit = false;
  while (!quit) {
    SDL_Event event{};
    while (SDL_PollEvent(&event)) {
      switch (event.type) {
        case SDL_QUIT:
          quit = true;
          break;
      }
    }
  }
  
  return 0;
}