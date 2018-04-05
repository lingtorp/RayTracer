#include <iostream>
#include <fstream>
#include <random>
#include <SDL2/SDL.h>
#include <SDL2/SDL_video.h>
#include "vector.h"
#include "world.h"
#include "sphere.h"
#include "camera.h"
#include "material.h"

/// dot((p(t) - c, p(t) - c)) = R*R sphere equation in vector form
Vec3<> color(const Ray& r, World& world, int depth) {
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

int main() {
  SDL_Init(SDL_INIT_EVERYTHING);

  size_t nx = 400;
  size_t ny = 200;
  size_t ns = 10; // Number of samples per px
  Camera cam{{-2, 2, 1}, {0, 0, -1}, {0, 1, 0}, 90, double(nx) / double(ny)};

  SDL_Window* window = SDL_CreateWindow("RayTracer", 0, 0, nx, ny, 0);
  SDL_Surface* scr = SDL_GetWindowSurface(window);
  uint32_t* pixels = (uint32_t*) scr->pixels;

  std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
  std::default_random_engine engine{static_cast<unsigned int>(tp.time_since_epoch().count())};
  std::uniform_real_distribution<double> distribution(0.0, 1.0); // 0.0 <= x < 1.0
  auto rand = std::bind(distribution, engine);

  World world;
  world.hitables.push_back(new Sphere{Vec3<>{0, 0, -1}, 0.5, new Lambertian{0.8, 0.3, 0.3}});
  world.hitables.push_back(new Sphere{Vec3<>{0, -100.5, -1}, 100, new Lambertian{0.8, 0.8, 0.0}});
  world.hitables.push_back(new Sphere{Vec3<>{1, 0, -1}, 0.5, new Metal{0.8, 0.8, 0.0, 0.3}});
  world.hitables.push_back(new Sphere{Vec3<>{-1, 0, -1}, 0.5, new Dielectric{1.5}});

  time_t start = std::time(nullptr);
  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {
      Vec3<> t_color{};
      for (int s = 0; s < ns; s++) {
        double u = double(i + rand()) / double(nx);
        double v = double(j + rand()) / double(ny);
        Ray r = cam.get_ray(u, v);

        Vec3<> p = r(2.0);
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
      pixels[(j * nx) + i] = pixel;
    }
  }
  
  time_t end = std::time(nullptr);
  double diff = std::difftime(end, start);
  std::cout << "Finished in: " << diff << " seconds." << std::endl;

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