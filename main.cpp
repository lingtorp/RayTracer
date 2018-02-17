#include <iostream>
#include <fstream>
#include <random>
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
  std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
  std::default_random_engine engine{static_cast<unsigned int>(tp.time_since_epoch().count())};
  std::uniform_real_distribution<double> distribution(0.0, 1.0); // 0.0 <= x < 1.0
  auto rand = std::bind(distribution, engine);
  
  time_t start = std::time(nullptr);
  
  std::ofstream image;
  image.open("gradient.ppm");
  image << "P3" << '\n';
  size_t nx = 400;
  size_t ny = 200;
  size_t ns = 100; // Number of samples per px
  image << nx << " " << ny << '\n' << 255 << '\n';
  Camera cam{120, double(nx) / double(ny)};
  World world;
  world.hitables.push_back(new Sphere{Vec3<>{0, 0, -1}, 0.5, new Lambertian{0.8, 0.3, 0.3}});
  world.hitables.push_back(new Sphere{Vec3<>{0, -100.5, -1}, 100, new Lambertian{0.8, 0.8, 0.0}});
  world.hitables.push_back(new Sphere{Vec3<>{1, 0, -1}, 0.5, new Metal{0.8, 0.8, 0.0, 0.3}});
  world.hitables.push_back(new Sphere{Vec3<>{-1, 0, -1}, 0.5, new Dielectric{1.5}});
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
      int ir = int(t_color.x * 255);
      int ig = int(t_color.y * 255);
      int ib = int(t_color.z * 255);
      image << ir << " " << ig << " " << ib << '\n';
    }
  }
  
  time_t end = std::time(nullptr);
  double diff = std::difftime(end, start);
  std::cout << "Finished in: " << diff << " seconds." << std::endl;
  
  return 0;
}