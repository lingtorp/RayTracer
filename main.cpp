#include <iostream>
#include <fstream>
#include <random>
#include "vector.h"
#include "world.h"
#include "sphere.h"
#include "camera.h"

/// dot((p(t) - c, p(t) - c)) = R*R sphere equation in vector form
Vec3<> color(const Ray& r, World& world) {
  Hit hit;
  if (world.hit(r, 0.0, MAXFLOAT, hit)) {
    return 0.5*(hit.normal + 1.0);
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
  
  std::ofstream image;
  image.open("gradient.ppm");
  image << "P3" << '\n';
  size_t nx = 200;
  size_t ny = 100;
  size_t ns = 100; // Number of samples per px
  image << nx << " " << ny << '\n' << 255 << '\n';
  Camera cam;
  World world;
  world.hitables.push_back(new Sphere{Vec3<>{0, 0, -1}, 0.5});
  world.hitables.push_back(new Sphere{Vec3<>{0, -100.5, -1}, 100});
  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {
      Vec3<> t_color{};
      for (int s = 0; s < ns; s++) {
        double u = double(i + rand()) / double(nx);
        double v = double(j + rand()) / double(ny);
        Ray r = cam.get_ray(u, v);
  
        Vec3<> p = r(2.0);
        t_color += color(r, world);
      }
      t_color /= double(ns);
      
      int ir = int(t_color.x * 255);
      int ig = int(t_color.y * 255);
      int ib = int(t_color.z * 255);
      image << ir << " " << ig << " " << ib << '\n';
    }
  }
  return 0;
}