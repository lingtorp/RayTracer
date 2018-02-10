#include <iostream>
#include <fstream>
#include "vector.h"
#include "world.h"
#include "sphere.h"

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
  std::ofstream image;
  image.open("gradient.ppm");
  image << "P3" << '\n';
  size_t nx = 200;
  size_t ny = 100;
  image << nx << " " << ny << '\n' << 255 << '\n';
  Vec3<> low_left_corner{-2.0, -1.0, -1.0}; // Frame a bit in front of camera
  Vec3<> horizontal{4.0, 0.0, 0.0};
  Vec3<> vertical{0.0, 2.0, 0.0};
  Vec3<> origin{};
  World world;
  world.hitables.push_back(new Sphere{Vec3<>{0, 0, -1}, 0.5});
  world.hitables.push_back(new Sphere{Vec3<>{0, -100.5, -1}, 100});
  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {
      double u = double(i) / double(nx);
      double v = double(j) / double(ny);
      Ray r{origin, low_left_corner + u * horizontal + v * vertical};
      
      Vec3<> p = r(2.0);
      Vec3<> col = color(r, world);
      
      int ir = int(col.x * 255);
      int ig = int(col.y * 255);
      int ib = int(col.z * 255);
      image << ir << " " << ig << " " << ib << '\n';
    }
  }
  return 0;
}