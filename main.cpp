#include <iostream>
#include <fstream>
#include "vector.h"

Vec3<> color(const Ray& ray) {
  Vec3<> dir = ray.direction().normalized();
  double t = 0.5 * (dir.y + 1.0);
  return (1.0 - t) * Vec3<>{1.0, 1.0, 1.0} + t * Vec3<>{0.5, 0.7, 1.0};
}

int main() {
  std::ofstream image;
  image.open("gradient.ppm");
  image << "P3" << '\n';
  size_t nx = 200;
  size_t ny = 100;
  image << nx << " " << ny << '\n' << 255 << '\n';
  Vec3<> low_left_corner{-2.0, -1.0, -1.0};
  Vec3<> horizontal{4.0, 0.0, 0.0};
  Vec3<> vertical{0.0, 2.0, 0.0};
  Vec3<> origin{};
  for (int j = ny - 1; j >= 0; j--) {
    for (int i = 0; i < nx; i++) {
      double u = double(i) / double(nx);
      double v = double(j) / double(ny);
      Ray r{origin, low_left_corner + u * horizontal + v * vertical};
      Vec3<> col = color(r);
      int ir = int(col.x * 255);
      int ig = int(col.y * 255);
      int ib = int(col.z * 255);
      image << ir << " " << ig << " " << ib << '\n';
    }
  }
  return 0;
}