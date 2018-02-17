#ifndef RAYTRACER_RAY_H
#define RAYTRACER_RAY_H

#include "vector.h"

class Material;

struct Hit {
  double t = 0.0;
  Vec3<> p{}; // Position
  Vec3<> normal{};
  Material* mat = nullptr;
};

class Hitable {
public:
  virtual bool hit(const Ray& r, double t_min, double t_max, Hit& hit) const = 0;
};

#endif // RAYTRACER_RAY_H
