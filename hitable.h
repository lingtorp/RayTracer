#ifndef RAYTRACER_RAY_H
#define RAYTRACER_RAY_H

#include "vector.h"
#include <utility>

class Material;

struct Hit {
  float t = 0.0f;
  Vec3f p; // Position
  float u = 0.0f;
  float v = 0.0f;
  Vec3f normal;
  Material* mat = nullptr;
};

class Hitable {
public:
  virtual bool hit(const Ray& r, float t_min, double t_max, Hit& hit) const = 0;
};

#endif // RAYTRACER_RAY_H
