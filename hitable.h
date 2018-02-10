#ifndef RAYTRACER_RAY_H
#define RAYTRACER_RAY_H

#include "vector.h"

struct Hit {
  double t = 0.0;
  Vec3<> p;
  Vec3<> normal;
};

class Hitable {
public:
  virtual bool hit(const Ray& r, double t_min, double t_max, Hit& hit) const = 0;
};

#endif // RAYTRACER_RAY_H
