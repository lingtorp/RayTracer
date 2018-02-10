#ifndef RAYTRACER_SPHERE_H
#define RAYTRACER_SPHERE_H

#include "hitable.h"

class Sphere: public Hitable {
public:
  Vec3<> center;
  double radius = 1.0;
  
  Sphere() = default;
  Sphere(Vec3<> center, float radius): center{center}, radius(radius) {};

  bool hit(const Ray& r, double t_min, double t_max, Hit& hit) const override {
    Vec3<> oc = r.origin() - center; // Origin to center
    double a = dot(r.direction(), r.direction());
    double b = 2 * dot(r.direction(), oc);
    double c = dot(oc, oc) - radius*radius;
    double discriminant = b*b - 4*a*c;
    if (discriminant > 0) { // Hit sphere
      double tmp = (-b - std::sqrt(discriminant)) / (2.0*a);
      if (t_min < tmp && tmp < t_max) {
        hit.t = tmp;
        hit.p = r(hit.t);
        hit.normal = (hit.p - center) / radius;
        return true;
      }
      tmp = (-b + std::sqrt(discriminant)) / (2.0*a);
      if (t_min < tmp && tmp < t_max) {
        hit.t = tmp;
        hit.p = r(hit.t);
        hit.normal = (hit.p - center) / radius;
        return true;
      }
      return true;
    }
    return false; // Missed sphere
  }
};

#endif // RAYTRACER_SPHERE_H
