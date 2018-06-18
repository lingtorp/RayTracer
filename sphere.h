#ifndef RAYTRACER_SPHERE_H
#define RAYTRACER_SPHERE_H

#include "hitable.h"
#include "material.h"

class Sphere: public Hitable {
public:
  Vec3<> center;
  double radius = 1.0;
  Material* material = nullptr;
  
  Sphere() = default;
  Sphere(Vec3<> center, float radius): center{center}, radius(radius) {};
  Sphere(Vec3<> center, float radius, Material* material): center{center}, radius(radius), material(material) {};

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
        hit.mat = material;
        return true;
      }
      tmp = (-b + std::sqrt(discriminant)) / (2.0*a);
      if (t_min < tmp && tmp < t_max) {
        hit.t = tmp;
        hit.p = r(hit.t);
        hit.normal = (hit.p - center) / radius;
        hit.mat = material;
        return true;
      }
    }
    return false; // Missed sphere
  }
  
  bool bounding_box(double t0, double t1, AABB& box) const override {
    box = AABB{center - Vec3<>{radius}, center + Vec3<>{radius}};
    return true;
  }
};

class Rectxy : public Hitable {
public:
  Material* mat = nullptr;
  double x0, x1, y0, y1, k;
  Rectxy(double x0, double x1, double y0, double y1, double k, Material* mat):
          x0(x0), x1(x1), y0(y0), y1(y1), k(k), mat(mat) {}
  
  bool hit(const Ray &r, double t_min, double t_max, Hit &hit) const override {
    float t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max) { return false; }
    const float x = r.origin().x + r.direction().x * t;
    const float y = r.origin().y + r.direction().y * t;
    if (x < x0 || x > x1 || y < y0 || y > y1) { return false; }
    
    hit.t = t;
    hit.mat = mat;
    hit.p = r(t);
    hit.normal = Vec3<>{0, 0, 1};
    return false;
  }
  
  bool bounding_box(double t0, double t1, AABB &box) const override {
    box = AABB{Vec3<>{x0, y0, k - 0.0001}, Vec3<>{x1, y1, k + 0.0001}};
    return true;
  }
};

#endif // RAYTRACER_SPHERE_H
