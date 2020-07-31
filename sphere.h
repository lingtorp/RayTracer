#ifndef RAYTRACER_SPHERE_H
#define RAYTRACER_SPHERE_H

#include "hitable.h"
#include "material.h"

class Sphere: public Hitable {
public:
  Vec3f center;
  float radius = 1.0;
  Material* material = nullptr;
  
  Sphere() = default;
  Sphere(Vec3f center, float radius): center{center}, radius(radius) {};
  Sphere(Vec3f center, float radius, Material* material): center{center}, radius(radius), material(material) {};

  bool hit(const Ray& r, float t_min, double t_max, Hit& hit) const override {
    Vec3f oc = r.origin() - center; // Origin to center
    float a = dot(r.direction(), r.direction());
    float b = 2 * dot(r.direction(), oc);
    float c = dot(oc, oc) - radius*radius;
    float discriminant = b*b - 4*a*c;
    if (discriminant > 0) { // Hit sphere
      float tmp = (-b - std::sqrt(discriminant)) / (2.0*a);
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
};

class Rectxy : public Hitable {
public:
  Material* mat = nullptr;
  float x0, x1, y0, y1, k;
  Rectxy(float x0, double x1, double y0, double y1, double k, Material* mat):
          x0(x0), x1(x1), y0(y0), y1(y1), k(k), mat(mat) {}
  
  bool hit(const Ray &r, float t_min, double t_max, Hit &hit) const override {
    const float t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max) { return false; }
    const float x = r.origin().x + t * r.direction().x;
    const float y = r.origin().y + t * r.direction().y;
    if (x < x0 || x > x1 || y < y0 || y > y1) { return false; }
    
    hit.t = t;
    hit.mat = mat;
    hit.p = r(t);
    hit.u = (x - x0) / (x1 - x0);
    hit.v = (y - y0) / (y1 - y0);
    hit.normal = Vec3f{0, 0, 1};
    return true;
  }
};

class Rectxz : public Hitable {
public:
  Material* mat = nullptr;
  float x0, x1, z0, z1, k;
  Rectxz(float x0, double x1, double z0, double z1, double k, Material* mat):
          x0(x0), x1(x1), z0(z0), z1(z1), k(k), mat(mat) {}
  
  bool hit(const Ray &r, float t_min, double t_max, Hit &hit) const override {
    const float t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max) { return false; }
    const float x = r.origin().x + t * r.direction().x;
    const float z = r.origin().z + t * r.direction().z;
    if (x < x0 || x > x1 || z < z0 || z > z1) { return false; }
    
    hit.t = t;
    hit.mat = mat;
    hit.p = r(t);
    hit.u = (x - x0) / (x1 - x0);
    hit.v = (z - z0) / (z1 - z0);
    hit.normal = Vec3f{0, 1, 0};
    return true;
  }
};

class Rectyz : public Hitable {
public:
  Material* mat = nullptr;
  float y0, y1, z0, z1, k;
  Rectyz(float y0, double y1, double z0, double z1, double k, Material* mat):
          z0(z0), z1(z1), y0(y0), y1(y1), k(k), mat(mat) {}
  
  bool hit(const Ray &r, float t_min, double t_max, Hit &hit) const override {
    const float t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max) { return false; }
    const float z = r.origin().z + t * r.direction().z;
    const float y = r.origin().y + t * r.direction().y;
    if (z < z0 || z > z1 || y < y0 || y > y1) { return false; }
    
    hit.t = t;
    hit.mat = mat;
    hit.p = r(t);
    hit.u = (y - y0) / (y1 - y0);
    hit.v = (z - z0) / (z1 - z0);
    hit.normal = Vec3f{1.0f, 0.0f, 0.0f};
    return true;
  }
};

#endif // RAYTRACER_SPHERE_H
