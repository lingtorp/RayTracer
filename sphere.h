#ifndef RAYTRACER_SPHERE_H
#define RAYTRACER_SPHERE_H

#include "hitable.h"
#include "material.h"

void get_sphere_uv(const Vec3<>& p, double& u, double& v) {
  double phi = atan2(p.z, p.x);
  double theta = asin(p.y);
  u = 1 - (phi + M_PI) / (2*M_PI);
  v = (theta + M_PI/2) / M_PI;
}

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
        get_sphere_uv(hit.p, hit.u, hit.v);
        hit.normal = (hit.p - center) / radius;
        hit.mat = material;
        return true;
      }
      tmp = (-b + std::sqrt(discriminant)) / (2.0*a);
      if (t_min < tmp && tmp < t_max) {
        hit.t = tmp;
        hit.p = r(hit.t);
        get_sphere_uv(hit.p, hit.u, hit.v);
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
    const double t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max) { return false; }
    const double x = r.origin().x + t * r.direction().x;
    const double y = r.origin().y + t * r.direction().y;
    if (x < x0 || x > x1 || y < y0 || y > y1) { return false; }
    
    hit.t = t;
    hit.mat = mat;
    hit.p = r(t);
    hit.u = (x - x0) / (x1 - x0);
    hit.v = (y - y0) / (y1 - y0);
    hit.normal = Vec3<>{0, 0, 1};
    return true;
  }
  
  bool bounding_box(double t0, double t1, AABB &box) const override {
    box = AABB{Vec3<>{x0, y0, k - 0.0001}, Vec3<>{x1, y1, k + 0.0001}};
    return true;
  }
};

class Rectxz : public Hitable {
public:
  Material* mat = nullptr;
  double x0, x1, z0, z1, k;
  Rectxz(double x0, double x1, double z0, double z1, double k, Material* mat):
          x0(x0), x1(x1), z0(z0), z1(z1), k(k), mat(mat) {}
  
  bool hit(const Ray &r, double t_min, double t_max, Hit &hit) const override {
    const double t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max) { return false; }
    const double x = r.origin().x + t * r.direction().x;
    const double z = r.origin().z + t * r.direction().z;
    if (x < x0 || x > x1 || z < z0 || z > z1) { return false; }
    
    hit.t = t;
    hit.mat = mat;
    hit.p = r(t);
    hit.u = (x - x0) / (x1 - x0);
    hit.v = (z - z0) / (z1 - z0);
    hit.normal = Vec3<>{0, 0, 1};
    return true;
  }
  
  bool bounding_box(double t0, double t1, AABB &box) const override {
    box = AABB{Vec3<>{x0, z0, k - 0.0001}, Vec3<>{x1, z1, k + 0.0001}};
    return true;
  }
};

class Rectyz : public Hitable {
public:
  Material* mat = nullptr;
  double y0, y1, z0, z1, k;
  Rectyz(double y0, double y1, double z0, double z1, double k, Material* mat):
          z0(z0), z1(z1), y0(y0), y1(y1), k(k), mat(mat) {}
  
  bool hit(const Ray &r, double t_min, double t_max, Hit &hit) const override {
    const double t = (k - r.origin().z) / r.direction().z;
    if (t < t_min || t > t_max) { return false; }
    const double z = r.origin().z + t * r.direction().z;
    const double y = r.origin().y + t * r.direction().y;
    if (z < z0 || z > z1 || y < y0 || y > y1) { return false; }
    
    hit.t = t;
    hit.mat = mat;
    hit.p = r(t);
    hit.u = (y - y0) / (y1 - y0);
    hit.v = (z - z0) / (z1 - z0);
    hit.normal = Vec3<>{0, 0, 1};
    return true;
  }
  
  bool bounding_box(double t0, double t1, AABB &box) const override {
    box = AABB{Vec3<>{z0, y0, k - 0.0001}, Vec3<>{z1, y1, k + 0.0001}};
    return true;
  }
};

#endif // RAYTRACER_SPHERE_H
