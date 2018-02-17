#ifndef RAYTRACER_MATERIAL_H
#define RAYTRACER_MATERIAL_H

#include <random>
#include "hitable.h"

class Material {
public:
  virtual bool scatter(const Ray& r, const Hit& hit, Vec3<>& attenuation, Ray& scattered) const = 0;
};

class Lambertian: public Material {
public:
  Vec3<> albedo;
  Lambertian(double r, double g, double b): albedo{r, g, b} {};
  explicit Lambertian(const Vec3<>& albedo): albedo{albedo} {};
  
  bool scatter(const Ray& r, const Hit& hit, Vec3<>& attenuation, Ray& scattered) const override {
    Vec3<> target = hit.p + hit.normal + random_in_unit_sphere();
    scattered = Ray{hit.p, target - hit.p};
    attenuation = albedo;
    return true;
  }
};

class Metal: public Material {
public:
  Vec3<> albedo;
  double fuzz = 0.0; // Fuzziness of the scattered rays
  Metal(double r, double g, double b): albedo{r, g, b} {};
  Metal(double r, double g, double b, double fuzz): albedo{r, g, b}, fuzz(fuzz) {};
  explicit Metal(const Vec3<>& albedo): albedo{albedo} {};
  explicit Metal(const Vec3<>& albedo, double fuzz): albedo{albedo}, fuzz(fuzz) {};
  
  bool scatter(const Ray& r, const Hit& hit, Vec3<>& attenuation, Ray& scattered) const override {
    Vec3<> reflected = reflect(r.direction().normalized(), hit.normal);
    scattered = Ray{hit.p, reflected + fuzz * random_in_unit_sphere()};
    attenuation = albedo;
    return (dot(scattered.direction(), hit.normal) > 0);
  }
};

#endif // RAYTRACER_MATERIAL_H
