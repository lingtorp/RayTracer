#ifndef RAYTRACER_MATERIAL_H
#define RAYTRACER_MATERIAL_H

#include <random>
#include "hitable.h"

class Material {
public:
  virtual bool scatter(const Ray& r, const Hit& hit, Vec3<>& attenuation, Ray& scattered) const = 0;
};

class Lambertian: public Material {
private:
  Vec3<> random_in_unit_sphere() const {
    std::chrono::system_clock::time_point tp = std::chrono::system_clock::now();
    std::default_random_engine engine{static_cast<unsigned int>(tp.time_since_epoch().count())};
    std::uniform_real_distribution<double> distribution(0.0, 1.0); // 0.0 <= x < 1.0
    auto rand = std::bind(distribution, engine);
    Vec3<> p;
    do {
      p = 2.0 * Vec3<>{rand(), rand(), rand()} - Vec3<>{1.0};
    } while (p.squared_length() >= 1.0);
    return p;
  }
  
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
  Metal(double r, double g, double b): albedo{r, g, b} {};
  explicit Metal(const Vec3<>& albedo): albedo{albedo} {};
  
  bool scatter(const Ray& r, const Hit& hit, Vec3<>& attenuation, Ray& scattered) const override {
    Vec3<> reflected = reflect(r.direction().normalized(), hit.normal);
    scattered = Ray{hit.p, reflected};
    attenuation = albedo;
    return (dot(scattered.direction(), hit.normal) > 0);
  }
};

#endif // RAYTRACER_MATERIAL_H
