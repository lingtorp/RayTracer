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

class Dielectric: public Material {
public:
  double ref_idx; // Refraction index (see Snell's law)
  explicit Dielectric(double ref_idx): ref_idx{ref_idx} {};
  
  /// Approximiation by Christophe Schlick
  double schlick(double cosine, double ref_idx) const {
    double r0 = (1.0 - ref_idx) / (1.0 - ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0 - r0) / std::pow((1.0 - cosine), 5.0);
  }
  
  bool refract(const Vec3<> v, const Vec3<>& n, double ni_over_nt, Vec3<>& refracted) const {
    Vec3<> uv = v.normalized();
    double dt = dot(uv, n);
    double discriminant = 1.0 - ni_over_nt*ni_over_nt*(1 - dt*dt);
    if (discriminant > 0) {
      refracted = ni_over_nt*(uv - n*dt) - n*std::sqrt(discriminant);
      return true;
    } else {
      return false;
    }
  }
  
  // FIXME: When there is a reflection ray the functions returns false so there are no reflections
  bool scatter(const Ray& r, const Hit& hit, Vec3<>& attenuation, Ray& scattered) const override {
    Vec3<> out_n;
    Vec3<> reflected = reflect(r.direction(), hit.normal);
    double ni_over_nt;
    attenuation = Vec3<>{1.0, 1.0, 1.0};
    double cosine;
    if (dot(r.direction(), hit.normal) > 0) {
      out_n = -hit.normal;
      ni_over_nt = ref_idx;
      cosine = ref_idx * dot(r.direction(), hit.normal) / r.direction().length();
    } else {
      out_n = hit.normal;
      ni_over_nt = 1.0 / ref_idx;
      cosine = -dot(r.direction(), hit.normal) / r.direction().length();
    }
    Vec3<> refracted;
    double reflect_prob;
    if (refract(r.direction(), out_n, ni_over_nt, refracted)) {
      reflect_prob = schlick(cosine, ref_idx);
    } else {
      reflect_prob = 1.0;
    }
    if (drand48() < reflect_prob) {
      scattered = Ray{hit.p, reflected};
    } else {
      scattered = Ray{hit.p, refracted};
    }
    return true;
  }
};

#endif // RAYTRACER_MATERIAL_H
