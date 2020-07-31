#ifndef RAYTRACER_MATERIAL_H
#define RAYTRACER_MATERIAL_H

#include "random.h"
#include "hitable.h"
#include "texture.h"

class Material {
public:
  /// Reflected ray
  virtual bool scatter(const Ray& r, const Hit& hit, Vec3f& attenuation, Ray& scattered) const = 0;
  /// Emitted color
  virtual Vec3f emitted(float u, float v, const Vec3f& p) const {
    return Vec3f{0.0};
  }
};

class Lambertian: public Material {
public:
  Texture* albedo;
  explicit Lambertian(Texture* texture): albedo{texture} {};
  
  bool scatter(const Ray& r, const Hit& hit, Vec3f& attenuation, Ray& scattered) const override {
    Vec3f target = hit.p + hit.normal + random_in_unit_sphere();
    scattered = Ray{hit.p, target - hit.p};
    attenuation = albedo->value(0.0, 0.0, hit.p);
    return true;
  }
};

class Metal: public Material {
public:
  Vec3f albedo;
  float fuzz = 0.0; // Fuzziness of the scattered rays
  Metal(float r, float g, float b): albedo{r, g, b} {};
  Metal(float r, float g, float b, float fuzz): albedo{r, g, b}, fuzz(fuzz) {};
  explicit Metal(const Vec3f& albedo): albedo{albedo} {};
  explicit Metal(const Vec3f& albedo, float fuzz): albedo{albedo}, fuzz(fuzz) {};
  
  bool scatter(const Ray& r, const Hit& hit, Vec3f& attenuation, Ray& scattered) const override {
    Vec3f reflected = reflect(r.direction().normalized(), hit.normal);
    scattered = Ray{hit.p, reflected + fuzz * random_in_unit_sphere()};
    attenuation = albedo;
    return (dot(scattered.direction(), hit.normal) > 0);
  }
};

class Dielectric: public Material {
public:
  float ref_idx; // Refraction index (see Snell's law)
  explicit Dielectric(float ref_idx): ref_idx{ref_idx} {};
  
  /// Approximiation by Christophe Schlick
  float schlick(float cosine, float ref_idx) const {
    float r0 = (1.0f - ref_idx) / (1.0f - ref_idx);
    r0 = r0 * r0;
    return r0 + (1.0f - r0) / std::pow((1.0f - cosine), 5.0f);
  }
  
  bool refract(const Vec3f v, const Vec3f& n, float ni_over_nt, Vec3f& refracted) const {
    Vec3f uv = v.normalized();
    float dt = dot(uv, n);
    float discriminant = 1.0f - ni_over_nt*ni_over_nt*(1.0f - dt*dt);
    if (discriminant > 0) {
      refracted = ni_over_nt*(uv - n*dt) - n*std::sqrt(discriminant);
      return true;
    } else {
      return false;
    }
  }
  
  // FIXME: When there is a reflection ray the functions returns false so there are no reflections
  bool scatter(const Ray& r, const Hit& hit, Vec3f& attenuation, Ray& scattered) const override {
    Vec3f out_n;
    Vec3f reflected = reflect(r.direction(), hit.normal);
    float ni_over_nt;
    attenuation = Vec3f{1.0, 1.0, 1.0};
    float cosine;
    if (dot(r.direction(), hit.normal) > 0) {
      out_n = -hit.normal;
      ni_over_nt = ref_idx;
      cosine = ref_idx * dot(r.direction(), hit.normal) / r.direction().length();
    } else {
      out_n = hit.normal;
      ni_over_nt = 1.0 / ref_idx;
      cosine = -dot(r.direction(), hit.normal) / r.direction().length();
    }
    Vec3f refracted;
    float reflect_prob;
    if (refract(r.direction(), out_n, ni_over_nt, refracted)) {
      reflect_prob = schlick(cosine, ref_idx);
    } else {
      reflect_prob = 1.0;
    }
    if (rand_0_1() < reflect_prob) {
      scattered = Ray{hit.p, reflected};
    } else {
      scattered = Ray{hit.p, refracted};
    }
    return true;
  }
};

class Emission : public Material {
private:
  Texture* texture;
  
public:
  explicit Emission(Texture* t): texture(t) {}
  
  bool scatter(const Ray &r, const Hit &hit, Vec3f &attenuation, Ray &scattered) const override {
    return false;
  }
  
  Vec3f emitted(float u, float v, const Vec3f &p) const override {
    return texture->value(u, v, p);
  }
};

#endif // RAYTRACER_MATERIAL_H
