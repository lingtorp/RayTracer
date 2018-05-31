#ifndef RAYTRACER_TEXTURE_H
#define RAYTRACER_TEXTURE_H

#include "vector.h"
#include "Procedural-Noise/noise.hpp"

class Texture {
public:
  virtual Vec3<> value(double u, double v, const Vec3<>& p) const = 0;
};

class ConstantTexture : public Texture {
public:
  Vec3<> rgb;
  ConstantTexture(double r, double g, double b): rgb{r, g, b} {}
  explicit ConstantTexture(Vec3<> rgb_color): rgb{rgb_color} {}
  
  Vec3<> value(double u, double v, const Vec3<> &p) const override {
    return rgb;
  }
};

class CheckerTexture : public Texture {
public:
  Texture* tex0;
  Texture* tex1;
  CheckerTexture(Texture* tex0, Texture* tex1): tex0{tex0}, tex1{tex1} {}
  
  Vec3<> value(double u, double v, const Vec3<> &p) const override {
    double sines = sin(10*p.x)*sin(10.0*p.y)*sin(10.0*p.z);
    if (sines < 0) {
      return tex0->value(u, v, p);
    } else {
      return tex1->value(u, v, p);
    }
  }
};

class PerlinTexture : public Texture {
  Perlin::Original noise_gen;
public:
  PerlinTexture(): noise_gen{1337} {};
  Vec3<> value(double u, double v, const Vec3<> &p) const override {
    const uint32_t zoom_factor = 128;
    Vec3<> color = Vec3<>{noise_gen.fbm(1000 * p, zoom_factor)};
    if (color <= 0.05) {
      return Vec3<>{0.05, 0.05, 0.05};
    }
    return color;
  }
};

#endif // RAYTRACER_TEXTURE_H
