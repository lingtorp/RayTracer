#ifndef RAYTRACER_TEXTURE_H
#define RAYTRACER_TEXTURE_H

#include "vector.h"

class Texture {
public:
  virtual Vec3<> value(double u, double v, const Vec3<>& p) const = 0;
};

class ConstantTexture : public Texture {
public:
  Vec3<> rgb;
  ConstantTexture(double r, double g, double b): rgb{r, g, b} {}
  ConstantTexture(Vec3<> rgb_color): rgb{rgb_color} {}
  
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

#endif // RAYTRACER_TEXTURE_H
