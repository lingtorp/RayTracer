#ifndef RAYTRACER_TEXTURE_H
#define RAYTRACER_TEXTURE_H

#include "vector.h"
#include "include/noise.h"

class Texture {
public:
  virtual Vec3f value(float u, float v, const Vec3f& p) const = 0;
};

class ConstantTexture : public Texture {
public:
  Vec3f rgb;
  ConstantTexture(float r, float g, float b): rgb{r, g, b} {}
  explicit ConstantTexture(Vec3f rgb_color): rgb{rgb_color} {}
  
  Vec3f value(float u, float v, const Vec3f &p) const override {
    return rgb;
  }
};

class CheckerTexture : public Texture {
public:
  Texture* tex0 = nullptr;
  Texture* tex1 = nullptr;
  CheckerTexture(Texture* tex0, Texture* tex1): tex0{tex0}, tex1{tex1} {}
  
  Vec3f value(float u, float v, const Vec3f &p) const override {
    float sines = sin(10.0f*p.x) * sin(10.0f*p.y) * sin(10.0f*p.z);
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
  Vec3f value(float u, float v, const Vec3f &p) const override {
    const uint32_t zoom_factor = 128;
    Vec3f color = Vec3f{noise_gen.fbm(1000 * p, zoom_factor)};
    if (color <= 0.05) {
      return Vec3f{0.05, 0.05, 0.05};
    }
    return color;
  }
};

#define STB_IMAGE_IMPLEMENTATION
#include "include/stb_image.h"

class ImageTexture : public Texture {
public:
  int nx = 0;
  int ny = 0;
  int nn = 0;
  uint8_t* data = nullptr;
  
  explicit ImageTexture(const std::string& filename) {
    data = stbi_load(filename.c_str(), &nx, &ny, &nn, 0);
    if (data == nullptr) {
      std::cerr << "ImageTexture failed with error: " << stbi_failure_reason() << std::endl;
    }
  }
  
  ~ImageTexture() { stbi_image_free(data); }
  
  Vec3f value(float u, float v, const Vec3f &p) const override {
    if (data == nullptr) { return Vec3f{}; }
    int i = int(u*nx);
    int j = int((1-v)*ny-0.001f);
    if (i < 0) i = 0;
    if (j < 0) j = 0;
    if (i > nx - 1) i = nx - 1;
    if (j > ny - 1) j = ny - 1;
    float r = int(data[3*i + 3*nx*j])   / 255.0f;
    float g = int(data[3*i + 3*nx*j+1]) / 255.0f;
    float b = int(data[3*i + 3*nx*j+2]) / 255.0f;
    return Vec3f{r, g, b};
  }
};

#endif // RAYTRACER_TEXTURE_H
