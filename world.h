#ifndef RAYTRACER_WORLD_H
#define RAYTRACER_WORLD_H

#include <vector>
#include "hitable.h"
#include "sphere.h"

struct World {
  std::vector<Hitable*> hitables;
  
  World() = default;
  
  bool hit(const Ray& r, float t_min, float t_max, Hit& hit) const {
    bool hit_anything = false;
    float closest_so_far = t_max;
    for (const auto& hitable : hitables) {
      if (hitable->hit(r, t_min, closest_so_far, hit)) {
        closest_so_far = hit.t;
        hit_anything = true;
      }
    }
    return hit_anything;
  }
  
  /// Scenes
  void default_scene() {
    auto checker_texture = new CheckerTexture{new ConstantTexture{0.2, 0.3, 0.1}, new ConstantTexture{0.9, 0.9, 0.9}};
    hitables.push_back(new Sphere{Vec3f{0, -100.5, -1}, 100, new Lambertian{checker_texture}});
    hitables.push_back(new Sphere{Vec3f{0, 0, -1}, 0.5, new Lambertian{new PerlinTexture{}}});
    hitables.push_back(new Sphere{Vec3f{1, 0, -1}, 0.5, new Metal{0.8, 0.8, 0.0, 0.3}});
    hitables.push_back(new Sphere{Vec3f{-1, 0, -1}, 0.5, new Dielectric{1.5}});
    const std::string filepath = "/Users/lingtorp/Desktop/earth.jpg";
    hitables.push_back(new Sphere{Vec3f{0.75, 0.5, -1.0}, 0.25, new Lambertian{new ImageTexture{filepath}}});
    hitables.push_back(new Rectxy{3, 5, 1, 3, -2, new Emission{new ConstantTexture{4.0, 4.0, 4.0}}});
  }
  
  void two_spheres() {
    Texture* checker = new CheckerTexture{new ConstantTexture{Vec3f{0.2, 0.3, 0.1}}, new ConstantTexture{Vec3f{0.9, 0.9, 0.9}}};
    hitables.push_back(new Sphere{{0, -10, 0}, 10, new Lambertian{checker}});
    hitables.push_back(new Sphere{{0,  10, 0}, 10, new Lambertian{checker}});
  }
  
  void cornell_box() {
    Material* red = new Lambertian{new ConstantTexture{0.65, 0.05, 0.05}};
    Material* white = new Lambertian{new ConstantTexture{0.73, 0.73, 0.73}};
    Material* green = new Lambertian{new ConstantTexture{0.12, 0.45, 0.15}};
    Material* light = new Emission{new ConstantTexture{15, 15, 15}};
    hitables.push_back(new Rectyz{0, 555, 0, 555, 555, green});
    hitables.push_back(new Rectyz{0, 555, 0, 555, 0, red});
    hitables.push_back(new Rectxz{213, 343, 227, 332, 554, light});
    hitables.push_back(new Rectxz{0, 555, 0, 555, 0, white});
    hitables.push_back(new Rectxz{0, 555, 0, 555, 555, white});
  }
  
  void simple_light() {
    Texture* pertext = new PerlinTexture{};
    hitables.push_back(new Sphere{Vec3f{0, -1000, 0}, 1000, new Lambertian{pertext}});
    hitables.push_back(new Sphere{Vec3f{0, 2, 0}, 1, new Lambertian{pertext}});
    hitables.push_back(new Sphere{Vec3f{0, 2, 1}, 0.5, new Emission{new ConstantTexture{1, 1, 1}}});
    hitables.push_back(new Rectxy{3, 5, 1, 3, -2, new Emission{new ConstantTexture{0.94, 0.5, 0.95}}});
  }
};

#endif // RAYTRACER_WORLD_H
