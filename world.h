#ifndef RAYTRACER_WORLD_H
#define RAYTRACER_WORLD_H

#include <vector>
#include "hitable.h"
#include "sphere.h"

struct World {
  std::vector<Hitable*> hitables;
  BVHNode* bvh;
  
  World() = default;
  
  /// Called whenever the world should not change and just render
  void bake_world() {
    delete bvh;
    bvh = new BVHNode{hitables.data(), hitables.size(), 0.0, 0.0};
  }
  
  bool hit(const Ray& r, double t_min, double t_max, Hit& hit) const {
    return bvh->hit(r, t_min, t_max, hit);
  }
  
  /// Scenes
  void default_scene() {
    auto checker_texture = new CheckerTexture{new ConstantTexture{0.2, 0.3, 0.1}, new ConstantTexture{0.9, 0.9, 0.9}};
    hitables.push_back(new Sphere{Vec3<>{0, -100.5, -1}, 100, new Lambertian{checker_texture}});
    hitables.push_back(new Sphere{Vec3<>{0, 0, -1}, 0.5, new Lambertian{new PerlinTexture{}}});
    hitables.push_back(new Sphere{Vec3<>{1, 0, -1}, 0.5, new Metal{0.8, 0.8, 0.0, 0.3}});
    hitables.push_back(new Sphere{Vec3<>{-1, 0, -1}, 0.5, new Dielectric{1.5}});
    hitables.push_back(new Sphere{Vec3<>{0.5, 0.5, -1.5}, 0.2, new Emission{new ConstantTexture{1.0, 1.0, 1.0}}});
    const std::string filepath = "/home/lingtorp/Desktop/earth.jpg";
    hitables.push_back(new Sphere{Vec3<>{0.75, 0.5, -1.0}, 0.25, new Lambertian{new ImageTexture{filepath}}});
  }
  
  void simple_light_scene() {
    // default_scene();
    Texture* pertext = new PerlinTexture{};
    hitables.push_back(new Sphere{Vec3<>{0, -1000, 0}, 1000, new Lambertian{pertext}});
    //hitables.push_back(new Sphere{Vec3<>{0, 2, 0}, 2, new Emission{new ConstantTexture{Vec3<>{1, 1, 1}}}});
    //hitables.push_back(new Sphere{Vec3<>{0, 7, 0}, 2, new Emission{new ConstantTexture{Vec3<>{1, 1, 1}}}});
    hitables.push_back(new Rectxy{3, 5, 1, 3, -2, new Emission{new ConstantTexture{4, 4, 4}}});
  }
};

#endif // RAYTRACER_WORLD_H
