#ifndef RAYTRACER_WORLD_H
#define RAYTRACER_WORLD_H

#include <vector>
#include "hitable.h"

struct World {
  std::vector<Hitable*> hitables;
  BVHNode bvh;
  
  World() = default;
  
  /// Called whenever the world should not change and just render
  void bake_world() {
    bvh = BVHNode{hitables.data(), hitables.size(), 0.0, 0.0};
  }
  
  bool hit(const Ray& r, double t_min, double t_max, Hit& hit) const {
    return bvh.hit(r, t_min, t_max, hit);
  }
};

#endif // RAYTRACER_WORLD_H
