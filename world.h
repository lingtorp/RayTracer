#ifndef RAYTRACER_WORLD_H
#define RAYTRACER_WORLD_H

#include <vector>
#include "hitable.h"

struct World {
  std::vector<Hitable*> hitables;
  
  World() = default;
  
  bool hit(const Ray& r, double t_min, double t_max, Hit& hit) const {
    bool hit_anything = false;
    double closest_so_far = t_max;
    for (auto& hitable : hitables) {
      if (hitable->hit(r, t_min, closest_so_far, hit)) {
        closest_so_far = hit.t;
        hit_anything = true;
      }
    }
    return hit_anything;
  }
};

#endif // RAYTRACER_WORLD_H
