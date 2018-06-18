#ifndef RAYTRACER_RAY_H
#define RAYTRACER_RAY_H

#include "vector.h"
#include <utility>

class Material;

struct AABB {
  Vec3<> min{};
  Vec3<> max{};
  AABB() = default;
  AABB(const Vec3<>& a, const Vec3<>& b): min(a), max(b) {}
  
  bool hit(const Ray& r, double tmin, double tmax) const {
    // Slab method in 3D
    for (size_t a = 0; a < 3; a++) {
      double inv_dir = 1.0 / r.direction()[a];
      double t0 = (min[a] - r.origin()[a]) * inv_dir;
      double t1 = (max[a] - r.origin()[a]) * inv_dir;
      if (inv_dir < 0.0) { std::swap(t0, t1); }
      tmin = t0 > tmin ? t0 : tmin;
      tmax = t1 < tmax ? t1 : tmax;
      if (tmax <= tmin) { return false; }
    }
    return true;
  }
  
  static AABB surrounding_box(AABB b0, AABB b1) {
    Vec3<> small{fmin(b0.min.x, b1.min.x), fmin(b0.min.y, b1.min.y), fmin(b0.min.z, b1.min.z)};
    Vec3<> big{fmax(b0.max.x, b1.max.x), fmax(b0.max.y, b1.max.y), fmax(b0.max.z, b1.max.z)};
    return AABB{small, big};
  }
};

struct Hit {
  double t = 0.0;
  Vec3<> p{}; // Position
  Vec3<> normal{};
  Material* mat = nullptr;
};

class Hitable {
public:
  virtual bool hit(const Ray& r, double t_min, double t_max, Hit& hit) const = 0;
  virtual bool bounding_box(double t0, double t1, AABB& box) const = 0;
};

static int box_x_compare(const void* a, const void* b) {
  AABB left, right;
  Hitable* ah = *(Hitable**)a;
  Hitable* bh = *(Hitable**)b;
  if (!ah->bounding_box(0, 0, left) || !bh->bounding_box(0, 0, right)) {
    std::cerr << "No bounding box in BVHNode constructor." << std::endl;
  }
  if (left.min.x - right.min.x < 0.0) {
    return -1;
  }
  return 1;
}

static int box_y_compare(const void* a, const void* b) {
  AABB left, right;
  Hitable* ah = *(Hitable**)a;
  Hitable* bh = *(Hitable**)b;
  if (!ah->bounding_box(0, 0, left) || !bh->bounding_box(0, 0, right)) {
    std::cerr << "No bounding box in BVHNode constructor." << std::endl;
  }
  if (left.min.y - right.min.y < 0.0) {
    return -1;
  }
  return 1;
}

static int box_z_compare(const void* a, const void* b) {
  AABB left, right;
  Hitable* ah = *(Hitable**)a;
  Hitable* bh = *(Hitable**)b;
  if (!ah->bounding_box(0, 0, left) || !bh->bounding_box(0, 0, right)) {
    std::cerr << "No bounding box in BVHNode constructor." << std::endl;
  }
  if (left.min.z - right.min.z < 0.0) {
    return -1;
  }
  return 1;
}

class BVHNode : public Hitable {
public:
  Hitable* left  = nullptr;
  Hitable* right = nullptr;
  AABB aabb;  // Box surrounding the nodes
  BVHNode(Hitable** l, uint64_t n, double time0, double time1) {
    int axis = int(3*rand_s(1.0));
    if (axis == 0) {
      qsort(l, n, sizeof(Hitable*), box_x_compare);
    } else if (axis == 1) {
      qsort(l, n, sizeof(Hitable*), box_y_compare);
    } else {
      qsort(l, n, sizeof(Hitable*), box_z_compare);
    }
    if (n == 1) {
      left = right = l[0];
    } else if (n == 2) {
      left = l[0]; right = l[1];
    } else {
      left = new BVHNode{l, n / 2, time0, time1};
      right = new BVHNode{l + n / 2, n - n / 2, time0, time1};
    }
    AABB box_left{}, box_right{};
    if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right)) {
      std::cerr << "Error during bounding volume hierarchy construction." << std::endl;
    }
    aabb = AABB::surrounding_box(box_left, box_right);
  };
  
  bool hit(const Ray &r, double t_min, double t_max, Hit& hit) const override {
    if (aabb.hit(r, t_min, t_max)) {
      Hit left_hit, right_hit;
      bool did_left_hit  = left->hit(r, t_min, t_max, left_hit);
      bool did_right_hit = right->hit(r, t_min, t_max, right_hit);
      if (did_left_hit && did_right_hit) {
        if (left_hit.t < right_hit.t) {
          hit = left_hit;
        } else {
          hit = right_hit;
        }
      } else if (did_left_hit) {
        hit = left_hit;
      } else if (did_right_hit) {
        hit = right_hit;
      } else {
        return false;
      }
      return true;
    }
    return false;
  }
  
  bool bounding_box(double t0, double t1, AABB& box) const override {
    box = aabb;
    return true;
  }
};

#endif // RAYTRACER_RAY_H
