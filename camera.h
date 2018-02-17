#ifndef RAYTRACER_CAMERA_H
#define RAYTRACER_CAMERA_H

#include "vector.h"

struct Camera {
  Vec3<> frame_anchor{};
  Vec3<> origin{0.0};
  // Direction vectors
  Vec3<> vertical{};
  Vec3<> horizontal{};
  
  Camera(double vfov, double aspect) {
    double theta = vfov * M_PI/180.0;
    double half_height = std::tan(theta / 2.0);
    double half_width = aspect * half_height;
    frame_anchor = {-half_width, -half_height, -1.0};
    horizontal = {2*half_width, 0.0, 0.0};
    vertical = {0.0, 2*half_height, 0.0};
  };
  
  Ray get_ray(double u, double v) { return Ray{origin, frame_anchor + u*horizontal + v*vertical - origin}; }
};

#endif // RAYTRACER_CAMERA_H
