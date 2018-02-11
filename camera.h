#ifndef RAYTRACER_CAMERA_H
#define RAYTRACER_CAMERA_H

#include "vector.h"

struct Camera {
  Vec3<> frame_anchor;
  Vec3<> origin;
  // Direction vectors
  Vec3<> vertical;
  Vec3<> horizontal;
  
  Camera(): frame_anchor{-2.0, -1.0, -1.0}, origin{}, vertical{0.0, 2.0, 0.0}, horizontal{4.0, 0.0, 0.0} {};
  
  Ray get_ray(double u, double v) { return Ray{origin, frame_anchor + u*horizontal + v*vertical - origin}; }
};

#endif // RAYTRACER_CAMERA_H
