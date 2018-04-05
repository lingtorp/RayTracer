#ifndef RAYTRACER_CAMERA_H
#define RAYTRACER_CAMERA_H

#include "vector.h"

struct Camera {
  /// Frame anchor is lower left corner of the Camera 'sensor/view'
  Vec3<> frame_anchor{};
  Vec3<> origin{0.0};
  // Direction vectors
  Vec3<> vertical{};
  Vec3<> horizontal{};
  
  Camera(Vec3<> lookfrom, Vec3<> lookat, Vec3<> vup, double vfov, double aspect) {
    double theta = vfov * M_PI/180.0;
    double half_height = std::tan(theta / 2.0);
    double half_width = aspect * half_height;
    origin = lookfrom;
    Vec3<> w = {lookfrom - lookat}; w.normalize();
    Vec3<> u = {cross(vup, w)}; u.normalize();
    Vec3<> v = cross(w, u);
    frame_anchor = origin - half_width*u - half_height*v - w;
    horizontal = 2*half_width*u;
    vertical = 2*half_height*v;
  };
  
  Ray get_ray(double u, double v) { return Ray{origin, frame_anchor + u*horizontal + v*vertical - origin}; }
};

#endif // RAYTRACER_CAMERA_H
