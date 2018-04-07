#ifndef RAYTRACER_CAMERA_H
#define RAYTRACER_CAMERA_H

#include "vector.h"

Vec3<> rand_in_unit_disc() {
  Vec3<> v;
  do {
    v = 2.0*Vec3<>{drand48(), drand48(), 0.0} - Vec3<>{1.0, 1.0, 0.0};
  } while (dot(v, v) >= 1.0);
  return v;
}

struct Camera {
  /// Frame anchor is lower left corner of the Camera 'sensor/view'
  Vec3<> frame_anchor{};
  Vec3<> origin{0.0};
  // Direction vectors
  Vec3<> vertical{};
  Vec3<> horizontal{};
  Vec3<> u, v, w;
  double lens_radius;
  
  Camera(Vec3<> lookfrom, Vec3<> lookat, double vfov, double aspect, double aperature, double focus_dist) {
    lens_radius = aperature / 2.0;
    constexpr Vec3<> vup = {0.0, 1.0, 0.0};
    double theta = vfov * M_PI/180.0;
    double half_height = std::tan(theta / 2.0);
    double half_width = aspect * half_height;
    origin = lookfrom;
    w = {lookfrom - lookat}; w.normalize();
    u = {cross(vup, w)}; u.normalize();
    v = cross(w, u);
    frame_anchor = origin - half_width*focus_dist*u - half_height*focus_dist*v - focus_dist*w;
    horizontal = 2*half_width*focus_dist*u;
    vertical = 2*half_height*focus_dist*v;
  };
  
  Ray get_ray(double s, double t) {
    Vec3<> rd = lens_radius*rand_in_unit_disc();
    Vec3<> offset = u*rd.x + v*rd.y;
    return Ray{origin + offset, frame_anchor + s*horizontal + t*vertical - origin - offset};
  }
};

#endif // RAYTRACER_CAMERA_H
