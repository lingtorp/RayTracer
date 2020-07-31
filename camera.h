#ifndef RAYTRACER_CAMERA_H
#define RAYTRACER_CAMERA_H

#include "vector.h"
#include "random.h"

struct Camera {
  /// Frame anchor is lower left corner of the Camera 'sensor/view'
  Vec3f frame_anchor;
  Vec3f origin;
  // Direction vectors
  Vec3f vertical;
  Vec3f horizontal;
  Vec3f u, v, w;
  float lens_radius;
  
  Camera(Vec3f lookfrom, Vec3f lookat, float vfov, float aspect, float aperature, float focus_dist) {
    lens_radius = aperature / 2.0f;
    constexpr Vec3f vup = {0.0f, 1.0f, 0.0f};
    float theta = vfov * M_PI / 180.0f;
    float half_height = std::tan(theta / 2.0f);
    float half_width = aspect * half_height;
    origin = lookfrom;
    w = {lookfrom - lookat}; w.normalize();
    u = {cross(vup, w)}; u.normalize();
    v = cross(w, u);
    frame_anchor = origin - half_width*focus_dist*u - half_height*focus_dist*v - focus_dist*w;
    horizontal = 2.0f * half_width * focus_dist * u;
    vertical = 2.0f * half_height * focus_dist * v;
  };
  
  Ray get_ray(float s, float t) const {
    Vec3f rd = lens_radius*rand_in_unit_disc();
    Vec3f offset = u*rd.x + v*rd.y;
    return Ray{origin + offset, frame_anchor + s*horizontal + t*vertical - origin - offset};
  }
};

#endif // RAYTRACER_CAMERA_H
