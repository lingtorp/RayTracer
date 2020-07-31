#ifndef RAYTRACER_RANDOM_H
#define RAYTRACER_RANDOM_H

#include "vector.h"

#ifdef __APPLE__
#include <stdlib.h>

void init_random_generator() {}
#elif __LINUX__
// TODO: Implement
#include <cstdlib>

thread_local drand48_data* random_state;

void init_random_generator() {
  random_state = (struct drand48_data*) calloc(sizeof(drand48_data), 1);
}
#elif __WIN__
  #error "Not implemented yet ..."
#endif

/// Thread safe (no global state in the called library function) uniform double number function
inline double rand_s(double upper_bound) {
#ifdef __APPLE__
  return arc4random_uniform(upper_bound); // thread safe according to Apple source
#elif __LINUX__
  // Docs: https://linux.die.net/man/3/drand48_r
  double value;
  drand48_r(random_state, &value);
  return value;
#elif __WIN__
  #error "Not implemented yet ..."
#endif
}

/// Random float in interval [0, 1]
float rand_0_1() {
#ifdef __APPLE__
  return (float) drand48();
#elif __LINUX__
  #error "Not implemented yet ..."
#elif __WIN__
  #error "Not implemented yet ..."
#endif
}

Vec3f rand_in_unit_disc() {
  Vec3f v{0.0f};
  do {
    v = 2.0f * Vec3f{rand_0_1(), rand_0_1(), 0.0f} - Vec3f{1.0, 1.0, 0.0};
  } while (dot(v, v) >= 1.0f);
  return v;
}

Vec3f random_in_unit_sphere() {
  Vec3f p{0.0f};
  do {
    p = 2.0f * Vec3f{rand_0_1(), rand_0_1(), rand_0_1()} - Vec3f{1.0f};
  } while (p.squared_length() >= 1.0f);
  return p;
}

#endif // RAYTRACER_RANDOM_H
