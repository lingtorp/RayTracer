#ifndef RAYTRACER_RANDOM_H
#define RAYTRACER_RANDOM_H

#ifdef __APPLE__
#include <random>
#elif __LINUX__
// TODO: Implement
#else
#include <cstdlib>

thread_local drand48_data* random_state;

void init_random_generator() {
  random_state = (struct drand48_data*) calloc(sizeof(drand48_data), 1);
}
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
#endif
}

#endif // RAYTRACER_RANDOM_H
