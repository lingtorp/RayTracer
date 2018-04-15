#ifndef RAYTRACER_RANDOM_H
#define RAYTRACER_RANDOM_H

#ifdef __APPLE__
#include <random>
#elif __LINUX__
// TODO: Implement
#else
exit(1); // Unsupported
#endif

/// Thread safe (no global state in the called library function) uniform double number function
inline double rand_s(double upper_bound) {
#ifdef __APPLE__
  return arc4random_uniform(upper_bound); // thread safe according to Apple source
#elif __LINUX__
  // Docs: https://linux.die.net/man/3/drand48_r
  exit(1); // TODO: Implement with int drand48_r(struct drand48_data *buffer, double *result);
#endif
}

#endif // RAYTRACER_RANDOM_H
