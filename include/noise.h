#ifndef NOISE_H
#define NOISE_H

#include <random>
#include <algorithm>
#include "../vector.h"
#include <iostream>
#include <cstdint>
#include <array>

// TODO: Color mapping module (noise + color scheme --> vector or smt)
// TODO: Doxygen style comments

/**
 * Base class for noise generating classes
 */
class Noise {
public:
  /// 2D raw noise from the underlying noise algorithm
  virtual float get_value(float x, float y) const = 0;
  
  /// 3D raw noise from the underlying noise algorithm
  virtual float get_value(float x, float y, float z) const = 0;
  
  // FIXME: Is turbulence like defined here really from the original Perlin patent?
  // FIXME: Is it a visually useful effect?
  /// 3D turbulence noise which simulates fBm
  float turbulence(float x, float y, float zoom_factor) const {
    float value = 0.0f;
    float zoom = zoom_factor;
    while (zoom >= 1.0f) {
      value += std::abs(get_value(x / zoom, y / zoom) * zoom);
      zoom /= 2.0f;
    }
    return value / zoom_factor;
  }
  
  // FIXME: Is turbulence like defined here really from the original Perlin patent?
  // FIXME: Is it a visually useful effect?
  /// 3D turbulence noise which simulates fBm
  /// Reference: http://lodev.org/cgtutor/randomnoise.html & orignal Perlin noise paper
  float turbulence(float x, float y, float z, float zoom_factor) const {
    float value = 0.0f;
    float zoom = zoom_factor;
    while (zoom >= 1.0f) {
      value += std::abs(get_value(x / zoom, y / zoom, z / zoom) * zoom);
      zoom /= 2.0f;
    }
    return value / zoom_factor;
  }
  
  /// 2D turbulence noise which simulates fBm
  float fbm(Vec2f v, float zoom_factor) const {
    return fbm(v.x, v.y, zoom_factor);
  }
  
  /// 2D turbulence noise which simulates fBm
  float fbm(float x, float y, float zoom_factor) const {
    float value = 0.0f;
    float zoom = zoom_factor;
    while (zoom >= 1.0f) {
      value += get_value(x / zoom, y / zoom) * zoom;
      zoom /= 2.0f;
    }
    return value / zoom_factor;
  }
  
  /// 3D turbulence noise which simulates fBm
  float fbm(Vec3f v, float zoom_factor) const {
    return fbm(v.x, v.y, v.z, zoom_factor);
  }
  
  /// 3D turbulence noise which simulates fBm
  float fbm(float x, float y, float z, float zoom_factor) const {
    float value = 0.0f;
    float zoom = zoom_factor;
    while (zoom >= 1.0f) {
      value += get_value(x / zoom, y / zoom, z / zoom) * zoom;
      zoom /= 2.0f;
    }
    return value / zoom_factor;
  }
  
  /// 3D Billowy turbulence
  float turbulence_billowy(float x, float y, float z, float zoom_factor) const {
    float value = 0.0f;
    float zoom = zoom_factor;
    while (zoom >= 1.0f) {
      value += std::abs(get_value(x / zoom, y / zoom, z / zoom) * zoom);
      zoom /= 2.0f;
    }
    return value / zoom_factor;
  }
  
  /// 3D Ridged turbulence
  float turbulence_ridged(float x, float y, float z, float zoom_factor) const {
    float value = 0.0f;
    float zoom = zoom_factor;
    while (zoom >= 1.0f) {
      value += (1.0f - std::abs(get_value(x / zoom, y / zoom, z / zoom) * zoom));
      zoom /= 2.0f;
    }
    return value / zoom_factor;
  }
  
  // FIXME: Octaves, is the implementation correct?
  // FIXME: Visually pleasing effect?
  /// 2D fractional Brownian motion noise of the underlying noise algorithm
  float octaves(float x, float y, uint32_t octaves, float persistance = 1.0f, float amplitude = 1.0) const {
    float total = 0.0f;
    float max_value = 0.0f;
    float frequency = 1.0f;
    for (size_t i = 0; i < octaves; ++i) {
      total += get_value(x / frequency, y / frequency) * amplitude;
      max_value += amplitude;
      
      amplitude *= persistance;
      frequency *= 2.0f;
    }
    
    // Dividing by the max amplitude sum brings it into [-1, 1] range
    return total / max_value;
  }
  
  /// 3D fractional Brownian motion noise of the underlying noise algorithm
  float octaves(float x, float y, float z, uint32_t octaves, float persistance = 1.0f, float amplitude = 1.0) const {
    float total = 0.0f;
    float max_value = 0.0f;
    float frequency = 1.0f;
    for (size_t i = 0; i < octaves; ++i) {
      total += get_value(x / frequency, y / frequency, z / frequency) * amplitude;
      max_value += amplitude;
      
      amplitude *= persistance;
      frequency *= 2;
    }
    
    // Dividing by the max amplitude sum brings it into [-1, 1] range
    return total / max_value;
  }
  
  /// 3D fractional Brownian motion noise in which each octave gets its own amplitude
  float octaves(float x, float y, float z, const std::vector<float> &amplitudes) const {
    float total = 0.0f;
    float max_value = 0.0f;
    float frequency = 1.0f;
    for (float amplitude : amplitudes) {
      total += get_value(x / frequency, y / frequency, z / frequency) * amplitude;
      max_value += amplitude;
      frequency *= 2;
    }
    
    // Dividing by the max amplitude sum brings it into [-1, 1] range
    return total / max_value;
  }
  
  /// Warps the domain of the noise function creating more natural looking features
  float domain_wrapping(float x, float y, float z, float scale) const {
    Vec3f p{x, y, z};
    Vec3f offset{50.2f, 10.3f, 10.5f};
    
    Vec3f q{fbm(p + offset, scale), fbm(p + offset, scale), fbm(p + offset, scale)};
    Vec3f qq{100.0f*q.x, 100.0f*q.y, 100.0f*q.z};
    
    /// Adjusting the scales in r makes a cool ripple effect through the noise
    Vec3f r{fbm(p + qq + Vec3f{1.7f, 9.2f, 5.1f}, scale * 1.0f),
            fbm(p + qq + Vec3f{8.3f, 2.8f, 2.5f}, scale * 1.0f),
            fbm(p + qq + Vec3f{1.2f, 6.9f, 8.4f}, scale * 1.0f)};
    Vec3f rr{100.0f*r.x, 100.0f*r.y, 100.0f*r.z};
    
    return fbm(p + rr, scale);
  }

protected:
  static inline float clamp(float in, float lo, float hi) {
    return std::max(lo, std::min(hi, in));
  }
  
  // TODO: Document
  static inline float smoothstep(float t) { return t * t * (3.0f - 2.0f * t); }
  // TODO: Document
  static inline float quintic_fade(float t) { return t * t * t * (t * (t * 6.0f - 15.0f) + 10.0f); }
  /// Linear interpolation between a and b with t as a variable
  static inline float lerp(float t, float a, float b) { return (1.0f - t) * a + t * b; }
};

namespace Simplex {
  /**
   * Simplex noise/Improved Perlin noise from the 'Improved noise' patent
   * - Gradient creation on-the-fly using bit manipulation
   * - Gradient selection uses bit manipulation (from the above point)
   */
  class Patent : public Noise {
    /// Bit patterns for the creation of the gradients
    std::vector<u_char> bit_patterns;
    
    /// Returns the n'th bit of num
    inline u_char bit(int num, int n) const {
      return (u_char) ((num >> n) & 0b1);
    }
  
  public:
    explicit Patent(uint64_t seed) : bit_patterns{0x15, 0x38, 0x32, 0x2C, 0x0D, 0x13, 0x07, 0x2A} {}
    
    /********************************** Simplex 2D Noise **********************************/
    
    /// Skews the coordinate to normal Euclidean coordinate system
    Vec2f skew(Vec2f v) const {
      const float F = (std::sqrt(1.0f + 2.0f) - 1.0f) / 2.0f;
      float s = (v.x + v.y) * F;
      return {v.x + s, v.y + s};
    }
    
    /// Unskews the coordinate back to the simpletic coordinate system
    Vec2f unskew(Vec2f v) const {
      const float G = (1.0f - (1.0f / sqrt(2.0f + 1.0f))) / 2.0f;
      float s = (v.x + v.y) * G;
      return {v.x - s, v.y - s};
    }
    
    /// Given a coordinate (i, j) selects the B'th bit
    u_char b(int i, int j, int B) const {
      auto bit_index = 2 * (i & (0b1 << B)) + (j & (0b1 << B));
      return bit_patterns[bit_index];
    }
    
    /// Given a coordinate (i, j) generates a gradient vector
    Vec2f grad(int i, int j) const {
      uint32_t bit_sum = b(i, j, 0) + b(j, i, 1) + b(i, j, 2) + b(j, i, 3);
      float u = (bit_sum & 0b01) ? 1.0f : 0.0f;
      float v = (bit_sum & 0b10) ? 1.0f : 0.0f;
      u = (bit_sum & 0b1000) ? -u : u;
      v = (bit_sum & 0b0100) ? -v : v;
      return {u, v};
    }
    
    // FIXME: Floatcheck against the patent
    float get_value(float x, float y) const override {
      /// Skew
      const float F = (std::sqrt(2.0f + 1.0f) - 1.0f) / 2.0f;
      float s = (x + y) * F;
      float xs = x + s;
      float ys = y + s;
      int i = (int) std::floor(xs);
      int j = (int) std::floor(ys);
      
      /// Unskew - find first vertex of the simplex
      const float G = (3.0f - std::sqrt(2.0f + 1.0f)) / 6.0f;
      float t = (i + j) * G;
      Vec2f cell_origin{i - t, j - t};
      Vec2f vertex_a = Vec2f{x, y} - cell_origin;
      
      // Figure out which vertex is next
      float x_step = 0.0f;
      float y_step = 0.0f;
      if (vertex_a.x > vertex_a.y) { // Lower triangle
        x_step = 1.0f;
      } else {
        y_step = 1.0f;
      }
      
      // A change of one unit step is; x = x' + (x' + y') * G <--> x = 1.0f + (1.0 + 1.0) * G <--> x = 1.0 + 2.0 * G
      Vec2f vertex_b{vertex_a.x - x_step + G, vertex_a.y - y_step + G};
      Vec2f vertex_c{vertex_a.x - 1.0f + 2.0f * G, vertex_a.y - 1.0f + 2.0f * G};
      
      auto grad_a = grad(i, j);
      auto grad_b = grad(i + x_step, j + y_step);
      auto grad_c = grad(i + 1, j + 1);
      
      /// Calculate contribution from the vertices in a circle
      // max(0, r^2 - d^2)^4 * gradient.dot(vertex)
      const float radius = 0.6f * 0.6f; // Radius of the surflet circle (0.6 in patent)
      float sum = 0.0f;
      
      float t0 = radius - vertex_a.length() * vertex_a.length();
      if (t0 > 0.0f) {
        sum += std::pow(t0, 4.0f) * grad_a.dot(vertex_a);
      }
      
      float t1 = radius - vertex_b.length() * vertex_b.length();
      if (t1 > 0.0f) {
        sum += std::pow(t1, 4.0f) * grad_b.dot(vertex_b);
      }
      
      float t2 = radius - vertex_c.length() * vertex_c.length();
      if (t2 > 0.0f) {
        sum += std::pow(t2, 4.0f) * grad_c.dot(vertex_c);
      }
      
      return 220.0f * sum;
    }
    
    /********************************** Simplex 3D Noise **********************************/
    
    /// Hashes a coordinate (i, j, k) then selects one of the bit patterns
    u_char b(int i, int j, int k, int B) const {
      u_char bit_index = bit(i, B) << 2 | bit(j, B) << 2 | bit(k, B);
      return bit_patterns[bit_index];
    }
    
    /**
     * Generates the gradient vector from the vertex and the relative vector.
     * @param vertex Skewed coordinate of the vertex, used to generate the bit sum.
     * @param rel Relative vector of (x, y, z) and the vertex in the unskewed coordinate system.
     * @return Gradient vector
     */
    Vec3f grad(Vec3f vertex, Vec3f rel) const {
      int i = (int) vertex.x;
      int j = (int) vertex.y;
      int k = (int) vertex.z;
      int sum = b(i, j, k, 0) + b(j, k, i, 1) + b(k, i, j, 2) + b(i, j, k, 3) + b(j, k, i, 4) + b(k, i, j, 5) +
                b(i, j, k, 6) + b(j, k, i, 7);
      
      // Magnitude computation based on the three lower bits of the bit sum
      Vec3f pqr = rel;
      if (bit(sum, 0) == !bit(sum, 1)) { // xor on bit 0, 1 --> rotation and zeroing
        if (bit(sum, 0)) { // Rotation
          pqr.x = rel.y;
          pqr.y = rel.z;
          pqr.z = rel.x;
        } else {
          pqr.x = rel.z;
          pqr.y = rel.x;
          pqr.z = rel.y;
        }
        
        if (bit(sum, 2)) { // Zeroing out
          pqr.y = 0.0f;
        } else {
          pqr.z = 0.0f;
        }
      } else if (bit(sum, 0) && bit(sum, 1)) {
        if (bit(sum, 2)) { // Zeroing out
          pqr.y = 0.0f;
        } else {
          pqr.z = 0.0f;
        }
      }
      
      // Octant computation based on the three upper bits of the bit sum
      if (bit(sum, 5) == bit(sum, 3)) { pqr.x = -pqr.x; }
      if (bit(sum, 5) == bit(sum, 4)) { pqr.y = -pqr.y; }
      if (bit(sum, 5) != (bit(sum, 4) == !bit(sum, 3))) { pqr.z = -pqr.z; }
      
      return pqr;
    }
    
    /// Skews the coordinate to normal Euclidean coordinate system
    Vec3f skew(Vec3f v) const {
      const float F = (std::sqrt(1.0f + 3.0f) - 1.0f) / 3.0f;
      float s = (v.x + v.y + v.z) * F;
      return {v.x + s, v.y + s, v.z + s};
    }
    
    /// Unskews the coordinate back to the simpletic coordinate system
    Vec3f unskew(Vec3f v) const {
      const float G = (1.0f - (1.0f / sqrt(3.0f + 1.0f))) / 3.0f;
      float s = (v.x + v.y + v.z) * G;
      return {v.x - s, v.y - s, v.z - s};
    }
    
    /**
     * Computes the spherical kernel contribution from one vertex in a simpletic cell offseted by the vector ijk.
     * @param uvw Position within the simplex cell (unskewed)
     * @param ijk First vertex in the simplex cell (unskewed)
     * @param vertex Vertex in the unit simplex cell (unskewed)
     * @return Contribution from the vertex
     */
    float kernel(Vec3f uvw, Vec3f ijk, Vec3f vertex) const {
      float sum = 0.0f;
      Vec3f rel = uvw - vertex; // Relative simplex cell vertex
      float t = 0.6f - rel.length() * rel.length(); // 0.6 - x*x - y*y - z*z
      if (t > 0.0f) {
        Vec3f pqr = grad(ijk + vertex, rel); // Generate gradient vector for vertex
        t *= t;
        sum += 8.0f * t * t * pqr.sum();
      }
      return sum;
    }
    
    float get_value(float x, float y, float z) const override {
      /// Skew in the coordinate to the euclidean coordinate system
      Vec3f xyz = {x, y, z};
      Vec3f xyzs = skew(xyz);
      /// Skewed unit simplex cell
      Vec3f ijks = xyzs.floor(); // First vertex in euclidean coordinates
      Vec3f ijk = unskew(ijks); // First vertex in the simpletic cell
      
      /// Finding the traversal order of vertices of the unit simplex in which (x,y,z) is in.
      Vec3f uvw = xyz - ijk; // Relative unit simplex cell origin
      std::array<Vec3f, 4> vertices; // n + 1 is the number of vertices in a n-dim. simplex
      vertices[0] = unskew({0.0f, 0.0f, 0.0f});
      if (uvw.x > uvw.y) {
        if (uvw.y > uvw.z) {
          // u, v, w
          vertices[1] = unskew({1.0f, 0.0f, 0.0f});
          vertices[2] = unskew({1.0f, 1.0f, 0.0f});
        } else {
          if (uvw.x > uvw.z) {
            // u, w, v
            vertices[1] = unskew({1.0f, 0.0f, 0.0f});
            vertices[2] = unskew({1.0f, 0.0f, 1.0f});
          } else {
            // w, u, v
            vertices[1] = unskew({0.0f, 0.0f, 1.0f});
            vertices[2] = unskew({1.0f, 0.0f, 1.0f});
          }
        }
      } else {
        if (uvw.y > uvw.z) {
          if (uvw.z > uvw.x) {
            // v, w, u
            vertices[1] = unskew({0.0f, 1.0f, 0.0f});
            vertices[2] = unskew({0.0f, 1.0f, 1.0f});
          } else {
            // v, u, w
            vertices[1] = unskew({0.0f, 1.0f, 0.0f});
            vertices[2] = unskew({1.0f, 1.0f, 0.0f});
          }
        } else {
          // w, v, u
          vertices[1] = unskew({0.0f, 0.0f, 1.0f});
          vertices[2] = unskew({0.0f, 1.0f, 1.0f});
        }
      }
      vertices[3] = unskew({1.0f, 1.0, 1.0});
      
      /// Spherical kernel summation - contribution from each vertex
      float sum = kernel(uvw, ijk, vertices[0]) + kernel(uvw, ijk, vertices[1]) +
                  kernel(uvw, ijk, vertices[2]) + kernel(uvw, ijk, vertices[3]);
      
      return clamp(sum, -1.0f, 1.0f);
    }
  };
  
  /**
   * Simplex noise implementation using the a hybrid approach with features from simplex & Perlin noise
   * - Gradient table instead of on-the-fly gradient creation, as the original Perlin noise algorithm
   * - Permutation table instead of bit manipulation, unlike the patented algorithm
   * - Using modulo hashing to select the gradients via the permutation table
   */
  template<int num_grads = 256>
  class Tables : public Noise {
    /// 2D Normalized gradients table
    std::array<Vec2f, num_grads> grads2;
    
    /// 3D Normalized gradients table
    std::array<Vec3f, num_grads> grads3;
    
    /// Permutation table for indices to the gradients
    std::array<u_char, num_grads> perms;
  public:
    /// Perms size is float that of grad to avoid index wrapping
    explicit Tables(uint64_t seed) {
      std::mt19937 engine(seed);
      std::uniform_real_distribution<float> distr(-1.0f, 1.0f);
      /// Fill the gradients list with random normalized vectors
      for (int i = 0; i < grads2.size(); i++) {
        float x = distr(engine);
        float y = distr(engine);
        float z = distr(engine);
        Vec2f grad_vector = Vec2f{x, y}.normalized();
        grads2[i] = grad_vector;
        Vec3f grad3_vector = Vec3f{x, y, z}.normalized();
        grads3[i] = grad3_vector;
      }
      
      /// Fill gradient lookup array with random indices to the gradients list
      /// Fill with indices from 0 to perms.size()
      std::iota(perms.begin(), perms.end(), 0);
      
      /// Randomize the order of the indices
      std::shuffle(perms.begin(), perms.end(), engine);
    }
    
    float get_value(float x, float y) const override {
      const float F = (std::sqrt(2.0f + 1.0f) - 1.0f) / 2.0f; // F = (sqrt(n + 1) - 1) / n
      float s = (x + y) * F;
      float xs = x + s;
      float ys = y + s;
      int i = (int) std::floor(xs);
      int j = (int) std::floor(ys);
      
      const float G = (3.0f - std::sqrt(2.0f + 1.0f)) / 6.0f; // G = (1 - (1 / sqrt(n + 1)) / n
      float t = (i + j) * G;
      Vec2f cell_origin{i - t, j - t};
      Vec2f vertex_a = Vec2f{x, y} - cell_origin;
      
      float x_step = 0.0f;
      float y_step = 0.0f;
      if (vertex_a.x > vertex_a.y) { // Lower triangle
        x_step = 1.0f;
      } else {
        y_step = 1.0f;
      }
      
      Vec2f vertex_b{vertex_a.x - x_step + G, vertex_a.y - y_step + G};
      Vec2f vertex_c{vertex_a.x - 1.0f + 2.0f * G, vertex_a.y - 1.0f + 2.0f * G};
      
      const uint32_t ii = i % 255; // FIXME: Bit mask instead? Measure speedup
      const uint32_t jj = j % 255;
      Vec2f grad_a = grads2[perms[ii + perms[jj]]];
      Vec2f grad_b = grads2[perms[ii + x_step + perms[jj + y_step]]];
      Vec2f grad_c = grads2[perms[ii + 1 + perms[jj + 1]]];
      
      /// Calculate contribution from the vertices in a circle
      const float radius = 0.6f; // Radius of the surflet circle (0.6 in patent)
      float sum = 0.0f;
      
      float t0 = radius - vertex_a.length() * vertex_a.length();
      if (t0 > 0.0f) {
        sum += 8.0f * std::pow(t0, 4.0f) * grad_a.dot(vertex_a);
      }
      
      float t1 = radius - vertex_b.length() * vertex_b.length();
      if (t1 > 0.0f) {
        sum += 8.0f * std::pow(t1, 4.0f) * grad_b.dot(vertex_b);
      }
      
      float t2 = radius - vertex_c.length() * vertex_c.length();
      if (t2 > 0.0f) {
        sum += 8.0f * std::pow(t2, 4.0f) * grad_c.dot(vertex_c);
      }
      
      return clamp(sum, -1.0f, 1.0f);
    }
    
    // TODO: Implement
    float get_value(float x, float y, float z) const override { exit(EXIT_FAILURE); }
  };
}

namespace Perlin {
  /**
   * Improved Perlin noise from 2002
   * ACM: http://dl.acm.org/citation.cfm?id=566636
   *
   * The improvements made by Kenneth Perlin were twofold;
   * 1) Randomly generated gradients changed to static gradients
   * 2) Changed interpolation function
   */
  template<int num_grads = 256>
  class Improved : public Noise {
    /// 2D Normalized gradients table
    std::array<Vec2f, 4> grads;
    
    /// 3D Normalized gradients table
    std::array<Vec3f, 16> grads3;
    
    /// Permutation table for indices to the gradients (2D)
    std::array<u_char, num_grads> perms;
    
    /// Permutation table for indices to the gradients (3D)
    std::array<u_char, num_grads> perms3;
  
  public:
    explicit Improved(uint64_t seed) {
      std::mt19937 engine(seed);
      std::uniform_real_distribution<float> distr(-1.0f, 1.0);
      /// 4 gradients for each edge of a unit square, no need for padding, is power of 2
      grads = {
              Vec2f{ 1.0f,  0.0f},
              Vec2f{ 0.0f,  1.0f},
              Vec2f{-1.0f,  0.0f},
              Vec2f{ 0.0f, -1.0f}
      };
      // FIXME: Is all of the vectors inside grads?
      /// 12 gradients from the center to each edge of a unit cube, 4 duplicated vectors for padding so that the modulo is on a power of 2 (faster)
      grads3 = {
              Vec3f{ 1.0f,  1.0f,  0.0f},
              Vec3f{-1.0f,  1.0f,  0.0f},
              Vec3f{ 1.0f, -1.0f,  0.0f},
              Vec3f{-1.0f, -1.0f,  0.0f},
              Vec3f{ 1.0f,  0.0f,  1.0f},
              Vec3f{-1.0f,  0.0f,  1.0f},
              Vec3f{ 1.0f,  0.0f, -1.0f},
              Vec3f{-1.0f,  0.0f, -1.0f},
              Vec3f{ 0.0f,  1.0f,  1.0f},
              Vec3f{ 0.0f, -1.0f,  1.0f},
              Vec3f{ 0.0f,  1.0f, -1.0f},
              Vec3f{ 0.0f, -1.0f, -1.0f},
              Vec3f{ 1.0f,  1.0f,  0.0f},
              Vec3f{-1.0f,  1.0f,  0.0f},
              Vec3f{ 0.0f, -1.0f,  1.0f},
              Vec3f{ 0.0f, -1.0f, -1.0f}
      };
      
      /// Fill gradient lookup array with random indices to the gradients list
      for (size_t i = 0; i < perms.size(); i++) { perms[i] = i % grads.size(); }
      for (size_t i = 0; i < perms3.size(); i++) { perms3[i] = i % grads3.size(); }
      
      /// Randomize the order of the indices
      std::shuffle(perms.begin(), perms.end(), engine);
    }
    
    float get_value(float X, float Y) const override {
      /// Compress the coordinates inside the chunk; float part + int part = point coordinate
      X += 0.1f;
      Y += 0.1f; // Skew coordinates to avoid integer lines becoming zero
      /// Grid points from the chunk in the world
      int X0 = (int) std::floor(X);
      int Y0 = (int) std::floor(Y);
      int X1 = (int) std::ceil(X);
      int Y1 = (int) std::ceil(Y);
      
      /// Gradients using hashed indices from lookup list
      Vec2f x0y0 = grads[perms[(X0 + perms[Y0 % perms.size()]) % perms.size()]];
      Vec2f x1y0 = grads[perms[(X1 + perms[Y0 % perms.size()]) % perms.size()]];
      Vec2f x0y1 = grads[perms[(X0 + perms[Y1 % perms.size()]) % perms.size()]];
      Vec2f x1y1 = grads[perms[(X1 + perms[Y1 % perms.size()]) % perms.size()]];
      
      /// Vectors from gradients to point in unit square
      Vec2f v00 = Vec2f{X - X0, Y - Y0};
      Vec2f v10 = Vec2f{X - X1, Y - Y0};
      Vec2f v01 = Vec2f{X - X0, Y - Y1};
      Vec2f v11 = Vec2f{X - X1, Y - Y1};
      
      /// Contribution of gradient vectors by dot product between relative vectors and gradients
      float d00 = x0y0.dot(v00);
      float d10 = x1y0.dot(v10);
      float d01 = x0y1.dot(v01);
      float d11 = x1y1.dot(v11);
      
      /// Interpolate dot product values at sample point using polynomial interpolation 6x^5 - 15x^4 + 10x^3
      float yf = Y - Y0; // Float offset inside the square [0, 1]
      float xf = X - X0; // Float offset inside the square [0, 1]
      
      float wx = quintic_fade(xf);
      float wy = quintic_fade(yf);
      
      /// Interpolate along x for the contributions from each of the gradients
      float xa = lerp(wx, d00, d10);
      float xb = lerp(wx, d01, d11);
      
      float val = lerp(wy, xa, xb);
      
      return clamp(val, -1.0f, 1.0);
    }
    
    float get_value(float X, float Y, float Z) const override {
      /// Compress the coordinates inside the chunk; float part + int part = point coordinate
      /// Grid points from the chunk in the world
      int X0 = (int) std::floor(X);
      int Y0 = (int) std::floor(Y);
      int X1 = (int) std::ceil(X);
      int Y1 = (int) std::ceil(Y);
      int Z0 = (int) std::floor(Z);
      int Z1 = (int) std::ceil(Z);
      
      /// Gradients using hashed indices from lookup list
      Vec3f x0y0z0 = grads3[perms[(X0 + perms[(Y0 + perms[Z0 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x1y0z0 = grads3[perms[(X1 + perms[(Y0 + perms[Z0 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x0y1z0 = grads3[perms[(X0 + perms[(Y1 + perms[Z0 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x1y1z0 = grads3[perms[(X1 + perms[(Y1 + perms[Z0 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      
      Vec3f x0y0z1 = grads3[perms[(X0 + perms[(Y0 + perms[Z1 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x1y0z1 = grads3[perms[(X1 + perms[(Y0 + perms[Z1 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x0y1z1 = grads3[perms[(X0 + perms[(Y1 + perms[Z1 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x1y1z1 = grads3[perms[(X1 + perms[(Y1 + perms[Z1 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      
      /// Vectors from gradients to point in unit cube
      Vec3f v000 = Vec3f{X - X0, Y - Y0, Z - Z0};
      Vec3f v100 = Vec3f{X - X1, Y - Y0, Z - Z0};
      Vec3f v010 = Vec3f{X - X0, Y - Y1, Z - Z0};
      Vec3f v110 = Vec3f{X - X1, Y - Y1, Z - Z0};
      
      Vec3f v001 = Vec3f{X - X0, Y - Y0, Z - Z1};
      Vec3f v101 = Vec3f{X - X1, Y - Y0, Z - Z1};
      Vec3f v011 = Vec3f{X - X0, Y - Y1, Z - Z1};
      Vec3f v111 = Vec3f{X - X1, Y - Y1, Z - Z1};
      
      /// Contribution of gradient vectors by dot product between relative vectors and gradients
      float d000 = dot(x0y0z0, v000);
      float d100 = dot(x1y0z0, v100);
      float d010 = dot(x0y1z0, v010);
      float d110 = dot(x1y1z0, v110);
      
      float d001 = dot(x0y0z1, v001);
      float d101 = dot(x1y0z1, v101);
      float d011 = dot(x0y1z1, v011);
      float d111 = dot(x1y1z1, v111);
      
      /// Interpolate dot product values at sample point using polynomial interpolation 6x^5 - 15x^4 + 10x^3
      float yf = Y - Y0; // Float offset inside the cube [0, 1]
      float xf = X - X0; // Float offset inside the cube [0, 1]
      float zf = Z - Z0; // Float offset inside the cube [0, 1]
      
      float wx = quintic_fade(xf);
      float wy = quintic_fade(yf);
      float wz = quintic_fade(zf);
      
      /// Interpolate along x for the contributions from each of the gradients
      float xa = lerp(wx, d000, d100);
      float xb = lerp(wx, d010, d110);
      
      float xc = lerp(wx, d001, d101);
      float xd = lerp(wx, d011, d111);
      
      /// Interpolate along y for the contributions from each of the gradients
      float ya = lerp(wy, xa, xb);
      float yb = lerp(wy, xc, xd);
      
      /// Interpolate along z for the contributions from each of the gradients
      float za = lerp(wz, ya, yb);
      
      return clamp(za, -1.0f, 1.0f);
    }
  };
  
  /**
   * Original Perlin noise from 1985
   * ACM: http://dl.acm.org/citation.cfm?id=325247&CFID=927914208&CFTOKEN=31672107
   */
  class Original : public Noise {
    /// 2D Normalized gradients table
    std::vector<Vec2f> grads;
    
    /// 3D Normalized gradients table
    std::vector<Vec3f> grads3;
    
    /// Permutation table for indices to the gradients
    std::vector<u_char> perms;
  
  public:
    Original(uint64_t seed) : grads(256), grads3(256), perms(256) {
      std::mt19937 engine(seed);
      std::uniform_real_distribution<float> distr(-1.0f, 1.0f);
      /// Fill the gradients list with random normalized vectors
      for (size_t i = 0; i < grads.size(); i++) {
        float x = distr(engine);
        float y = distr(engine);
        float z = distr(engine);
        Vec2f grad_vector = Vec2f{x, y}.normalized();
        grads[i] = grad_vector;
        Vec3f grad3_vector = Vec3f{x, y, z}.normalized();
        grads3[i] = grad3_vector;
      }
      
      /// Fill gradient lookup array with random indices to the gradients list
      /// Fill with indices from 0 to perms.size()
      std::iota(perms.begin(), perms.end(), 0);
      
      /// Randomize the order of the indices
      std::shuffle(perms.begin(), perms.end(), engine);
    }
    
    float get_value(float X, float Y) const override {
      /// Compress the coordinates inside the chunk; float part + int part = point coordinate
      X += 0.1f;
      Y += 0.1f; // Skew coordinates to avoid integer lines becoming zero
      /// Grid points from the chunk in the world
      int X0 = (int) std::floor(X);
      int Y0 = (int) std::floor(Y);
      int X1 = (int) std::ceil(X);
      int Y1 = (int) std::ceil(Y);
      
      /// Gradients using hashed indices from lookup list
      // FIXME: Implement variation where perms.size() is a power of two in order to do a bit masking instead, measure speedup.
      Vec2f x0y0 = grads[perms[(X0 + perms[Y0 % perms.size()]) % perms.size()]];
      Vec2f x1y0 = grads[perms[(X1 + perms[Y0 % perms.size()]) % perms.size()]];
      Vec2f x0y1 = grads[perms[(X0 + perms[Y1 % perms.size()]) % perms.size()]];
      Vec2f x1y1 = grads[perms[(X1 + perms[Y1 % perms.size()]) % perms.size()]];
      
      /// Vectors from gradients to point in unit square
      Vec2f v00 = Vec2f{X - X0, Y - Y0};
      Vec2f v10 = Vec2f{X - X1, Y - Y0};
      Vec2f v01 = Vec2f{X - X0, Y - Y1};
      Vec2f v11 = Vec2f{X - X1, Y - Y1};
      
      /// Contribution of gradient vectors by dot product between relative vectors and gradients
      float d00 = x0y0.dot(v00);
      float d10 = x1y0.dot(v10);
      float d01 = x0y1.dot(v01);
      float d11 = x1y1.dot(v11);
      
      /// Interpolate dot product values at sample point using polynomial interpolation 6x^5 - 15x^4 + 10x^3
      float yf = Y - Y0; // Float offset inside the square [0, 1]
      float xf = X - X0; // Float offset inside the square [0, 1]
      
      float wx = smoothstep(xf);
      float wy = smoothstep(yf);
      
      /// Interpolate along x for the contributions from each of the gradients
      float xa = lerp(wx, d00, d10);
      float xb = lerp(wx, d01, d11);
      
      float val = lerp(wy, xa, xb);
      
      return clamp(val, -1.0f, 1.0);
    }
    
    float get_value(float X, float Y, float Z) const override {
      /// Compress the coordinates inside the chunk; float part + int part = point coordinate
      /// Grid points from the chunk in the world
      int X0 = (int) std::floor(X);
      int Y0 = (int) std::floor(Y);
      int X1 = (int) std::ceil(X);
      int Y1 = (int) std::ceil(Y);
      int Z0 = (int) std::floor(Z);
      int Z1 = (int) std::ceil(Z);
      
      /// Gradients using hashed indices from lookup list
      Vec3f x0y0z0 = grads3[perms[(X0 + perms[(Y0 + perms[Z0 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x1y0z0 = grads3[perms[(X1 + perms[(Y0 + perms[Z0 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x0y1z0 = grads3[perms[(X0 + perms[(Y1 + perms[Z0 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x1y1z0 = grads3[perms[(X1 + perms[(Y1 + perms[Z0 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      
      Vec3f x0y0z1 = grads3[perms[(X0 + perms[(Y0 + perms[Z1 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x1y0z1 = grads3[perms[(X1 + perms[(Y0 + perms[Z1 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x0y1z1 = grads3[perms[(X0 + perms[(Y1 + perms[Z1 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      Vec3f x1y1z1 = grads3[perms[(X1 + perms[(Y1 + perms[Z1 % perms.size()]) % perms.size()]) %
                                         perms.size()]];
      
      /// Vectors from gradients to point in unit cube
      Vec3f v000 = Vec3f{X - X0, Y - Y0, Z - Z0};
      Vec3f v100 = Vec3f{X - X1, Y - Y0, Z - Z0};
      Vec3f v010 = Vec3f{X - X0, Y - Y1, Z - Z0};
      Vec3f v110 = Vec3f{X - X1, Y - Y1, Z - Z0};
      
      Vec3f v001 = Vec3f{X - X0, Y - Y0, Z - Z1};
      Vec3f v101 = Vec3f{X - X1, Y - Y0, Z - Z1};
      Vec3f v011 = Vec3f{X - X0, Y - Y1, Z - Z1};
      Vec3f v111 = Vec3f{X - X1, Y - Y1, Z - Z1};
      
      /// Contribution of gradient vectors by dot product between relative vectors and gradients
      float d000 = dot(x0y0z0, v000);
      float d100 = dot(x1y0z0, v100);
      float d010 = dot(x0y1z0, v010);
      float d110 = dot(x1y1z0, v110);
      
      float d001 = dot(x0y0z1, v001);
      float d101 = dot(x1y0z1, v101);
      float d011 = dot(x0y1z1, v011);
      float d111 = dot(x1y1z1, v111);
      
      /// Interpolate dot product values at sample point using polynomial interpolation 6x^5 - 15x^4 + 10x^3
      float yf = Y - Y0; // Float offset inside the cube [0, 1]
      float xf = X - X0; // Float offset inside the cube [0, 1]
      float zf = Z - Z0; // Float offset inside the cube [0, 1]
      
      float wx = smoothstep(xf);
      float wy = smoothstep(yf);
      float wz = smoothstep(zf);
      
      /// Interpolate along x for the contributions from each of the gradients
      float xa = lerp(wx, d000, d100);
      float xb = lerp(wx, d010, d110);
      
      float xc = lerp(wx, d001, d101);
      float xd = lerp(wx, d011, d111);
      
      /// Interpolate along y for the contributions from each of the gradients
      float ya = lerp(wy, xa, xb);
      float yb = lerp(wy, xc, xd);
      
      /// Interpolate along z for the contributions from each of the gradients
      float za = lerp(wz, ya, yb);
      
      return clamp(za, -1.0f, 1.0f);
    }
  };
}

#endif // NOISE_H
