#ifndef RAYTRACER_VECTOR_H
#define RAYTRACER_VECTOR_H

#include <complex>

template<typename T = double>
struct Vec3 {
  T x, y, z;
  
  Vec3(): x(0.0f), y(0.0f), z(0.0f) {};
  explicit Vec3(T value): x(value), y(value), z(value) {};
  Vec3(T x, T y, T z): x(x), y(y), z(z) {};
  
  /// Operators
  inline Vec3<T> operator-() const { return Vec3<T>{-x, -y, -z}; }
  
  inline double length() const { return std::sqrt(x*x + y*y + z*z); }
  inline void normalize() { auto lng = length(); x /= lng; y /= lng; z /= lng;  }
  inline Vec3<T> normalized() { auto lng = length(); return {x /= lng, y /= lng, z /= lng};  }
};

// Vec3 operators
template<typename T>
inline Vec3<T> operator*(const Vec3<T> &v, double r) { return {v.x * r, v.y * r, v.z * r}; }
template<typename T>
inline Vec3<T> operator*(double l, const Vec3<T>& v) { return {v.x * l, v.y * l, v.z * l}; }
template<typename T>
inline Vec3<T> operator+(const Vec3<T>& l, const Vec3<T>& r) { return {l.x + r.x, l.y + r.y, l.z + r.z}; }
template<typename T>
inline Vec3<T> operator-(const Vec3<T>& l, const Vec3<T>& r) { return {l.x - r.x, l.y - r.y, l.z - r.z}; }

/// p(t) = A + t * B
struct Ray {
  // Origin
  Vec3<> A;
  // Direction
  Vec3<> B;
  
  Ray() = default;
  Ray(Vec3<> a, Vec3<> b): A(a), B(b) {};
  
  /// Operators
  Vec3<> operator()(double t) const { return A + t * B; }
  
  Vec3<> origin() const { return A; };
  Vec3<> direction() const { return B; };
  
};

#endif // RAYTRACER_VECTOR_H
