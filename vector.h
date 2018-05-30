#ifndef RAYTRACER_VECTOR_H
#define RAYTRACER_VECTOR_H

#include <complex>

template<typename T = double>
struct Vec3 {
  T x, y, z;
  
  constexpr Vec3(): x(0.0), y(0.0), z(0.0) {};
  constexpr explicit Vec3(T value): x(value), y(value), z(value) {};
  constexpr Vec3(T x, T y, T z): x(x), y(y), z(z) {};
  
  /// Operators
  inline Vec3<T> operator-() const { return Vec3<T>{-x, -y, -z}; }
  inline T operator[](size_t i) const {
    switch (i) { // Hopefully optimised into a jump table
      case 0:
        return x;
      case 1:
        return y;
      case 2:
        return z;
      default:
        return T{};
    }
  };
  
  inline double length() const { return std::sqrt(x*x + y*y + z*z); }
  inline double squared_length() const { return x*x + y*y + z*z; }
  inline void normalize() { auto lng = length(); x /= lng; y /= lng; z /= lng; }
  inline Vec3<T> normalized() const { auto lng = length(); return {x / lng, y / lng, z / lng}; }
};

/// Vec3 operators
template<typename T>
inline Vec3<T> operator*(const Vec3<T> &v, double r) { return {v.x * r, v.y * r, v.z * r}; }

template<typename T>
inline Vec3<T> operator*(double l, const Vec3<T>& v) { return {v.x * l, v.y * l, v.z * l}; }

template<typename T>
inline Vec3<T> operator*(const Vec3<T>& l, const Vec3<T>& r) { return {l.x * r.x, l.y * r.y, l.z * r.z}; }

template<typename T>
inline Vec3<T> operator+(const Vec3<T>& l, const Vec3<T>& r) { return {l.x + r.x, l.y + r.y, l.z + r.z}; }

template<typename T>
inline Vec3<T> operator+(const Vec3<T>& l, T r) { return {l.x + r, l.y + r, l.z + r}; }

template<typename T>
inline void operator+=(Vec3<T>& l, const Vec3<T>& r) { l.x += r.x; l.y += r.y; l.z += r.z; }

template<typename T>
inline Vec3<T> operator-(const Vec3<T>& l, const Vec3<T>& r) { return {l.x - r.x, l.y - r.y, l.z - r.z}; }

template<typename T>
inline Vec3<T> operator/(const Vec3<T>& l, T r) { return {l.x / r, l.y / r, l.z / r}; }

template<typename T>
inline void operator/=(Vec3<T>& l, const Vec3<T>& r) { l.x /= r.x; l.y /= r.y; l.z /= r.z; }

template<typename T>
inline void operator/=(Vec3<T>& l, T r) { l.x /= r; l.y /= r; l.z /= r; }

template<typename T>
inline double dot(const Vec3<T>& l, const Vec3<T>& r) { return l.x * r.x + l.y * r.y + l.z * r.z; }

template<typename T>
inline Vec3<T> cross(const Vec3<T>& l, const Vec3<T>& r) {
  return Vec3<T>{l.y*r.z - l.z*r.y, -(l.x*r.z - l.z*r.x), l.x*r.y - l.y*r.x};
}

/// Linear Algebra
inline Vec3<> reflect(const Vec3<>& v, const Vec3<>& n) {
  return v - 2*dot(v, n)*n;
}

Vec3<> random_in_unit_sphere() {
  Vec3<> p;
  do {
    p = 2.0 * Vec3<>{drand48(), drand48(), drand48()} - Vec3<>{1.0};
  } while (p.squared_length() >= 1.0);
  return p;
}

/// p(t) = A + t * B
struct Ray {
  // Origin
  Vec3<> A;
  // Direction
  Vec3<> B;
  
  Ray() = default;
  Ray(Vec3<> a, Vec3<> b): A(a), B(b) {};
  
  /// Operators
  inline Vec3<> operator()(double t) const { return A + t * B; }
  
  inline const Vec3<>& origin() const { return A; };
  inline const Vec3<>& direction() const { return B; };
  
};

#endif // RAYTRACER_VECTOR_H
