#ifndef RAYTRACER_VECTOR_H
#define RAYTRACER_VECTOR_H

#include <complex>

template<typename T>
struct Vec2 {
  T x, y;
  Vec2(T x, T y): x(x), y(y) {};
  Vec2(): x(0.0f), y(0.0f) {};
  
  /// Sum of the components of the vector
  inline T sum() const { return x + y; }
  
  /// Floors the components and returns a copy
  inline Vec2<T> floor() const { return {std::floor(x), std::floor(y)}; }

  /// Rounds the components up and returns a copy
  inline Vec2<T> ceil() const { return {std::ceil(x), std::ceil(y)}; }
  
  /// Dot product
  inline T dot(Vec2<T> u) const { return x * u.x + y * u.y; }
  
  /// Operators
  Vec2<T> operator+(const Vec2 &rhs) const { return {x + rhs.x, y + rhs.y}; }
  
  Vec2<T> operator-(const Vec2 &rhs) const { return {x - rhs.x, y - rhs.y}; }
  
  bool operator==(const Vec2 &rhs) const { return x == rhs.x && y == rhs.y; }
  
  /// Returns a copy of this vector normalized
  inline Vec2<T> normalized() const { float lng = length(); return {x / lng, y / lng}; }
  
  /// Normalisation in place
  inline void normalize() { float lng = length(); x /= lng; y /= lng; }
  
  /// Length of the vector
  inline float length() const { return std::sqrt(std::pow(x, 2.0f) + std::pow(y, 2.0f)); }
  
  friend std::ostream &operator<<(std::ostream& os, const Vec2<T> &v) {
    return os << "(x:" << v.x << " y:" << v.y << ")";
  }
};

template<typename T>
struct Vec3 {
  T x, y, z;
  
  constexpr Vec3(): x(0.0f), y(0.0f), z(0.0f) {};
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
  
  /// Element wide operators
  inline bool operator==(T rhs) const { return x == rhs && y == rhs && z == rhs; }
  inline bool operator<=(T rhs) const { return x <= rhs && y <= rhs && z <= rhs; }
  
  Vec3<T> floor() const { return {std::floor(x), std::floor(y), std::floor(z)}; }
  inline float sum() const { return x + y + z; }
  inline float length() const { return std::sqrt(x*x + y*y + z*z); }
  inline float squared_length() const { return x*x + y*y + z*z; }
  inline void normalize() { auto lng = length(); x /= lng; y /= lng; z /= lng; }
  inline Vec3<T> normalized() const { auto lng = length(); return {x / lng, y / lng, z / lng}; }
  
  friend std::ostream &operator<<(std::ostream& os, const Vec3<T>& v) {
    return os << "(x:" << v.x << " y:" << v.y << " z:" << v.z << ")";
  }
};

/// Type helpers
using Vec2f = Vec2<float>;
using Vec3f = Vec3<float>;
using Vec3d = Vec3<float>;

/// Vec3 operators
template<typename T>
inline Vec3<T> operator*(const Vec3<T> &v, float r) { return {v.x * r, v.y * r, v.z * r}; }

template<typename T>
inline Vec3<T> operator*(float l, const Vec3<T>& v) { return {v.x * l, v.y * l, v.z * l}; }

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
inline float dot(const Vec3<T>& l, const Vec3<T>& r) { return l.x * r.x + l.y * r.y + l.z * r.z; }

template<typename T>
inline Vec3<T> cross(const Vec3<T>& l, const Vec3<T>& r) {
  return Vec3<T>{l.y*r.z - l.z*r.y, -(l.x*r.z - l.z*r.x), l.x*r.y - l.y*r.x};
}

/// Linear Algebra
inline Vec3f reflect(const Vec3f& v, const Vec3f& n) {
  return v - 2*dot(v, n)*n;
}

/// p(t) = A + t * B
struct Ray {
  // Origin
  Vec3f A;
  // Direction
  Vec3f B;
  
  Ray() = default;
  Ray(Vec3f a, Vec3f b): A(a), B(b) {};
  
  /// Operators
  inline Vec3f operator()(float t) const { return A + t * B; }
 
  inline const Vec3f& origin() const { return A; };
  inline const Vec3f& direction() const { return B; };
};

#endif // RAYTRACER_VECTOR_H
