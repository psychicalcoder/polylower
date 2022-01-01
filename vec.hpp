#ifndef POLYLOWER_VEC_H
#define POLYLOWER_VEC_H

#include <utility>
#include <cmath>

#define PL_EPS 1e-7

struct vec {
  double x, y, z;
  vec() : x(0.0), y(0.0), z(0.0) {}
  ~vec() {}
  vec(double x, double y, double z) : x(x), y(y), z(z) {}
  vec(const vec &v) : x(v.x), y(v.y), z(v.z) {}
  vec(vec &&v) : x(v.x), y(v.y), z(v.z) {}
  
  vec& operator=(const vec &rhs) {
    x = rhs.x; y = rhs.y; z = rhs.z;
    return *this;
  }

  vec& operator=(vec &&rhs) {
    x = rhs.x; y = rhs.y; z = rhs.z;
    return *this;
  }

  bool operator==(const vec &rhs) const {
    return
      fabs(x-rhs.x) < PL_EPS
      && fabs(y-rhs.y) < PL_EPS
      && fabs(z-rhs.z) < PL_EPS;
  }

  vec operator+(const vec &rhs) const {
    return vec(x+rhs.x, y+rhs.y, z+rhs.z);
  }

  vec& operator+=(const vec &rhs) {
    x += rhs.x; y += rhs.y; z += rhs.z;
    return *this;
  }

  vec operator-(const vec &rhs) const {
    return vec(x-rhs.x, y-rhs.y, z-rhs.z);
  }

  vec& operator-=(const vec &rhs) {
    x -= rhs.x; y -= rhs.y; z -= rhs.z;
    return *this;
  }

  /*
   * Cross Product
   */
  vec operator*(const vec &rhs) const {
    return vec(y * rhs.z - z * rhs.y,
               z * rhs.x - x * rhs.z,
               x * rhs.y - y * rhs.x);
  }

  /*
   * Dot Product
   */
  double operator^(const vec &rhs) const {
    return x * rhs.x + y * rhs.y + z * rhs.z;
  }

  vec cross(const vec &rhs) const {
    return (*this)*rhs;
  }

  double dot(const vec &rhs) const {
    return (*this)^rhs;
  }

  double len() const {
    return sqrt(x * x + y * y + z * z);
  }

  vec unit() const {
    double l = this->len();
    return vec(x / l, y / l, z / l);
  }

  vec &make_unit() {
    double l = this->len();
    x /= l; y /= l; z /= l;
    return *this;
  }
};

inline vec unit(const vec& v) {
  double l = v.len();
  return vec(v.x / l, v.y / l, v.z / l);
}

#endif
