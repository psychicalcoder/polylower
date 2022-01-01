#ifndef POLYLOWER_UTIL_H
#define POLYLOWER_UTIL_H

#include "vec.hpp"
#include <vector>
#include <array>

using std::vector;
using std::array;

struct Mat4 {
  double *m;
  Mat4() : m(new double[16]) {
    for (int i = 0; i < 16; i++)
      m[i] = 0.0;
  }
  Mat4(const Mat4 &M) {
    for (int i = 0; i < 16; i++)
      m[i] = M.m[i];
  }
  
  Mat4(Mat4 &&M) : m(M.m) {
    M.m = nullptr;
  }

  Mat4& operator=(const Mat4 &rhs) {
    for (int i = 0; i < 16; i++)
      m[i] = rhs.m[i];
    return *this;
  }

  Mat4& operator=(Mat4 &&rhs) {
    if (this != &rhs) {
      delete [] m;
      m = rhs.m;
      rhs.m = nullptr;
    }
    return *this;
  }
  
  Mat4 operator+(const Mat4 &rhs) const {
    Mat4 n;
    for (int i = 0; i < 16; i++) {
      m[i] = rhs.m[i];
    }
    return n;
  }

  Mat4& operator+=(const Mat4 &rhs) {
    for (int i = 0; i < 16; i++)
      m[i] += rhs.m[i];
    return *this;
  }

  int solve(vec &ans);
};

void triangularize(const vector<vec> &vertices, const vector<vector<int>> &faces, vector<array<int, 3>> &tfaces);

#endif
