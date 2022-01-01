#include <cstddef>
#include <cstring>
#include <vector>
#include <cstdio>
#include "triangulation2d.hpp"

namespace polylower {

  inline double cross2d(const Point &p1, const Point &p2) {
    return (p1.x * p2.y) - (p1.y * p2.x);
  }

  inline Point sub(const Point &p1, const Point &p2) {
    return (Point){p1.x - p2.x, p1.y - p2.y};
  }
  
  bool point_in_triangle(const Point &p, const Point &p1, const Point &p2, const Point &p3) {
    Point v1 = (Point){p1.x - p.x, p1.y - p.y};
    Point v2 = (Point){p2.x - p.x, p2.y - p.y};
    Point v3 = (Point){p3.x - p.x, p3.y - p.y};
    double c1 = cross2d(v1, v2);
    double c2 = cross2d(v2, v3);
    double c3 = cross2d(v3, v1);
    return (c1 > 0 && c2 > 0 && c3 > 0) || (c1 < 0 && c2 < 0 && c3 < 0);
  }

  bool is_ear(int x, Point *p, int *prev, int *next) {
    if (cross2d(sub(p[x], p[prev[x]]), sub(p[next[x]], p[x])) <= 0) {
      // printf("spike %d\n", x);
      return false;
    }
    for (int i = next[next[x]]; i != prev[x]; i = next[i]) {
      if (point_in_triangle(p[i], p[prev[x]], p[x], p[next[x]])) {
        return false;
      }
    }
    return true;
  }

  /* return a triangle list vector<int> = {t1a, t1b, t1c, t2a, t2b, t2c, ...} */
  std::vector<int> triangulation2d(Point *plist, size_t p_num) {
    int *prev = new int[p_num];
    int *next = new int[p_num];
    int *ear  = new int[p_num];
    for (int i = 0; i < p_num; i++) {
      prev[i] = i - 1;
      next[i] = i + 1;
    }
    prev[0] = p_num-1;
    next[p_num-1] = 0;

    memset(ear, -1, p_num * sizeof(int));

    std::vector<int> tv(3 * (p_num - 2));

    int i = 0;
    int k;
    for (k = 0; k < p_num-3; k++) {
      while (true) {
        if (ear[i] == -1) ear[i] = is_ear(i, plist, prev, next);
        if (ear[i] == 1) break;
        i = next[i];
      }

      // printf("found ear %d\n", i);
      // make prev[i], i, next[i] a triangle
      tv[k*3] = prev[i];
      tv[k*3+1] = i;
      tv[k*3+2] = next[i];
      
      ear[prev[i]] = ear[next[i]] = -1;
      next[prev[i]] = next[i];
      prev[next[i]] = prev[i];
      i = next[i];
    }
    tv[k * 3] = prev[i];
    tv[k * 3 + 1] = i;
    tv[k * 3 + 2] = next[i];
    delete [] prev;
    delete [] next;
    delete [] ear;
    return tv;
  }
  
};
