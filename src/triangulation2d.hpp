#ifndef _POLYLOWER_TRI2D_H
#define _POLYLOWER_TRI2D_H
#include <cstddef>
#include <vector>

namespace polylower {
  struct Point { double x, y; };

    /* return a triangle list vector<int> = {t1a, t1b, t1c, t2a, t2b, t2c, ...} */
  std::vector<int> triangulation2d(Point *plist, size_t p_num);
}

#endif
