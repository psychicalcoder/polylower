#ifndef _POLYLOWER_OBJ_H
#define _POLYLOWER_OBJ_H
#include "mesh.hpp"

namespace polylower {

  Mesh load_mesh(const char *path);
  void save_mesh(const char *path, const Mesh &mesh);

}
#endif
