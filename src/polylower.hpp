#ifndef _POLYLOWER_H
#define _POLYLOWER_H

#include "mesh.hpp"

namespace polylower {

  Mesh simplify(const Mesh &mesh, double ratio);
  
}

#endif
