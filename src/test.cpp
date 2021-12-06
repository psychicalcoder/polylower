#include "obj.hpp"
#include <iostream>

using namespace polylower;

int main() {
  Mesh mesh = load_mesh("bunny.obj");
  std::cout << "mesh loaded" << std::endl;
  save_mesh("mod_bunny.obj", mesh);
  std::cout << "mesh saved" << std::endl;
  return 0;
}
