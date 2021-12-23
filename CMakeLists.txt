cmake_minimum_required(VERSION 3.4)
project(polylower)

add_subdirectory(pybind11)
pybind11_add_module(polylower src/polylower.cpp src/mesh.cpp src/triangulation2d.cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS False)

add_compile_options(-Wall)


enable_testing()
add_test(NAME pytest COMMAND test/test.sh)
