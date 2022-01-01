cmake_minimum_required(VERSION 3.4)
project(polylower)

add_subdirectory(pybind11)
pybind11_add_module(polylower ./mesh.cpp ./util.cpp ./triangulation2d.cpp ./vec.hpp)
