cmake_minimum_required(VERSION 3.4)
project(polylower)

add_subdirectory(pybind11)
pybind11_add_module(polylower src/polylower.cpp)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_EXPORT_COMPILE_COMMANDS True)
