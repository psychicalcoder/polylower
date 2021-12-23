#!/bin/bash

set -x

rm -rf build
mkdir build
cd build
cmake ../
cmake --build ./
python3 -m pytest ../test/test.py
