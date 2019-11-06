#!/bin/bash

# rm -rf build
mkdir build
cd build

cmake ..
make -j 8
cd ..
# python3 main.py