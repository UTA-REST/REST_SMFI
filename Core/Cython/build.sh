#!/usr/bin/env bash

# build the code
python3 setup_cython.py clean

python3 setup_cython.py build_ext --inplace


