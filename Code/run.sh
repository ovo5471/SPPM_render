#!/usr/bin/env bash
ITERATIONS=5000
PHOTONS=1000000
CKPT_INTERVAL=1
# If project not ready, generate cmake file.
if [[ ! -d build ]]; then
    mkdir -p build
    cd build
    cmake ..
    cd ..
fi

# Build project.
cd build
make -j
cd ..

# Run all testcases. 
# You can comment some lines to disable the run of specific examples.
mkdir -p output
mkdir -p output/scene_test
time bin/PA1 testcases/scene_test.txt output/scene_test $ITERATIONS $PHOTONS $CKPT_INTERVAL
