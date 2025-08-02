#!/bin/bash

# Options for TPX3 simulation compilation
multiThread=0

# Process input options/parameters
while getopts "pmqnadt:j:e:cli" arg
do
    case $arg in
        m)
            # Run multi-threaded simulation
            multiThread=1
        ;;
        ?)
            echo "Unknown argument $arg"
        ;;
    esac
done

rm -rf build
mkdir build
cd build || exit
if [ $multiThread -eq 1 ]; then
    cmake -DCMAKE_CXX_STANDARD="17" -DUSE_MT=ON  ..
else
    cmake -DCMAKE_CXX_STANDARD="17" ..
fi
make -j8
