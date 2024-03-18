#!/bin/bash

# Serial simulations for Newton's Second Law
echo "Starting serial simulation for Newton's Second Law: Sun-Earth system"
./nbody-s 60 31557600 1000 examples/sun-earth.npy output/serial/2ndLaw/sun-earthActual.npy

echo "Starting serial simulation for Newton's Second Law: Solar System Inners"
./nbody-s 60 31557600 1000 examples/solar-system-inner.npy output/serial/2ndLaw/solar-system-innersActual.npy

echo "Starting serial simulation for Newton's Second Law: Solar System"
./nbody-s 60 31557600 1000 examples/solar-system.npy output/serial/2ndLaw/solar-system.npy

echo "Starting serial simulation for Newton's Second Law: Pluto-Charon system"
./nbody-s 1 302400 1000 examples/pluto-charon.npy output/serial/2ndLaw/pluto-charonActual.npy

echo "Starting serial simulation for Newton's Second Law: Figure 8"
./nbody-s 1e-4 2.1 1000 examples/figure8.npy output/serial/2ndLaw/figure8Actual.npy

echo "Starting serial simulation for Newton's Second Law: Figure 8 Rotate"
./nbody-s 1e-4 75 1000 examples/figure8-rotate.npy output/serial/2ndLaw/figure8-rotateActual.npy

echo "Starting serial simulation for Newton's Second Law: Random 25"
./nbody-s 0.1 500000 1000 examples/random25.npy output/serial/2ndLaw/random25Actual.npy

echo "Starting serial simulation for Newton's Second Law: Random 100"
./nbody-s 0.01 1000 1000 examples/random100.npy output/serial/2ndLaw/random100Actual.npy

echo "Starting serial simulation for Newton's Second Law: Random 1000"
./nbody-s 0.01 1000 1000 examples/random1000.npy output/serial/2ndLaw/random1000Actual.npy

# echo "Starting serial simulation for Newton's Second Law: Random 10000"
# ./nbody-s 0.01 1000 1000 examples/random10000.npy output/serial/2ndLaw/random10000Actual.npy

# Serial simulations for Newton's Third Law
echo "Starting serial simulation for Newton's Third Law: Sun-Earth system"
./nbody-s3 60 31557600 1000 examples/sun-earth.npy output/serial/3rdLaw/sun-earthActual.npy

echo "Starting serial simulation for Newton's Third Law: Solar System Inners"
./nbody-s3 60 31557600 1000 examples/solar-system-inners.npy output/serial/3rdLaw/solar-system-innersActual.npy

echo "Starting serial simulation for Newton's Third Law: Solar System"
./nbody-s3 60 31557600 1000 examples/solar-system.npy output/serial/3rdLaw/solar-system.npy

echo "Starting serial simulation for Newton's Third Law: Pluto-Charon system"
./nbody-s3 1 302400 1000 examples/pluto-charon.npy output/serial/3rdLaw/pluto-charonActual.npy

echo "Starting serial simulation for Newton's Third Law: Figure 8"
./nbody-s3 1e-4 2.1 1000 examples/figure8.npy output/serial/3rdLaw/figure8Actual.npy

echo "Starting serial simulation for Newton's Third Law: Figure 8 Rotate"
./nbody-s3 1e-4 75 1000 examples/figure8-rotate.npy output/serial/3rdLaw/figure8-rotateActual.npy

echo "Starting serial simulation for Newton's Third Law: Random 25"
./nbody-s3 0.1 500000 1000 examples/random25.npy output/serial/3rdLaw/random25Actual.npy

echo "Starting serial simulation for Newton's Third Law: Random 100"
./nbody-s3 0.01 1000 1000 examples/random100.npy output/serial/3rdLaw/random100Actual.npy

echo "Starting serial simulation for Newton's Third Law: Random 1000"
./nbody-s3 0.01 1000 1000 examples/random1000.npy output/serial/3rdLaw/random1000Actual.npy

# echo "Starting serial simulation for Newton's Third Law: Random 10000"
# ./nbody-s3 0.01 1000 1000 examples/random10000.npy output/serial/3rdLaw/random10000Actual.npy

# Parallel simulations for Newton's Second Law
echo "Starting parallel simulation for Newton's Second Law: Sun-Earth system"
./nbody-p 60 31557600 1000 examples/sun-earth.npy output/parallel/2ndLaw/sun-earthActual.npy 8

echo "Starting parallel simulation for Newton's Second Law: Solar System Inners"
./nbody-p 60 31557600 1000 examples/solar-system-inners.npy output/parallel/2ndLaw/solar-system-innersActual.npy 8

echo "Starting parallel simulation for Newton's Second Law: Solar System"
./nbody-p 60 31557600 1000 examples/solar-system.npy output/parallel/2ndLaw/solar-system.npy 8

echo "Starting parallel simulation for Newton's Second Law: Pluto-Charon system"
./nbody-p 1 302400 1000 examples/pluto-charon.npy output/parallel/2ndLaw/pluto-charonActual.npy 8

echo "Starting parallel simulation for Newton's Second Law: Figure 8"
./nbody-p 1e-4 2.1 1000 examples/figure8.npy output/parallel/2ndLaw/figure8Actual.npy 8

echo "Starting parallel simulation for Newton's Second Law: Figure 8 Rotate"
./nbody-p 1e-4 75 1000 examples/figure8-rotate.npy output/parallel/2ndLaw/figure8-rotateActual.npy 8

echo "Starting parallel simulation for Newton's Second Law: Random 25"
./nbody-p 0.1 500000 1000 examples/random25.npy output/parallel/2ndLaw/random25Actual.npy 8

echo "Starting parallel simulation for Newton's Second Law: Random 100"
./nbody-p 0.01 1000 1000 examples/random100.npy output/parallel/2ndLaw/random100Actual.npy 8

echo "Starting parallel simulation for Newton's Second Law: Random 1000"
./nbody-p 0.01 1000 1000 examples/random1000.npy output/parallel/2ndLaw/random1000Actual.npy 8

echo "Starting parallel simulation for Newton's Second Law: Random 10000"
./nbody-p 0.01 1000 1000 examples/random10000.npy output/parallel/2ndLaw/random10000Actual.npy 8
