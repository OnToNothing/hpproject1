#!/bin/bash

# Define log files for each simulation type
LOG_FILE_SERIAL_2ND="Log/serial_2ndLaw.log"
LOG_FILE_SERIAL_3RD="Log/serial_3rdLaw.log"
LOG_FILE_PARALLEL_2ND="Log/parallel_2ndLaw.log"
LOG_FILE_PARALLEL_3RD="Log/parallel_3rdLaw.log"

# Function to run simulation and comparison
run_simulation_and_compare() {
    simulation_type=$1  # Command to run, e.g., "./nbody-s"
    log_file=$2         # Log file to write to
    simulation_args=$3  # Arguments for the simulation command
    expected_result_path=$4 # Path to the expected result .npy file
    actual_result_path=$5   # Path to the actual result .npy file

    echo "Starting $simulation_type simulation with args: $simulation_args" | tee -a $log_file
    $simulation_type $simulation_args 2>&1 | tee -a $log_file
    python scripts/compare_npy.py $expected_result_path $actual_result_path 2>&1 | tee -a $log_file
    echo "" >> $log_file
}

# Clear existing log files
> $LOG_FILE_SERIAL_2ND
> $LOG_FILE_SERIAL_3RD
> $LOG_FILE_PARALLEL_2ND
> $LOG_FILE_PARALLEL_3RD

# Serial simulations for Newton's Second Law
run_simulation_and_compare "./nbody-s" $LOG_FILE_SERIAL_2ND "60 31557600 1000 examples/sun-earth.npy output/serial/2ndLaw/sun-earthActual.npy" "examples/sun-earth-expected.npy" "output/serial/2ndLaw/sun-earthActual.npy"
run_simulation_and_compare "./nbody-s" $LOG_FILE_SERIAL_2ND "60 31557600 1000 examples/solar-system-inner.npy output/serial/2ndLaw/solar-system-innerActual.npy" "examples/solar-system-inner-expected.npy" "output/serial/2ndLaw/solar-system-innerActual.npy"
run_simulation_and_compare "./nbody-s" $LOG_FILE_SERIAL_2ND "60 31557600 1000 examples/solar-system.npy output/serial/2ndLaw/solar-systemActual.npy" "examples/solar-system-expected.npy" "output/serial/2ndLaw/solar-systemActual.npy"
run_simulation_and_compare "./nbody-s" $LOG_FILE_SERIAL_2ND "1 302400 1000 examples/pluto-charon.npy output/serial/2ndLaw/pluto-charonActual.npy" "examples/pluto-charon-expected.npy" "output/serial/2ndLaw/pluto-charonActual.npy"
run_simulation_and_compare "./nbody-s" $LOG_FILE_SERIAL_2ND "1e-4 2.1 1000 examples/figure8.npy output/serial/2ndLaw/figure8Actual.npy" "examples/figure8-expected.npy" "output/serial/2ndLaw/figure8Actual.npy"
run_simulation_and_compare "./nbody-s" $LOG_FILE_SERIAL_2ND "1e-4 75 1000 examples/figure8-rotate.npy output/serial/2ndLaw/figure8-rotateActual.npy" "examples/figure8-rotate-expected.npy" "output/serial/2ndLaw/figure8-rotateActual.npy"
run_simulation_and_compare "./nbody-s" $LOG_FILE_SERIAL_2ND "0.1 500000 1000 examples/random25.npy output/serial/2ndLaw/random25Actual.npy" "examples/random25-expected.npy" "output/serial/2ndLaw/random25Actual.npy"
run_simulation_and_compare "./nbody-s" $LOG_FILE_SERIAL_2ND "0.01 1000 1000 examples/random100.npy output/serial/2ndLaw/random100Actual.npy" "examples/random100-expected.npy" "output/serial/2ndLaw/random100Actual.npy"

# Serial simulations for Newton's Third Law using nbody-s3
run_simulation_and_compare "./nbody-s3" $LOG_FILE_SERIAL_3RD "60 31557600 1000 examples/sun-earth.npy output/serial/3rdLaw/sun-earthActual.npy" "examples/sun-earth-expected.npy" "output/serial/3rdLaw/sun-earthActual.npy"
run_simulation_and_compare "./nbody-s3" $LOG_FILE_SERIAL_3RD "60 31557600 1000 examples/solar-system-inner.npy output/serial/3rdLaw/solar-system-innerActual.npy" "examples/solar-system-inner-expected.npy" "output/serial/3rdLaw/solar-system-innerActual.npy"
run_simulation_and_compare "./nbody-s3" $LOG_FILE_SERIAL_3RD "60 31557600 1000 examples/solar-system.npy output/serial/3rdLaw/solar-systemActual.npy" "examples/solar-system-expected.npy" "output/serial/3rdLaw/solar-systemActual.npy"
run_simulation_and_compare "./nbody-s3" $LOG_FILE_SERIAL_3RD "1 302400 1000 examples/pluto-charon.npy output/serial/3rdLaw/pluto-charonActual.npy" "examples/pluto-charon-expected.npy" "output/serial/3rdLaw/pluto-charonActual.npy"
run_simulation_and_compare "./nbody-s3" $LOG_FILE_SERIAL_3RD "1e-4 2.1 1000 examples/figure8.npy output/serial/3rdLaw/figure8Actual.npy" "examples/figure8-expected.npy" "output/serial/3rdLaw/figure8Actual.npy"
run_simulation_and_compare "./nbody-s3" $LOG_FILE_SERIAL_3RD "1e-4 75 1000 examples/figure8-rotate.npy output/serial/3rdLaw/figure8-rotateActual.npy" "examples/figure8-rotate-expected.npy" "output/serial/3rdLaw/figure8-rotateActual.npy"
run_simulation_and_compare "./nbody-s3" $LOG_FILE_SERIAL_3RD "0.1 500000 1000 examples/random25.npy output/serial/3rdLaw/random25Actual.npy" "examples/random25-expected.npy" "output/serial/3rdLaw/random25Actual.npy"
run_simulation_and_compare "./nbody-s3" $LOG_FILE_SERIAL_3RD "0.01 1000 1000 examples/random100.npy output/serial/3rdLaw/random100Actual.npy" "examples/random100-expected.npy" "output/serial/3rdLaw/random100Actual.npy"

# Parallel simulations for Newton's Second Law
run_simulation_and_compare "./nbody-p" $LOG_FILE_PARALLEL_2ND "60 31557600 1000 examples/sun-earth.npy output/parallel/2ndLaw/sun-earthActual.npy 8" "examples/sun-earth-expected.npy" "output/parallel/2ndLaw/sun-earthActual.npy"
run_simulation_and_compare "./nbody-p" $LOG_FILE_PARALLEL_2ND "60 31557600 1000 examples/solar-system-inner.npy output/parallel/2ndLaw/solar-system-innerActual.npy 8" "examples/solar-system-inner-expected.npy" "output/parallel/2ndLaw/solar-system-innerActual.npy"
run_simulation_and_compare "./nbody-p" $LOG_FILE_PARALLEL_2ND "60 31557600 1000 examples/solar-system.npy output/parallel/2ndLaw/solar-systemActual.npy 8" "examples/solar-system-expected.npy" "output/parallel/2ndLaw/solar-systemActual.npy"
run_simulation_and_compare "./nbody-p" $LOG_FILE_PARALLEL_2ND "1 302400 1000 examples/pluto-charon.npy output/parallel/2ndLaw/pluto-charonActual.npy 8" "examples/pluto-charon-expected.npy" "output/parallel/2ndLaw/pluto-charonActual.npy"
run_simulation_and_compare "./nbody-p" $LOG_FILE_PARALLEL_2ND "1e-4 2.1 1000 examples/figure8.npy output/parallel/2ndLaw/figure8Actual.npy 8" "examples/figure8-expected.npy" "output/parallel/2ndLaw/figure8Actual.npy"
run_simulation_and_compare "./nbody-p" $LOG_FILE_PARALLEL_2ND "1e-4 75 1000 examples/figure8-rotate.npy output/parallel/2ndLaw/figure8-rotateActual.npy 8" "examples/figure8-rotate-expected.npy" "output/parallel/2ndLaw/figure8-rotateActual.npy"
run_simulation_and_compare "./nbody-p" $LOG_FILE_PARALLEL_2ND "0.1 500000 1000 examples/random25.npy output/parallel/2ndLaw/random25Actual.npy 8" "examples/random25-expected.npy" "output/parallel/2ndLaw/random25Actual.npy"
run_simulation_and_compare "./nbody-p" $LOG_FILE_PARALLEL_2ND "0.01 1000 1000 examples/random100.npy output/parallel/2ndLaw/random100Actual.npy 8" "examples/random100-expected.npy" "output/parallel/2ndLaw/random100Actual.npy"

# Parallel simulations for Newton's Third Law using nbody-p3
run_simulation_and_compare "./nbody-p3" $LOG_FILE_PARALLEL_3RD "60 31557600 1000 examples/sun-earth.npy output/parallel/3rdLaw/sun-earthActual.npy 8" "examples/sun-earth-expected.npy" "output/parallel/3rdLaw/sun-earthActual.npy"
run_simulation_and_compare "./nbody-p3" $LOG_FILE_PARALLEL_3RD "60 31557600 1000 examples/solar-system-inner.npy output/parallel/3rdLaw/solar-system-innerActual.npy 8" "examples/solar-system-inner-expected.npy" "output/parallel/3rdLaw/solar-system-innerActual.npy"
run_simulation_and_compare "./nbody-p3" $LOG_FILE_PARALLEL_3RD "60 31557600 1000 examples/solar-system.npy output/parallel/3rdLaw/solar-systemActual.npy 8" "examples/solar-system-expected.npy" "output/parallel/3rdLaw/solar-systemActual.npy"
run_simulation_and_compare "./nbody-p3" $LOG_FILE_PARALLEL_3RD "1 302400 1000 examples/pluto-charon.npy output/parallel/3rdLaw/pluto-charonActual.npy 8" "examples/pluto-charon-expected.npy" "output/parallel/3rdLaw/pluto-charonActual.npy"
run_simulation_and_compare "./nbody-p3" $LOG_FILE_PARALLEL_3RD "1e-4 2.1 1000 examples/figure8.npy output/parallel/3rdLaw/figure8Actual.npy 8" "examples/figure8-expected.npy" "output/parallel/3rdLaw/figure8Actual.npy"
run_simulation_and_compare "./nbody-p3" $LOG_FILE_PARALLEL_3RD "1e-4 75 1000 examples/figure8-rotate.npy output/parallel/3rdLaw/figure8-rotateActual.npy 8" "examples/figure8-rotate-expected.npy" "output/parallel/3rdLaw/figure8-rotateActual.npy"
run_simulation_and_compare "./nbody-p3" $LOG_FILE_PARALLEL_3RD "0.1 500000 1000 examples/random25.npy output/parallel/3rdLaw/random25Actual.npy 8" "examples/random25-expected.npy" "output/parallel/3rdLaw/random25Actual.npy"
run_simulation_and_compare "./nbody-p3" $LOG_FILE_PARALLEL_3RD "0.01 1000 1000 examples/random100.npy output/parallel/3rdLaw/random100Actual.npy 8" "examples/random100-expected.npy" "output/parallel/3rdLaw/random100Actual.npy"

echo "All simulations and comparisons completed."
