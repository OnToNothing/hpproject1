/**
 * Runs a simulation of the n-body problem in 3D.
 * 
 * To compile the program:
 *   gcc -Wall -O3 -march=native nbody-s3.c matrix.c util.c -o nbody-s3 -lm
 * 
 * To run the program:
 *   ./nbody-s3 time-step total-time outputs-per-body input.npy output.npy
 * where:
 *   - time-step is the amount of time between steps (Δt, in seconds)
 *   - total-time is the total amount of time to simulate (in seconds)
 *   - outputs-per-body is the number of positions to output per body
 *   - input.npy is the file describing the initial state of the system (below)
 *   - output.npy is the output of the program (see below)
 * 
 * input.npy has a n-by-7 matrix with one row per body and the columns:
 *   - mass (in kg)
 *   - initial x, y, z position (in m)
 *   - initial x, y, z velocity (in m/s)
 * 
 * output.npy is generated and has a (outputs-per-body)-by-(3n) matrix with each
 * row containing the x, y, and z positions of each of the n bodies after a
 * given timestep.
 * 
 * See the PDF for implementation details and other requirements.
 * 
 * AUTHORS:
 */

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>

#include "matrix.h"
#include "util.h"

// Gravitational Constant in N m^2 / kg^2 or m^3 / kg / s^2
#define G 6.6743015e-11

// Softening factor to reduce divide-by-near-zero effects
#define SOFTENING 1e-9


// gravitational force between two bodies 
void forceCalulation(double* f, const double* m1, const double* m2) {
    double dx = m2[1] - m1[1]; // x position difference
    double dy = m2[2] - m1[2]; // y position difference
    double dz = m2[3] - m1[3]; // z position difference
    double r2 = dx*dx + dy*dy + dz*dz + SOFTENING;
    double r = sqrt(r2);
    double mag = G * m1[0] * m2[0] / r2 / r; // m1[0] and m2[0] are the masses
    f[0] = mag * dx / r;
    f[1] = mag * dy / r;
    f[2] = mag * dz / r;
}

double totalForce[3] = {0.0, 0.0, 0.0}; // Total force on the body

// superposition principle

void superpositionPrinciple(double* forces, const Matrix* input, size_t number_of_bodies) {
    // Assuming forces is a pre-allocated array for cumulative forces on each body
    memset(forces, 0, sizeof(double) * 3 * number_of_bodies); // Reset all forces to zero

    for (size_t i = 0; i < number_of_bodies; i++) {
        for (size_t j = 0; j < i; j++) { // Ensuring force calculation only for j < i
            double tempForce[3] = {0.0, 0.0, 0.0};
            forceCalulation(tempForce, input->data + i*7, input->data + j*7);

            // Apply the force from j to i
            forces[i*3 + 0] += tempForce[0];
            forces[i*3 + 1] += tempForce[1];
            forces[i*3 + 2] += tempForce[2];

            // Apply the opposite force from i to j, using Newton's Third Law
            forces[j*3 + 0] -= tempForce[0];
            forces[j*3 + 1] -= tempForce[1];
            forces[j*3 + 2] -= tempForce[2];
        }
    }
}


// calculate the acceleration 

void calculateAcceleration(double* acc, const double* m, const Matrix* input, size_t number_of_bodies) {
    double f[3] = {0.0, 0.0, 0.0}; // Reset net force to zero before calculation
    superpositionPrinciple(f, m, input, number_of_bodies);
    // Assuming m[3] is the mass of the body and it's not zero
    if (m[3] != 0) {
        acc[0] = f[0] / m[3]; // Acceleration in x-axis
        acc[1] = f[1] / m[3]; // Acceleration in y-axis
        acc[2] = f[2] / m[3]; // Acceleration in z-axis
    } else {
        // Handle the case where mass is zero or not provided to avoid division by zero
        acc[0] = acc[1] = acc[2] = 0.0;
    }
}

// calculate the velocity

double calculateFutureVelocityX(double vx_current, double ax, double delta_t) {
    // Calculate the future velocity along the x-axis
    double vx_future = vx_current + ax * delta_t;
    return vx_future;
}

double calculateFutureVelocityY(double vy_current, double ay, double delta_t) {
    // Calculate the future velocity along the y-axis
    double vy_future = vy_current + ay * delta_t;
    return vy_future;
}

double calculateFutureVelocityZ(double vz_current, double az, double delta_t) {
    // Calculate the future velocity along the z-axis
    double vz_future = vz_current + az * delta_t;
    return vz_future;
}




int main(int argc, const char* argv[]) {
    // parse arguments
    if (argc != 6 && argc != 7) { fprintf(stderr, "usage: %s time-step total-time outputs-per-body input.npy output.npy [num-threads]\n", argv[0]); return 1; }
    double time_step = atof(argv[1]), total_time = atof(argv[2]);
    if (time_step <= 0 || total_time <= 0 || time_step > total_time) { fprintf(stderr, "time-step and total-time must be positive with total-time > time-step\n"); return 1; }
    size_t num_outputs = atoi(argv[3]);
    if (num_outputs <= 0) { fprintf(stderr, "outputs-per-body must be positive\n"); return 1; }
    Matrix* input = matrix_from_npy_path(argv[4]);
    if (input == NULL) { perror("error reading input"); return 1; }
    if (input->cols != 7) { fprintf(stderr, "input.npy must have 7 columns\n"); return 1; }
    size_t n = input->rows;
    if (n == 0) { fprintf(stderr, "input.npy must have at least 1 row\n"); return 1; }
    size_t num_steps = (size_t)(total_time / time_step + 0.5);
    if (num_steps < num_outputs) { num_outputs = 1; }
    size_t output_steps = num_steps/num_outputs;
    num_outputs = (num_steps+output_steps-1)/output_steps;

    // variables available now:
    //   time_step    number of seconds between each time point
    //   total_time   total number of seconds in the simulation
    //   num_steps    number of time steps to simulate (more useful than total_time)
    //   num_outputs  number of times the position will be output for all bodies
    //   output_steps number of steps between each output of the position
    //   input        n-by-7 Matrix of input data
    //   n            number of bodies to simulate

    // start the clock
    struct timespec start, end;
    clock_gettime(CLOCK_MONOTONIC, &start);


    // get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f secs\n", time);

    // save results
    //matrix_to_npy_path(argv[5], output);

    // cleanup


    return 0;
}
