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
 * AUTHORS: Yousuf Kanan and Derek ALmon
 */

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "util.h"

typedef struct {
     double position[3];
     double velocity[3];
     double force[3];
     double mass;
} body;

// Gravitational Constant in N m^2 / kg^2 or m^3 / kg / s^2
#define G 6.6743015e-11

// Softening factor to reduce divide-by-near-zero effects
#define SOFTENING 1e-9

void calculate_force(body* body1, body* body2) {
    double dx = body2->position[0] - body1->position[0];
    double dy = body2->position[1] - body1->position[1];
    double dz = body2->position[2] - body1->position[2];
    double dist_sq = dx*dx + dy*dy + dz*dz + SOFTENING*SOFTENING;
    double dist = sqrt(dist_sq);
    double force_mag = G * body1->mass * body2->mass / dist_sq;

    double fx = force_mag * dx / dist;
    double fy = force_mag * dy / dist;
    double fz = force_mag * dz / dist;

    body1->force[0] += fx;
    body1->force[1] += fy;
    body1->force[2] += fz;
    body2->force[0] -= fx;
    body2->force[1] -= fy;
    body2->force[2] -= fz;
}

void updateForces(body* bodies, size_t n) {
    for (size_t i = 0; i < n; i++) {
        bodies[i].force[0] = bodies[i].force[1] = bodies[i].force[2] = 0.0;
    }

    for (size_t i = 0; i < n - 1; i++) {
        for (size_t j = i + 1; j < n; j++) {
            calculate_force(&bodies[i], &bodies[j]);
        }
    }
}

int main(int argc, const char *argv[])
{
    // parse arguments
    if (argc != 6 && argc != 7)
    {
        fprintf(stderr, "usage: %s time-step total-time outputs-per-body input.npy output.npy [num-threads]\n", argv[0]);
        return 1;
    }
    double time_step = atof(argv[1]), total_time = atof(argv[2]);
    if (time_step <= 0 || total_time <= 0 || time_step > total_time)
    {
        fprintf(stderr, "time-step and total-time must be positive with total-time > time-step\n");
        return 1;
    }
    size_t num_outputs = atoi(argv[3]);
    if (num_outputs <= 0)
    {
        fprintf(stderr, "outputs-per-body must be positive\n");
        return 1;
    }
    Matrix *input = matrix_from_npy_path(argv[4]);
    if (input == NULL)
    {
        perror("error reading input");
        return 1;
    }
    if (input->cols != 7)
    {
        fprintf(stderr, "input.npy must have 7 columns\n");
        return 1;
    }
    size_t n = input->rows;
    if (n == 0)
    {
        fprintf(stderr, "input.npy must have at least 1 row\n");
        return 1;
    }
    size_t num_steps = (size_t)(total_time / time_step + 0.5);
    if (num_steps < num_outputs)
    {
        num_outputs = 1;
    }
    size_t output_steps = num_steps / num_outputs;
    num_outputs = (num_steps + output_steps - 1) / output_steps;

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

    // allocate output matrix as a 3n-by-num_outputs matrix
    Matrix *output = matrix_create_raw(num_outputs, 3 * n);

    body* bodies = (body*)malloc(sizeof(body) * n);
    for (size_t i = 0; i < n; i++) {
        bodies[i].position[0] = input->data[i * 7 + 1];
        bodies[i].position[1] = input->data[i * 7 + 2];
        bodies[i].position[2] = input->data[i * 7 + 3];
        bodies[i].velocity[0] = input->data[i * 7 + 4];
        bodies[i].velocity[1] = input->data[i * 7 + 5];
        bodies[i].velocity[2] = input->data[i * 7 + 6];
        bodies[i].mass = input->data[i * 7 + 0];
    }

    

    // make row 0 the iniitial position
    for (size_t i = 0; i < n; i++)
    {
        output->data[i * 3] = input->data[i * 7 + 1];
        output->data[i * 3 + 1] = input->data[i * 7 + 2];
        output->data[i * 3 + 2] = input->data[i * 7 + 3];
        
    }
    
for (size_t i = 0; i < n; i++) {
    for (size_t d = 0; d < 3; d++) {
        // Convert initial forces to accelerations and apply half the timestep's acceleration to the velocity
        double acceleration = bodies[i].force[d] / bodies[i].mass;
        bodies[i].velocity[d] += acceleration * (time_step / 2.0);
    }
}

size_t outputIndex = 0;
for (size_t step = 0; step < num_steps; step++) {
    for (size_t i = 0; i < n; i++) {
        for (size_t d = 0; d < 3; d++) {
            // Update position based on velocity
            bodies[i].position[d] += bodies[i].velocity[d] * time_step;    
            }
    }

    // Force update based on new positions
    updateForces(bodies, n);

    // Final velocity update using the newly calculated forces (accelerations)
    for (size_t i = 0; i < n; i++) {
        for (size_t d = 0; d < 3; d++) {
            double newAcceleration = bodies[i].force[d] / bodies[i].mass;
            bodies[i].velocity[d] += newAcceleration * (time_step / 2.0);
        }
    }

    // Save positions to the output matrix at specified intervals
    if (step % output_steps == 0 || step == num_steps - 1) {
        for (size_t i = 0; i < n; i++) {
            output->data[outputIndex + i * 3 + 0] = bodies[i].position[0];
            output->data[outputIndex + i * 3 + 1] = bodies[i].position[1];
            output->data[outputIndex + i * 3 + 2] = bodies[i].position[2];
        }
        outputIndex += 3 * n;
    }
}

        

    // get the end and computation time
    clock_gettime(CLOCK_MONOTONIC, &end);
    double time = get_time_diff(&start, &end);
    printf("%f secs\n", time);

    // save results
    matrix_to_npy_path(argv[5], output);

    // cleanup
    matrix_free(input);
    matrix_free(output);
    free(bodies);

    return 0;
}
