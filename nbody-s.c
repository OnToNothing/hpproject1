/**
 * Runs a simulation of the n-body problem in 3D.
 * 
 * To compile the program:
 *   gcc -Wall -Ofast -march=native nbody-s.c matrix.c util.c -o nbody-s -lm
 * 
 * To run the program:
 *   ./nbody-s time-step total-time outputs-per-body input.npy output.npy
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
 * AUTHORS:Derek Allmon and yousuf kanan 
 */

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include "matrix.h"
#include "util.h"

// Gravitational Constant in N m^2 / kg^2 or m^3 / kg / s^2
#define G 6.6743015e-11

// Softening factor to reduce divide-by-near-zero effects
#define SOFTENING 1e-9

// double lawOfGravity(double massi , double massj , double xi , double xj , double yi , double yj 
//                     , double zi , double zj){
//     double delta_x = xj - xi;
//     double delta_y = yj - yi;
//     double delta_z = zj - zi;
//     double r = sqrt(delta_x*delta_x + delta_y*delta_y + delta_z*delta_z + SOFTENING);
//     double rSquared = r*r; 
//     double force = G * ((massi * massj) / rSquared);
//     return force;
// }


double superpositionPrinciple(double massi, double massj,double delta ,size_t n, double euclideanDistance) {
    
    const double distanceCubed = euclideanDistance * euclideanDistance * euclideanDistance;

    // Calculate superpositions
    double superposition = 0;

    for (size_t j = 0; j <= n; ++j) {
        superposition += massj * (delta) / (distanceCubed);
    }
    superposition *= G; // Apply gravitational constant
    return superposition;
}


double computeVelocity(double time_step, double accelartion, double velocity) {
  // Half-step acceleration update
  return velocity + time_step * accelartion;
}
//run the simulation

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

    // allocate output matrix as a 3n-by-num_outputs matrix
    Matrix* output = matrix_create_raw(num_outputs, 3*n);

    //make row 0 the iniitial position
    for (size_t i = 0; i < n; i++) {
        output->data[i*3] = input->data[i*7+1];
        output->data[i*3+1] = input->data[i*7+2];
        output->data[i*3+2] = input->data[i*7+3];
    }
    Matrix* inputCopy = matrix_copy(input); 
    //run the simulation
    //MAKE INPUT DATA USABLE BY COPYING IT TO A NEW ARRAY
    
    //HAVE INDEXED I AND J LOOPS TO CALL SUPERPOSITION PRINCIPLE
    //CALL SUPERPOSITION PRINCIPLE TO GET THE FORCE
    //CALL COMPUTE VELOCITY TO GET THE VELOCITY
    //UPDATE THE POSITION OF EACH BODY
    //PERIODICALLY COPY THE POSITION TO THE OUTPUT MATRIX
    //if i == j then skip the iteration
    //FREE THE MEMORY
    double velocityx = input->data[4];
    double velocityy = input->data[5];
    double velocityz = input->data[6];
    double accelerationX;
    double accelerationY;
    double accelerationZ;
    double massi;
    double massj;
    double xi;
    double xj;
    double yi;
    double yj;
    double zi;
    double zj;
    double euclideanDistance;
    for (size_t t = 1; t < num_steps; t++) { 
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                massi = input->data[j*7];
                massj = input->data[i*7]; 
                xi = input->data[i*7+1];
                xj = input->data[j*7+1];
                yi = input->data[i*7+2];
                yj = input->data[j*7+2];
                zi = input->data[i*7+3];
                zj = input->data[j*7+3];
                euclideanDistance = sqrt((xj - xi) * (xj - xi) + (yj - yi) * (yj - yi) + (zj - zi) * (zj - zi) + SOFTENING);
                accelerationX = superpositionPrinciple(massi, massj,xj - xi, n, euclideanDistance);
                accelerationY = superpositionPrinciple(massi, massj,yj - yi, n, euclideanDistance);
                accelerationZ = superpositionPrinciple(massi, massj,zj - zi, n, euclideanDistance);
                velocityx = computeVelocity(time_step, accelerationX, velocityx);
                velocityy = computeVelocity(time_step, accelerationY, velocityy);
                velocityz = computeVelocity(time_step, accelerationZ, velocityz);
                inputCopy->data[j*7+1] += velocityx;
                inputCopy->data[j*7+2] += velocityy;
                inputCopy->data[j*7+3] += velocityz;

                if (i % output_steps == 0) {
                output->data[i*3*n+j*3] = inputCopy->data[j*7+1];
                output->data[i*3*n+j*3+1] = inputCopy->data[j*7+2];
                output->data[i*3*n+j*3+2] = inputCopy->data[j*7+3];
            }
        
        
          
        }
    }
    }
        
    
    if (num_steps % output_steps != 0) {
        //Save positions to reow "num_output-1" of output 
        for (size_t i = 0; i < n; i++) {
            output->data[(num_outputs-1)*3*n+i*3] = input->data[i*7+1];
            output->data[(num_outputs-1)*3*n+i*3+1] = input->data[i*7+2];
            output->data[(num_outputs-1)*3*n+i*3+2] = input->data[i*7+3];
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
    matrix_free(inputCopy);
    matrix_free(output);


    return 0;
}