#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

const double G  = 6.67259e-7;  /* Gravitational constant (should be e-10 but modified to get more action */
const double dt = 1.0;         /* Length of timestep */

struct Body {
  double mass;
  double X;
  double Y;
  double X_old;
  double Y_old;
  double Vx;
  double Vy;
  double Fx;
  double Fy;
}

// Distance between points with coordinates (px,py) and (qx,qy)
double dist(double px, double py, double qx, double qy) {
  return sqrt (pow(px-qx,2)+pow(py-qy,2) );
}

/* Computes forces between bodies */
void ComputeForce(int M, Body *B, Body *inbuf) {

}

int main(int argc, char **argv) {
  const int N=100;                   // Number of bodies 
  const int timesteps = 1000;        // Number of timesteps
  const double size = 100.0;           // Initial positions are in the range [0,100]

  // While

    // t = 0
    // Initialize local mass, X_old and Y_old

    // Compute forces B x B

    // Send B to next

    // Receive inbuf from previous

    // Compute B x inbuf

    // Initialize local velocity

    // t > 0
    // Update local position

    // Compute forces B x B

    // Send B to next

    // Receive inbuf from previous

    // Compute B x inbuf

    // Update local velocity

    // Update local old position

  // All process send local body to 0
}

