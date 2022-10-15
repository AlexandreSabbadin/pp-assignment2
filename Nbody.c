// Compile with  gcc -O2 Nbody.c -o Nbody -lm

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

const double G  = 6.67259e-7;  /* Gravitational constant (should be e-10 but modified to get more action */
const double dt = 1.0;         /* Length of timestep */

/* Writes out positions (x,y) of N particles to the file fn 
   Returns zero if the file couldn't be opened, otherwise 1 */
int write_particles(int N, double *X, double *Y, char *fn) {
  FILE *fp;
  /* Open the file */
  if ((fp=fopen(fn, "w")) == NULL) {
    printf("Couldn't open file %s\n", fn);
    return 0;
  }
  /* Write the positions to the file fn */
  for (int i=0; i<N; i++) {
    fprintf(fp, "%3.2f %3.2f \n", X[i], Y[i] );
  }
  fprintf(fp, "\n" );
  fclose(fp);  /* Close the file */
  return(1);
}

// Distance between points with coordinates (px,py) and (qx,qy)
double dist(double px, double py, double qx, double qy) {
  return sqrt (pow(px-qx,2)+pow(py-qy,2) );
  // Could also be written as sqrt( (px-qx)*(px-qx) + (py-qy)*(py-qy) )
}


/* Computes forces between bodies */
void ComputeForce(int N, double *X, double *Y, double *mass, double *Fx, double *Fy){
  const double mindist  = 0.0001;  /* Minimal distance of two bodies of being in interaction*/
  
  for (int i = 0; i < N; i++) {      // Compute the force for all bodies
    Fx[i] = Fy[i] = 0.0;             // Initialize force vector to zero
    for ( int j = 0; j < N; j++) {   // The force on a body i is the sum of forces from all other bodies j
      if (i != j) {                  //     but not from it self
        // Distance between points i and j
        double r = dist(X[i], Y[i], X[j], Y[j]); 
        
        if (r>mindist) {        // Very near-distance forces are ignored
          double r3 = pow(r,3);     // Could also be written as r3=r*r*r;
          Fx[i] += G*mass[i]*mass[j]*(X[j]-X[i]) / r3;
          Fy[i] += G*mass[i]*mass[j]*(Y[j]-Y[i]) / r3;
        }
      }
    }
  }
}
	    

int main(int argc, char **argv) {
  
  const int N=400;                   // Number of bodies 
  const int timesteps = 1000;        // Number of timesteps
  const double size = 100.0;           // Initial positions are in the range [0,100]

  double *mass;          /* mass of bodies */
  double *X;             /* x-positions of bodies */
  double *Y;             /* y-positions of bodies */
  double *X_old;         /* old x-positions of bodies */
  double *Y_old;         /* old y-positions of bodies */
  double *Vx;            /* velocities on x-axis of bodies */
  double *Vy;            /* velocities on y-axis of bodies */
  double *Fx;            /* forces on x-axis of bodies */
  double *Fy;            /* forces on y-axis of bodies */

  printf ("N-body simulation, number of bodies = %d \n", N);
  
  /* Allocate space for variables  */ 
  mass     = (double *) calloc(N, sizeof(double));  // Mass 
  X        = (double *) calloc(N, sizeof(double));  // Position (x,y) at current time step
  Y        = (double *) calloc(N, sizeof(double));
  X_old    = (double *) calloc(N, sizeof(double));  // Position (x,y) at previous time step
  Y_old    = (double *) calloc(N, sizeof(double));
  Vx       = (double *) calloc(N, sizeof(double));  // Velocities
  Vy       = (double *) calloc(N, sizeof(double));
  Fx       = (double *) calloc(N, sizeof(double));  // Forces
  Fy       = (double *) calloc(N, sizeof(double));

  // Seed the random number generator so that it generates a fixed sequence
  short int seedval[3] = {7, 7, 7};
  seed48(seedval);

  /* Initialize mass and position of bodies */
  for (int i = 0; i<N; i++){
    mass[i]  = 1000.0*drand48();   // 0 <= mass < 1000
    X_old[i] = size*drand48();      // 0 <= x_old < 100
    Y_old[i] = size*drand48();      // 0 <= y_old < 100
  }

  /* DEBUG
  for (int i=0; i<N; i++) {
    printf("%3.1f ", Y_old[i]);
  }
  printf("\n");
  */

  // Write intial particle coordinates to a file
  write_particles(N, X_old, Y_old, "initial_pos.txt");

  // Save position of one body so we can see where it has moved
  double pos0x = X_old[0];
  double pos0y = Y_old[0];

  // Compute the initial forces that we get 
  ComputeForce(N, X_old, Y_old, mass, Fx, Fy);

  // Set up the velocity vectors caused by initial forces for Leapfrog method
  for(int i = 0; i<N; i++){
    Vx[i] = 0.5*dt*Fx[i]/mass[i];
    Vy[i] = 0.5*dt*Fy[i]/mass[i];
  }
  
  /* Main loop: 
    - Move the bodies
    - Calculate forces of the bodies with their new position
    - Calculate velocities of the bodies with the new forces
    - Copy the updated positions to the old positions (for use in next timestep)
  */
  int t=0;
  while (t<timesteps) {    // Loop for this many timesteps
    t++;
    printf("%d ", t); fflush(stdout);  // Print out the timestep
    
    // Calculate new positions 
    for (int i=0;i<N;i++){
      X[i] = X_old[i] + Vx[i]*dt;
      Y[i] = Y_old[i] + Vy[i]*dt;
    }
    
    /* Calculate forces for the new positions */
    ComputeForce(N, X_old, Y_old, mass, Fx, Fy);
    
    /* Update velocities of bodies */ 
    for (int i=0;i<N;i++){		
      Vx[i] = Vx[i] + Fx[i]*dt/mass[i];
      Vy[i] = Vy[i] + Fy[i]*dt/mass[i];
    }	
    
    /* Copy updated positions to (X_old, Y_old) before next time step begins */
    for (int i=0;i<N;i++){
      X_old[i] = X[i];
      Y_old[i] = Y[i];
    }

  }  /* end of while-loop */
  printf("\n");

  // Write final particle coordinates to a file
  write_particles(N, X, Y, "final_pos.txt");

  
    //DEBUG
  double pos0_finalx = X_old[0];
  double pos0_finaly = Y_old[0];

  printf("Original pos %3.8f %3.8f\n", pos0x, pos0y);
  printf("Final pos    %3.8f %3.8f\n", pos0_finalx, pos0_finaly);

  exit(0);
}

