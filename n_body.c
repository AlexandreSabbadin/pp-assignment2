#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <mpi.h>

const double G  = 6.67259e-7;  /* Gravitational constant (should be e-10 but modified to get more action */
const double dt = 1.0;         /* Length of timestep */

double pos0x;
double pos0y;

typedef struct Body {
  double mass;
  double X;
  double Y;
  double X_old;
  double Y_old;
  double Vx;
  double Vy;
  double Fx;
  double Fy;
} Body;

/* Writes out positions (x,y) of N particles to the file fn 
   Returns zero if the file couldn't be opened, otherwise 1 */
int write_particles(int N, Body *buffer, char *fn) {
  FILE *fp;
  /* Open the file */
  if ((fp=fopen(fn, "w")) == NULL) {
    printf("Couldn't open file %s\n", fn);
    return 0;
  }
  /* Write the positions to the file fn */
  for (int i = 0; i < N; i++) {
    fprintf(fp, "%3.2f %3.2f \n", buffer[i].X, buffer[i].Y);
  }
  fprintf(fp, "\n" );
  fclose(fp);  /* Close the file */
  return(1);
}

// Distance between points with coordinates (px,py) and (qx,qy)
double dist(double px, double py, double qx, double qy) {
  return sqrt (pow(px-qx,2)+pow(py-qy,2) );
}

// Computes forces between bodies
void ComputeForce(int M, Body *B, Body *inbuf) {
  const double mindist  = 0.0001;  // Minimal distance of two bodies of being in interaction
  
  for (int i = 0; i < M; i++) {      // Compute the force for all bodies
    for (int j = 0; j < M; j++) {   // The force on a body i is the sum of forces from all other bodies j
      // Distance between points i and j
      double r = dist(B[i].X_old, B[i].Y_old, inbuf[j].X_old, inbuf[j].Y_old); 
      
      if (r > mindist) {        // Very near-distance forces are ignored
        double r3 = pow(r, 3);     // Could also be written as r3=r*r*r;
        B[i].Fx += G * B[i].mass * inbuf[j].mass * (inbuf[j].X_old - B[i].X_old) / r3;
        B[i].Fy += G * B[i].mass * inbuf[j].mass * (inbuf[j].Y_old - B[i].Y_old) / r3;
      }
    }
  }
}

int main(int argc, char **argv) {
  const int N=400;                   // Number of bodies 
  const int timesteps = 1000;        // Number of timesteps
  const double size = 100.0;           // Initial positions are in the range [0,100]

  double start;

  // Initialize MPI
  int P, k;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &P);
  MPI_Comm_rank(MPI_COMM_WORLD, &k);

  // Body MPI_Datatype
  MPI_Datatype MPI_BODY;
  MPI_Type_contiguous(9, MPI_DOUBLE, &MPI_BODY);
  MPI_Type_commit(&MPI_BODY);

  /* Check that we run on more than one process */
  if (P < 2) {
    printf("You have to use at least 2 processes to run this program\n");
    MPI_Finalize();            // Quit if there is only one process
    exit(0);
  }

  // Local body array
  int M = N / P;
  Body *B = calloc(M, sizeof(Body));

  // Buffer body array
  Body *inbuf = calloc(M, sizeof(Body));

  // Process 0
  if (k == 0) {
    printf ("N-body simulation, N = %d, P = %d, M = %d, timesteps = %d\n\n", N, P, M, timesteps);

    // Start timer
    start = MPI_Wtime();

    Body *buffer = calloc(N, sizeof(Body));

    // Seed the random number generator so that it generates a fixed sequence
    short int seedval[3] = {7, 7, 7};
    seed48(seedval);

    for (int i = 0; i < N; i++) {
      buffer[i].mass = 1000.0 * drand48();
      buffer[i].X_old = buffer[i].X = size * drand48();
      buffer[i].Y_old = buffer[i].Y = size * drand48();
    }

    // Send initialized body arrays
    Body *inbuf = calloc(M, sizeof(Body));
    for (int i = 1; i < P; i++) {
      for (int j = 0; j < M; j++) {
        inbuf[j] = buffer[i * M + j];
      }
      MPI_Send(inbuf, M, MPI_BODY, i, 0, MPI_COMM_WORLD);
      printf("[0] Initialized B sent to process [%d]\n", i); fflush(stdout);
    }
    free(inbuf);

    // Initialize B
    for (int i = 0; i < M; i++) {
      B[i] = buffer[i];
    }

    // Write buffer to file
    write_particles(N, buffer, "initial_pos.txt");

    // Save position of one body so we can see where it has moved
    pos0x = buffer[0].X_old;
    pos0y = buffer[0].Y_old;

    free(buffer);
  } else {
    // Process > 0

    MPI_Recv(inbuf, M, MPI_BODY, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("[%d] B received from process [0]\n", k); fflush(stdout);

    for (int i = 0; i < M; i++) {
      B[i] = inbuf[i];
    }
  }

  // Compute forces B x B
  ComputeForce(M, B, B);

  for (int step = 1; step < P; step++) {
    int next = (k + step) % P;
    int previous = (k - step + P) % P;

    MPI_Request r;

    // Send B to next
    MPI_Isend(B, M, MPI_BODY, next, 0, MPI_COMM_WORLD, &r);

    // Receive inbuf from previous
    MPI_Recv(inbuf, M, MPI_BODY, previous, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // Compute B x inbuf
    ComputeForce(M, B, inbuf);
  }

  // Initialize local velocity
  for(int i = 0; i < M; i++) {
    B[i].Vx = 0.5 * dt * B[i].Fx / B[i].mass;
    B[i].Vy = 0.5 * dt * B[i].Fy / B[i].mass;
  }

  fflush(stdout);

  // While
  int t = 0;
  while (t < timesteps) {
    t++;
    /*if (k == 0) { //DEBUG
      printf("%d ", t); fflush(stdout);  // Print out the timestep
    }*/

    // Update local position
    for (int i = 0; i < M; i++) {
      B[i].X = B[i].X_old + B[i].Vx * dt;
      B[i].Y = B[i].Y_old + B[i].Vy * dt;
    }

    // Compute forces B x B
    for (int i = 0; i < M; i++) {
      B[i].Fx = B[i].Fy = 0.0;
    }
    ComputeForce(M, B, B);

    for (int step = 1; step < P; step++) {
      int next = (k + step) % P;
      int previous = (k - step + P) % P;

      // Send B to next
      MPI_Request r;
      MPI_Isend(B, M, MPI_BODY, next, 0, MPI_COMM_WORLD, &r);

      // Receive inbuf from previous
      MPI_Recv(inbuf, M, MPI_BODY, previous, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

      // Compute B x inbuf
      ComputeForce(M, B, inbuf);
    } 

    // Update local velocity and old position
    for (int i = 0; i < M; i++) {		
      B[i].Vx = B[i].Vx + B[i].Fx * dt / B[i].mass;
      B[i].Vy = B[i].Vy + B[i].Fy * dt / B[i].mass;
      B[i].X_old = B[i].X;
      B[i].Y_old = B[i].Y;
    }
  }

  // All process send local body to 0
  if (k == 0) {
    printf("\n"); fflush(stdout);

    Body *buffer = calloc(N, sizeof(Body));
    for (int i = 0; i < M; i++) {
      buffer[i] = B[i];
    }

    // Receive from all process
    for (int i = 1; i < P; i++) {
      MPI_Recv(inbuf, M, MPI_BODY, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("[0] Final B received from process [%d]\n", i); fflush(stdout);
      for (int j = 0; j < M; j++) {
        buffer[i*M+j] = inbuf[j];
      }
    }

    // Write final particle coordinates to a file
    write_particles(N, buffer, "final_pos.txt");

    printf("\nOriginal pos\t%3.8f\t%3.8f\n", pos0x, pos0y);
    printf("Final pos\t%3.8f\t%3.8f\n\n", buffer[0].X, buffer[0].Y);

    free(buffer);

    // End timer
    printf("Total time: %fs.\n", MPI_Wtime() - start);
  } else {
    // Send B to process 0
    MPI_Send(B, M, MPI_BODY, 0, 0, MPI_COMM_WORLD);
    printf("[%d] Final B sent to process [0]\n", k); fflush(stdout);
  }  

  free(inbuf);
  free(B);
  printf("[%d] Process terminated\n", k); fflush(stdout);
  MPI_Finalize();
  exit(0);
}

