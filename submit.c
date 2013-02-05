/*
Assignment 3 
Team Member 1 :
Team Member 2 :
*/

#include "nBody.h"
#include <math.h>
#include <string.h>

#define G 6.674e-11 // gravitational constant, m^3 kg^-1 s^-2
typedef struct {
  double s[3];
  double v[3];
  double m;
} body;

void readnbody(double** ss, double** vs, double* ms, int n) {
	int myrank;
	int nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int size = n / nprocs;
	
	// This is an example of reading the body parameters from the input file. 
	if (myrank == 0) {
    int i;
		for (i = 0; i < n; i++) {
      double s[3];
      double v[3];
      double m;

			int result = scanf(INPUT_BODY, 
          &s[0], &s[1], &s[2], &v[0], &v[1], &v[2], &m);
			if (result != 7) {
				fprintf(stderr, "error reading body %d. Check if the number of bodies is correct.\n", i);
				exit(0);
			}

      if(i/size == 0) {
        ss[i][0] = s[0];
        ss[i][1] = s[1];
        ss[i][2] = s[2];
        vs[i][0] = v[0];
        vs[i][1] = v[1];
        vs[i][2] = v[2];
        ms[i] = m;
      } else {
        // This limits size to approximately 32767/6 = 5461
        MPI_Send(s, 3, MPI_DOUBLE, i/size, i%size, MPI_COMM_WORLD);
        MPI_Send(v, 3, MPI_DOUBLE, i/size, i%size+size*2, MPI_COMM_WORLD);
        MPI_Send(&m, 1, MPI_DOUBLE, i/size, i%size+size*4, MPI_COMM_WORLD);
      }
		}
	} else {
    MPI_Status status;
    int i;
    for(i = 0; i < size; i++) {
      MPI_Recv(ss[i], 3, MPI_DOUBLE, 0, i, MPI_COMM_WORLD, &status);
      MPI_Recv(vs[i], 3, MPI_DOUBLE, 0, i + size*2, MPI_COMM_WORLD, &status);
      MPI_Recv(&ms[i], 1, MPI_DOUBLE, 0, i + size*4, MPI_COMM_WORLD, &status);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
}

void gennbody(double** s, double** v, double* m, int n) {
  
	printf("Generate nBody initial condition here.\n");
}

void nbody(double** ss, double** vs, double* ms, 
    int n, int iter, int timestep) {
	int myrank;
	int nprocs;
  MPI_Status status;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int size = n/nprocs;

  // Copy data to our MPI-friendly data structure.
  body* mine = (body *)malloc(sizeof(body) * size);
  body* theirs = (body *)malloc(sizeof(body) * size);
  body* buffer = (body *)malloc(sizeof(body) * size);
  for(int i = 0; i < size; i++) {
    mine[i].s[0] = ss[i][0];
    mine[i].s[1] = ss[i][1];
    mine[i].s[2] = ss[i][2];
    mine[i].v[0] = vs[i][0];
    mine[i].v[1] = vs[i][1];
    mine[i].v[2] = vs[i][2];
    mine[i].m = ms[i];
  }

  double* accel = (double *)malloc(sizeof(double) * size * 3);
  for(int i = 0; i < size; i++) {
    accel[i*3 + 0] = 0;
    accel[i*3 + 1] = 0;
    accel[i*3 + 2] = 0;
  }

  for(int i = 0; i < iter; i++) {
    // Calculate my own bodies' effects on each other
    for(int j = 0; j < size; j++) {
      for(int k = 0; k < size; k++) {
        double r = sqrt(
            (mine[j].s[0] - mine[k].s[0]) * (mine[j].s[0] - mine[k].s[0]) +
            (mine[j].s[1] - mine[k].s[1]) * (mine[j].s[1] - mine[k].s[1]) +
            (mine[j].s[2] - mine[k].s[2]) * (mine[j].s[2] - mine[k].s[2])
            );
        if(r != 0 && j != k) {
          accel[j*3 + 0] -= G * mine[k].m *
                            (mine[j].s[0] - mine[k].s[0]) / r / r / r;
          accel[j*3 + 1] -= G * mine[k].m *
                            (mine[j].s[1] - mine[k].s[1]) / r / r / r;
          accel[j*3 + 2] -= G * mine[k].m *
                            (mine[j].s[2] - mine[k].s[2]) / r / r / r;
        }
      }
    }

    if(nprocs > 1) {
      if(myrank % 2) { // avoid deadlock
        MPI_Send(mine, sizeof(body) * size, MPI_BYTE,
            (myrank + 1) % nprocs, 0, MPI_COMM_WORLD);
        MPI_Recv(theirs, sizeof(body) * size, MPI_BYTE,
            (myrank - 1) % nprocs, 0, MPI_COMM_WORLD, &status);
      } else {
        MPI_Recv(theirs, sizeof(body) * size, MPI_BYTE,
            (myrank - 1) % nprocs, 0, MPI_COMM_WORLD, &status);
        MPI_Send(mine, sizeof(body) * size, MPI_BYTE,
            (myrank + 1) % nprocs, 0, MPI_COMM_WORLD);
      }
    }

    for(int round = 1; round < nprocs; round++) {
      for(int j = 0; j < size; j++) {
        for(int k = 0; k < size; k++) {
          double r = sqrt(
              (mine[j].s[0] - theirs[k].s[0]) * (mine[j].s[0] - theirs[k].s[0]) +
              (mine[j].s[1] - theirs[k].s[1]) * (mine[j].s[1] - theirs[k].s[1]) +
              (mine[j].s[2] - theirs[k].s[2]) * (mine[j].s[2] - theirs[k].s[2])
              );
          if(r != 0) {
            accel[j*3 + 0] -= G * theirs[k].m *
                              (mine[j].s[0] - theirs[k].s[0]) / r / r / r;
            accel[j*3 + 1] -= G * theirs[k].m *
                              (mine[j].s[1] - theirs[k].s[1]) / r / r / r;
            accel[j*3 + 2] -= G * theirs[k].m *
                              (mine[j].s[2] - theirs[k].s[2]) / r / r / r;
          }
        }
      }

      if(round < nprocs - 1) {
        if(myrank % 2) { // avoid deadlock
          MPI_Send(theirs, sizeof(body) * size, MPI_BYTE,
              (myrank + 1) % nprocs, round, MPI_COMM_WORLD);
          MPI_Recv(theirs, sizeof(body) * size, MPI_BYTE,
              (myrank - 1) % nprocs, round, MPI_COMM_WORLD, &status);
        } else {
          MPI_Recv(buffer, sizeof(body) * size, MPI_BYTE,
              (myrank - 1) % nprocs, round, MPI_COMM_WORLD, &status);
          MPI_Send(theirs, sizeof(body) * size, MPI_BYTE,
              (myrank + 1) % nprocs, round, MPI_COMM_WORLD);
          memcpy(theirs, buffer, sizeof(body) * size);
        }
      }
    }

    // Let acceleration affect velocity and position
    for(int j = 0; j < size; j++) {
      mine[j].v[0] += timestep * accel[j*3 + 0];
      mine[j].v[1] += timestep * accel[j*3 + 1];
      mine[j].v[2] += timestep * accel[j*3 + 2];
      mine[j].s[0] += timestep * mine[j].v[0];
      mine[j].s[1] += timestep * mine[j].v[1];
      mine[j].s[2] += timestep * mine[j].v[2];
    }
  }
  // This is an example of printing the body parameters to the stderr. Your
  // code should print out the final body parameters in the exact order as the
  // input file. Since we are writing to the stderr in this case, rather than
  // the stdout, make sure you dont add extra debugging statements in stderr.

  MPI_Barrier(MPI_COMM_WORLD);
  for(int rank = 0; rank < nprocs; rank++) {
    if(myrank == 0) {
      if(rank != 0) {
        MPI_Recv(
            mine, sizeof(body) * size, MPI_BYTE, rank, 0,
            MPI_COMM_WORLD, &status);
      }
      for(int i = 0; i < size; i++) {
        fprintf(stderr, OUTPUT_BODY,
            mine[i].s[0], mine[i].s[1], mine[i].s[2],
            mine[i].v[0], mine[i].v[1], mine[i].v[2],
            mine[i].m);
      }
    } else {
      MPI_Send(
          mine, sizeof(body) * size, MPI_BYTE, 0, 0,
          MPI_COMM_WORLD);
    }
  }
}

