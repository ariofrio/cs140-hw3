/*
Assignment 3 
Team Member 1 :
Team Member 2 :
*/

#include "nBody.h"

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
}

void gennbody(double** s, double** v, double* m, int n) {
  
	printf("Generate nBody initial condition here.\n");
}

void nbody(double** s, double** v, double* m, int n, int iter, int timestep) {
	int myrank;
	int nprocs;
	int i;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int size = n/nprocs;
	
	// This is an example of printing the body parameters to the stderr. Your code should print out the final body parameters
	// in the exact order as the input file. Since we are writing to the stderr in this case, rather than the stdout, make
	// sure you dont add extra debugging statements in stderr.

  for (i = 0; i < size; i++) {
    fprintf(stderr, OUTPUT_BODY,
        s[i][0], s[i][1], s[i][2], v[i][0], v[i][1], v[i][2], m[i]);
  }
}

