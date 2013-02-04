#include "nBody.h"

int main(int argc, char *argv[]) {

	MPI_Init(&argc,&argv);
	int myrank;
	int nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	int n, iters, timestep;

	double** s; // position in 3D space for each body
	double** v; // velocity in 3D space for each body
	double* m;  // mass of each body
	int size; // # of bodies stored on each proc.

	// arguments: ./nBody r #n #iter #timestep
	// or         ./nBody g #n #iter #timestep
	if (argc != 5) {
		fprintf(stderr, "Usage: %s MODE BODIES ITERATIONS TIMESTEP\n", argv[0]);
    fprintf(stderr, "\nMODE must be either r[ead] or g[enerate]\n");
		exit(0);
	}
	n = atoi(argv[2]); // total number of bodies
	size = n / nprocs; // number of bodies in this processor
	iters = atoi(argv[3]);
	timestep = atoi(argv[4]); 
	size = n / nprocs;
	
	int i; 
	int j;

  // position vectors for bodies in this processor
	s = (double **)malloc(sizeof(double *) * size);
	for (i = 0; i < size; i++) {
		s[i] = (double*)malloc(sizeof(double) * 3);

		for(j = 0; j < 3; j++) {
			s[i][j] = 0;
		}
	}

  // velocity vectors for bodies in this processor
	v = (double **)malloc(sizeof(double *) * size);
	for (i = 0; i < size; i++) {
		v[i] = (double*)malloc(sizeof(double) * 3);
		
		for(j = 0; j < 3; j++) {
			v[i][j] = 0;
		}
	}

  // mass values for bodies in this processor
	m = (double *)malloc(sizeof(double) * size);

	for(i = 0; i < size; i++) {
		m[i] = 0;
	}

	if (strcmp(argv[1], "r") == 0) {
		readnbody(s, v, m, n); 
	} else {
		gennbody(s, v, m, n);
	}

	nbody(s, v, m, n, iters, timestep);

	for (i = 0; i < size; i++) {
		free(s[i]);
		free(v[i]);
	}

	free(s);
	free(v);
	free(m);

	MPI_Finalize();
	return 0;
}
