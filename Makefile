assn3: nBody.c submit.c 
	mpicc -std=gnu99 -o assn3 nBody.c submit.c -lm

clean: 
	rm *.o
	rm assn3

