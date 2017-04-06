//-------------------------------------------------------------------------
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University.
//-------------------------------------------------------------------------

// A sequential prime-finding algorithm.
//
// Usage: ./prime <N>
//
//
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

int main(int argc, char **argv) {
  int N, numOfThreads;

  /* check command line first */
  if (argc < 3) {
    printf ("Usage: ./prime <N> <Number of threads>\n");
    exit(0);
  }
  if ((N=atoi(argv[1])) < 2) {
    printf ("N must be greater than 1\n");
    exit(0);
  }
  if ((numOfThreads = atoi(argv[2])) < 1) {
    printf("Number of threads must be greater than 0\n");
    exit(0);
  }

  printf("Finding primes in range 1..%d\n", N);

  int *array = (int *) malloc(sizeof(int) * (N+1));

  omp_set_num_threads(numOfThreads);

  int i, j, tid;

  #pragma omp parallel for private(tid, i)
  for (i = 2; i <= N; i++) {
    //printf("Loop 1 - Thread ID: %d\n", tid = omp_get_thread_num());
    array[i] = 1;
  }

  int limit = (int) sqrt((double) N);

  #pragma omp parallel for private(tid, i, j)
  for (i = 2; i <= limit; i++) {
    //printf("Loop 2 - Thread ID: %d\n", tid = omp_get_thread_num());
    if (array[i] == 1) {
      for (j = i+i; j <= N; j += i)
        array[j] = 0;
    }
  }

  int cnt = 0;

  #pragma omp parallel for reduction(+:cnt) private(tid, i)
  for (i = 2; i <= N; i++) {
    //printf("Loop 3 - Thread ID: %d\n", tid = omp_get_thread_num());
    if (array[i] == 1)
      cnt++;
  }

  printf("Total %d primes found\n", cnt);
}


