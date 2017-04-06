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
  int N;

  /* check command line first */
  if (argc < 2) {
    printf ("Usage: ./prime <N>\n");
    exit(0);
  }
  if ((N=atoi(argv[1])) < 2) {
    printf ("N must be greater than 1\n");
    exit(0);
  }

  printf("Finding primes in range 1..%d\n", N);

  int *array = (int *) malloc(sizeof(int) * (N+1));
  for (int i = 2; i <= N; i++)
    array[i] = 1;

  int limit = (int) sqrt((double) N);
  for (int i = 2; i <= limit; i++) {
    if (array[i] == 1) {
      for (int j = i+i; j <= N; j += i) 
	array[j] = 0;
    }
  }

  int cnt = 0;
  for (int i = 2; i <= N; i++) {
    if (array[i] == 1)
      cnt++;
  }

  printf("Total %d primes found\n", cnt);
}

