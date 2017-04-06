//-------------------------------------------------------------------------
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University.
//-------------------------------------------------------------------------

// Jacobi method for solving a Laplace equation.
//
// Usage: ./laplace [N]
//
//
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define EPSILON 0.001   // convergence tolerance
#define VERBOSE 0       // printing control

// Initialize the mesh with a fixed set of boundary conditions.
//
void init_array(int n, double a[n][n])  {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      a[i][j] = 0;
  }
  for (i = 1; i < n-1; i++) {
    a[n-1][i] = 1.0;
    a[i][n-1] = 1.0;
  }
}

// Display the whole mesh.
//
void print_array(int n, double a[n][n])  {
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++)
      printf("%8.4f ", a[i][j]);
    printf("\n");
  }
}

// Jacobi iteration -- return the iteration count.
//
int jacobi(int n, double x[n][n], double epsilon) {
  double xnew[n][n];    // buffer for new values
  double delta;         // measure of convergence
  int cnt = 0;          // iteration counter
  int i, j;

  do {
    delta = 0.0;
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
        xnew[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
        delta = fmax(delta, fabs(xnew[i][j] - x[i][j]));
      }
    }
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
        x[i][j] = xnew[i][j];
      }
    }
    cnt++;
    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", cnt, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return cnt;
}

// natural ordering to allow race conditions
int gauss_seidel(int n, double x[n][n], double epsilon) {
  double delta;
  // temp variable to get rid of buffers
  double temp;
  int count = 0;
  int i, j;

  do {
    delta = 0.0;
    for (i = 1; i < n-1; i++) {
      for (j = 1; j < n-1; j++) {
        // temporarliy save the old value
        temp = x[i][j];
        x[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
        // apply the temp value to get
        delta = fmax(delta, fabs(x[i][j] - temp));
      }
    }
    count++;
    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", count, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return count;
}

// red-black ordering to avoid race condition
int red_black(int n, double x[n][n], double epsilon) {
  double delta;
  double temp;
  int count = 0;
  int i, j;

  do {
    delta = 0.0;
    for (i = 1; i < n-1; i++) {
      // iteration is modified to avoid race condition
      for (j = (((i + count) % 2) + 1); j < n-1; j += 2) {
        temp = x[i][j];
        x[i][j] = (x[i-1][j] + x[i][j-1] + x[i+1][j] + x[i][j+1]) / 4.0;
        delta = fmax(delta, fabs(x[i][j] - temp));
      }
    }
    count++;
    if (VERBOSE) {
      printf("Iter %d: (delta=%6.4f)\n", count, delta);
      print_array(n, x);
    }
  } while (delta > epsilon);
  return count / 2 + count % 2;
}

// Main routine.
//
int main(int argc, char **argv) {

  int n = 8;            // mesh size, default 8 x 8
  if (argc > 1) {       // check command line for overwrite
    if ((n = atoi(argv[1])) < 2) {
      printf("Mesh size must must be greater than 2, use default\n");
      n = 8;
    }
  }

  double a[n][n];       // mesh array
  init_array(n, a);

  // Jacobi iteration, return value is the total iteration number
  int cnt = jacobi(n, a, EPSILON);
  printf("Mesh size: %d x %d, epsilon=%6.4f, total Jacobi iterations: %d\n",
         n, n, EPSILON, cnt);
  if (VERBOSE)
    print_array(n, a);

  // gauss_siedel result
  double b[n][n];
  init_array(n, b);
  cnt = gauss_seidel(n, b, EPSILON);

  printf("Mesh size: %d x %d, epsilon=%6.4f, total gauss_seidel iterations: %d\n",
         n, n, EPSILON, cnt);
  if (VERBOSE)
    print_array(n, b);

  // red_black result
  double c[n][n];
  init_array(n, c);
  cnt = red_black(n, c, EPSILON);

  printf("Mesh size: %d x %d, epsilon=%6.4f, total red_black iterations: %d\n",
         n, n, EPSILON, cnt);
  if (VERBOSE)
    print_array(n, c);
}
