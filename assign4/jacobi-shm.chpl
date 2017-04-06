//-------------------------------------------------------------------------
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University.
//-------------------------------------------------------------------------

// Jacobi method for solving a Laplace equation.
// (shared-memory version)
//

config const epsilon = 0.001;   // convergence tolerance
config const n = 8;             // mesh size (including boundary)

// printing control
config const verbose = true;

// Jacobi iteration -- return the iteration count.
//
proc jacobi(D: domain(2), x: [D] real, epsilon: real) {
  var cnt = 0;                  // iteration counter

  // ... need code here ...
  const interior_domain = D.expand(-1, -1);

  // buffer for new values
  var xnew: [D] real;
  // measure of convergence
  var delta: real;

  do {
    // ij for matrix index
    forall ij in interior_domain do
      xnew[ij] = (x[ij+(-1, 0)] + x[ij+(0, -1)] + x[ij+(1, 0)] + x[ij+(0, 1)]) / 4.0;
    delta = max reduce abs(xnew[interior_domain] - x[interior_domain]);

    // copy new values to x
    x[interior_domain] = xnew[interior_domain];

    cnt += 1;
    if (verbose) {
      writeln("Iter: ", cnt, " (delta=", delta, ")\n");
      writeln(x);
    }
  } while (delta > epsilon);
  return cnt;
}

// Main routine.
//
proc main() {
  const D = {0..n-1, 0..n-1};   // domain including boundary points
  var a: [D] real = 0.0;        // mesh array
  a[n-1, 0..n-1] = 1.0;         // - setting boundary values
  a[0..n-1, n-1] = 1.0;
  var cnt = jacobi(D, a, epsilon);
  writeln("Mesh size: ", n, " x ", n, ", epsilon=", epsilon,
            ", total Jacobi iterations: ", cnt);
}
