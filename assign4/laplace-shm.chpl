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

// natural ordering to allow race conditions
proc gauss_seidel(D: domain(2), x: [D] real, epsilon: real) {
  var cnt = 0;
  const inner_domain = D.expand(-1,-1);
  var delta: real;
  var delta_temp: [n] real;

  do {
    delta_temp = 0;
    forall ij in inner_domain do {
      var temp: real = x[ij];
      x[ij] = (x[ij+(-1, 0)] + x[ij+(0, -1)] + x[ij+(1, 0)] + x[ij+(0, 1)]) / 4.0;
      delta_temp[here.id] = max(delta_temp[here.id], abs(x[ij] - temp));
    }

    // get the max value
    delta = max reduce delta_temp;

    cnt += 1;
    if (verbose) {
      writeln("Iter: ", cnt, " (delta=", delta, ")\n");
      writeln(x);
    }
  } while (delta > epsilon);
  return cnt;
}

// red-black ordering to avoid race condition
proc red_black(D: domain(2), x: [D] real, epsilon: real) {
  var cnt = 0;
  const inner_domain = D.expand(-1,-1);
  var delta: real;
  var delta_temp: [n] real;

  do {
    delta_temp = 0;
    forall ij in inner_domain do {
      if (ij[1] % 2 == ij[2] % 2) {
        var temp: real = x[ij];
        x[ij] = (x[ij+(-1, 0)] + x[ij+(0, -1)] + x[ij+(1, 0)] + x[ij+(0, 1)]) / 4.0;
        delta_temp[here.id] = max(delta_temp[here.id], abs(x[ij] - temp));
      }
    }

    forall ij in inner_domain do {
      if (ij[1] % 2 != ij[2] % 2) {
        var temp: real = x[ij];
        x[ij] = (x[ij+(-1, 0)] + x[ij+(0, -1)] + x[ij+(1, 0)] + x[ij+(0, 1)]) / 4.0;
        delta_temp[here.id] = max(delta_temp[here.id], abs(x[ij] - temp));
      }
    }

    delta = max reduce delta_temp;

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

  var b: [D] real = 0.0;
  b[n-1, 0..n-1] = 1.0;
  b[0..n-1, n-1] = 1.0;
  cnt = gauss_seidel(D, b, epsilon);
  writeln("Mesh size: ", n, " x ", n, ", epsilon=", epsilon,
            ", total gauss_seidel iterations: ", cnt);

  var c: [D] real = 0.0;
  c[n-1, 0..n-1] = 1.0;
  c[0..n-1, n-1] = 1.0;
  cnt = red_black(D, c, epsilon);
  writeln("Mesh size: ", n, " x ", n, ", epsilon=", epsilon,
            ", total red_black iterations: ", cnt);
}
