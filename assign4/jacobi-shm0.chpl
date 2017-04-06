//------------------------------------------------------------------------- 
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University.
//------------------------------------------------------------------------- 

// Jacobi method for solving a Laplace equation.  
// (shared-memory version)
//

config const epsilon = 0.001;	// convergence tolerance
config const n = 8; 	        // mesh size (including boundary)

// Jacobi iteration -- return the iteration count.
// 
proc jacobi(D: domain(2), x: [D] real, epsilon: real) { 
  var cnt = 0;			// iteration counter

  // ... need code here ...

  return cnt;
}

// Main routine.
//
proc main() {
  const D = {0..n-1, 0..n-1};   // domain including boundary points
  var a: [D] real = 0.0;	// mesh array
  a[n-1, 0..n-1] = 1.0;         // - setting boundary values
  a[0..n-1, n-1] = 1.0;
  var cnt = jacobi(D, a, epsilon);
  writeln("Mesh size: ", n, " x ", n, ", epsilon=", epsilon, 
            ", total Jacobi iterations: ", cnt);
}
