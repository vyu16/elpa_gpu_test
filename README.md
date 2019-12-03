## Background

This is a mini-app to test the performance of the ELPA eigensolver library. It
is adapted from the test program of the 2013 release of ELPA. The ELPA2 solver
features a unique two-stage diagonalization that solves a standard symmetric
(Hermitian) eigenproblem (A C = C lambda) in 5 computational steps.

1. Transform A to banded matrix B
2. Transform B to tridiagonal matrix T
3. Solve tridiagonal matrix T
4. Back-transform eigenvectors to B
5. Back-transform eigenvectors to A

## How to use

1. [Download](https://gitlab.mpcdf.mpg.de/elpa/elpa) and install ELPA. The
   following versions of ELPA are compatible with this app: 2018.05, 2018.11,
   and 2019.05. GPU (CUDA) must be enabled for ELPA. MPI may be enabled
   optionally.

2. Four test cases may be found in the subfolders, with simple makefiles
   provided. The test cases are:

  * serial real
  * serial complex
  * MPI real
  * MPI complex

3. Each test expects 5 command line arguments:

  * arg 1: Size of test matrix
  * arg 2: Number of eigenvectors to compute
  * arg 3: Block size of BLACS style 2D block-cyclic distribution
  * arg 4: 1 for ELPA 1-stage solver, 2 for ELPA 2-stage solver
  * arg 5: 1 for CPU-only, 2 for GPU acceleration

4. The mini-app reports total time to solution as well as two measures of error,
   namely (AC - C lambda) and (C^T C - I). Both should be small.

## Examples

1. Serial:
```
./test_serial_real.x 200 200 128 2 2
```

2. MPI:
```
mpirun -n 16 ./test_mpi_real.x 2000 1000 128 2 2
```
