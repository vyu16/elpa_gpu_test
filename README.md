## Background

It appears that there is a bug in the GPU-accelerated version of the ELPA
2-stage eigensolver (aka ELPA2). This is a mini reproducer adapted from the test
program of the 2013 version of ELPA.

ELPA2 features a unique 2-stage diagonalization that solves a standard symmetric
(Hermitian) eigenproblem (A C = C lambda) in 5 computational steps.

1. Transform A to banded matrix B
2. Transform B to tridiagonal matrix T
3. Solve tridiagonal matrix T
4. Back-transform eigenvectors to B
5. Back-transform eigenvectors to A

All 5 steps can be GPU-accelerated. To use the GPU kernel of step 4, the block
size of the BLACS style 2D block-cyclic distribution must be 128. When this is
not the case, step 4 will be executed on CPU.

We believe there is a bug in the GPU code of step 4, because

1. The eigenvalues are always correct regardless of the usage of GPU in each
individual step.
2. The eigenvalues and eigenvectors are always correct regardless of the usage
of GPU acceleration in steps 1-3 and 5, as long as GPU acceleration for step 4
is not enabled.
3. When GPU acceleration for step 4 is enabled, the eigenvectors are incorrect
for some specific combination of matrix size and process count, regardless of
the usage of GPU acceleration in steps 1-3 and 5.

## How to reproduce

1. [Download](https://gitlab.mpcdf.mpg.de/elpa/elpa) and install ELPA. This
reproducer should work with the following versions of ELPA: 2018.05.001,
2018.11.001, 2019.05.001.rc1. MPI and GPU (CUDA) must be enabled for ELPA.

2. Compile this test program. A simple Makefile is provided.

3. Run the `test_elpa.x` executable. It expects 5 command line arguments:
  * arg 1: Size of test matrix
  * arg 2: Number of eigenvectors to compute
  * arg 3: Block size of BLACS style 2D block-cyclic distribution
  * arg 4: 1 for ELPA 1-stage solver, 2 for ELPA 2-stage solver
  * arg 5: 1 for CPU-only, 2 for GPU acceleration

4. The test program reports total time to solution as well as two measures of
error. The first one corresponds to (AC - C lambda), the second one corresponds
to (C^T C - I).

5. The bug may be reproduced by, e.g.,
```
mpirun -n 16 ./test_elpa.x 2000 1000 128 2 2
```
