# SpMV using the CSR format

This is a sparse matrix-vector multiplication (SpMV)
program that reads a sparse matrix from a file,
and performs a Sparse Matrix-Vector Multiplication (SpMV)
operation using the CSR format for sparse matrices.

## How to Compile

Use `cmake` to build. E.g.:

```
~ $ cd spmv-CSR
~/spmv-CSR $ mkdir build
~/spmv-CSR $ cd build
~/spmv-CSR/build $ cmake -G Ninja ../src
~/spmv-CSR/build $ ninja
```

This will produce a main executable file, named `spMV`. 

## How to Run

```
./spMV <matrixName>
```

`<matrixName>` is the path to the `.mtx` file
(i.e. the matrix file as downloaded from the Matrix Market or the U. of Florida collection).
The name should be provided **without** the `.mtx` extention.
 
### Examples
Run for the [fidap005](http://math.nist.gov/MatrixMarket/data/SPARSKIT/fidap/fidap005.html)
matrix (assuming `fidap005.mtx` exists in the current directory):

```
./spMV fidap005
```
