#ifndef _SPMV_MATRIX_H
#define _SPMV_MATRIX_H

#include <cstdlib> // defines NULL
#include <vector>

class Matrix final {
 public:
  int *cols;
  int *rows;
  double *vals;
  unsigned long n;
  unsigned long nz;
  // For some representations, numRows, numCols, numVals
  // may not be the same as n, nz.
  unsigned long numRows, numCols, numVals;
  
  Matrix(int *rows, int *cols, double *vals, unsigned long n, unsigned long nz);
  
  ~Matrix();
  
  void print();
};

class MMElement {
public:
  int row; int col; double val;
  
  MMElement(int row, int col, double val):
    row(row), col(col), val(val) {}
  
  static bool compare(const MMElement &elt1, const MMElement &elt2);
  bool compareCol(int col);
};

class MMMatrix final {
 private:
  std::vector<MMElement> elts;
  unsigned long n;
  
 public:
  MMMatrix(unsigned long n);
  
  ~MMMatrix();
  
  void add(int row, int col, double val);
  void normalize();
  void print();
  void printMTX();
  
  // Return a matrix in the CSR format
  Matrix* toCSRMatrix();
  
 private:
  void sort();
};

#endif
