#include "matrix.h"
#include <iostream>
#include <algorithm>

using std::cout;

Matrix::Matrix(int *rows, int *cols, double *vals, unsigned long n, unsigned long nz):
  rows(rows), cols(cols), vals(vals), n(n), nz(nz) { 
  numRows = n;
  numCols = nz;
  numVals = nz;
}

Matrix::~Matrix() {
  delete[] rows;
  delete[] cols;
  delete[] vals;
}

void Matrix::print() {
  cout << "int numMatrixRows = " << n << ";\n";
  cout << "int numMatrixValues = " << nz << ";\n";
  if (numCols == 0) {
    cout << "int *matrixCols = 0;\n";
  } else {
    cout << "int matrixCols[" << numCols << "] = {\n";
    for(int i = 0; i < numCols; ++i) {
      cout << cols[i] << ", \n";
    }
    cout << "};\n";
  }
  if (numRows == 0) {
    cout << "int *matrixRows = 0;\n";
  } else {
    cout << "int matrixRows[" << numRows << "] = {\n";
    for(int i = 0; i < numRows; ++i) {
      cout << rows[i] << ", \n";
    }
    cout << "};\n";
  }
  if (numVals == 0) {    
    cout << "double *matrixVals = 0;\n";
  } else {
    cout << "double matrixVals[" << numVals << "] = {\n";
    for(int i = 0; i < numVals; ++i) {
      if(vals[i] == 0) 
        cout << "0.0, \n";
      else 
        cout << vals[i] << ", \n";
    }
    cout << "};\n";      
  }
}  

bool MMElement::compare(const MMElement &elt1, const MMElement &elt2) {
  if (elt1.row < elt2.row) return true;
  else if (elt2.row < elt1.row) return false;
  else return elt1.col < elt2.col;
}

bool MMElement::compareCol(int limit){
  return (col < limit);
}

MMMatrix::MMMatrix(unsigned long n) {
  this->n = n;
}

MMMatrix::~MMMatrix() {

}

void MMMatrix::sort() {
  std::sort(elts.begin(), elts.end(), MMElement::compare);
}

void MMMatrix::add(int row, int col, double val) {
  if(val != 0)
    elts.push_back(MMElement(row, col, val));
}

void MMMatrix::normalize() {
  sort();
}
    
void MMMatrix::print() {
  cout << n << " " << n << " " << elts.size() << "\n";
  for (MMElement &elt : elts) {
    cout << elt.row << " "
         << elt.col << " "
         << elt.val << "\n";
  }
}

void MMMatrix::printMTX() {
  cout << n << " " << n << " " << elts.size() << "\n";
  for (MMElement &elt : elts) {
    cout << (elt.row + 1) << " "
         << (elt.col + 1) << " "
         << elt.val << "\n";
  }
}

Matrix* MMMatrix::toCSRMatrix() {
  long sz = elts.size();
  int *rows = new int[n+1];
  rows[0] = 0;
  int *cols = new int[sz]; 
  double *vals = new double[sz]; 

  unsigned int eltIndex = 0;
  unsigned int rowIndex = 0;
  for (auto &elt : elts) {
    cols[eltIndex] = elt.col;
    vals[eltIndex] = elt.val;
    while (rowIndex != elt.row) {
      rows[rowIndex+1] = eltIndex;
      rowIndex++;
    }
    eltIndex++;
  }
  while (rowIndex < n) {
    rows[rowIndex+1] = eltIndex;
    rowIndex++;
  }
  rows[n] = eltIndex;
  
  Matrix *matrix = new Matrix(rows, cols, vals, n, sz);
  matrix->numRows = n + 1;
  return matrix;
}

