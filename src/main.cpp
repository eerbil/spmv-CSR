#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <chrono>
#include "matrix.h"

using namespace std;


void spMV_CSR(double *vals, int *cols, int *rows, int N, double *v, double *w) {
  for (int i = 0; i < N; i++) {
    double sum = 0.0;
    for (int k = rows[i]; k < rows[i+1]; k++) {
      sum += vals[k] * v[cols[k]];
    }
    w[i] += sum;
  }
}

// Unrolling inner loop for 4 times
void spMV_CSR4(double *vals, int *cols, int *rows, int N, double *v, double *w) {
  for (int i = 0; i < N; i++) {
    double sum = 0.0;
    int k;
    for (k = rows[i]; k < rows[i+1]-4; k+=4) {
      sum += vals[k] * v[cols[k]];
      sum += vals[k+1] * v[cols[k+1]];
      sum += vals[k+2] * v[cols[k+2]];
      sum += vals[k+3] * v[cols[k+3]];
    }
    // Left-over elements
    for (; k < rows[i+1]; k++) {
      sum += vals[k] * v[cols[k]];
    }
    w[i] += sum;
  }
}

// Caller of this method is responsible for destructing the
// returned matrix. 
Matrix* readMatrixFromFile(string fileName) {
  ifstream mmFile(fileName.c_str());
  if (!mmFile.is_open()) {
    std::cerr << "Problem with file " << fileName << ".\n";
    exit(1);
  }
  string headerLine;
  // consume the comments until we reach the size info
  while (mmFile.good()) {
    getline (mmFile, headerLine);
    if (headerLine[0] != '%') break;
  }
  
  // Read N, M, NZ
  stringstream header(headerLine, ios_base::in);
  int n, m, nz;
  header >> n >> m >> nz;
  if (n != m) {
    std::cerr << "Only square matrices are accepted.\n";
    exit(1);
  }
  
  // Read rows, cols, vals
  MMMatrix matrix(n);
  int row; int col; double val;

  string line;
  for (int i = 0; i < nz; ++i) {
    getline(mmFile, line);
    stringstream linestream(line, ios_base::in);
    linestream >> row >> col;
    // pattern matrices do not contain val entry.
    // Such matrices are filled in with consecutive numbers.
    linestream >> val;
    if (linestream.fail()) 
      val = i+1;
    // adjust to zero index
    matrix.add(row-1, col-1, val);
  }
  mmFile.close();
  matrix.normalize();
  Matrix *csrMatrix = matrix.toCSRMatrix();
  return csrMatrix;
}

int main(int argc, const char *argv[]) {
  if (argc != 2) {
    std::cout << "Usage: spMV <matrixName>\n";
    exit(1);
  }
  
  string matrixName(argv[1]);
  Matrix *csrMatrix = readMatrixFromFile(matrixName + ".mtx");

  unsigned long N = csrMatrix->n;
  unsigned long NZ = csrMatrix->nz;
  double *v = new double[N];
  double *w = new double[N];

  // Artificially populate input and output vectors
  for(int i = 0; i < N; ++i) {
    w[i] = 0;
    v[i] = i;    
  }

  unsigned int ITERS;
  if (NZ < 5000) {
    ITERS = 500000;
  } else if (NZ < 10000) {
    ITERS = 200000;
  } else if (NZ < 50000) {
    ITERS = 100000;
  } else if (NZ < 100000) {
    ITERS = 50000;
  } else if (NZ < 200000) {
    ITERS = 10000;
  } else if (NZ < 1000000) {
    ITERS = 5000;
  } else if (NZ < 2000000) {
    ITERS = 1000;
  } else if (NZ < 3000000) {
    ITERS = 500;
  } else {
    ITERS = 200;
  }

  auto startTime = std::chrono::high_resolution_clock::now();
  for (int i=0; i < ITERS; i++) {
    spMV_CSR4(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
  }
  auto endTime = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
  auto runtime = (duration / (double)ITERS);

  std::cout << "Runtime (usecs): " << runtime << "\n";

    
  // For debugging purposes:
  //for(int i = 0; i < N; ++i) w[i] = 0;
  //spMV_CSR4(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
  //for(int i = 0; i < N; ++i) printf("%g\n", w[i]);
  
  return 0;
}

