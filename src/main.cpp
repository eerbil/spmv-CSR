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

// Unrolling inner loop for 5 times
void spMV_CSR5(double *vals, int *cols, int *rows, int N, double *v, double *w) {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        int k;
        for (k = rows[i]; k < rows[i+1]-5; k+=5) {
            sum += vals[k] * v[cols[k]];
            sum += vals[k+1] * v[cols[k+1]];
            sum += vals[k+2] * v[cols[k+2]];
            sum += vals[k+3] * v[cols[k+3]];
            sum += vals[k+4] * v[cols[k+4]];
        }
        // Left-over elements
        for (; k < rows[i+1]; k++) {
            sum += vals[k] * v[cols[k]];
        }
        w[i] += sum;
    }
}

// Unrolling inner loop for 6 times
void spMV_CSR6(double *vals, int *cols, int *rows, int N, double *v, double *w) {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        int k;
        for (k = rows[i]; k < rows[i+1]-6; k+=6) {
            sum += vals[k] * v[cols[k]];
            sum += vals[k+1] * v[cols[k+1]];
            sum += vals[k+2] * v[cols[k+2]];
            sum += vals[k+3] * v[cols[k+3]];
            sum += vals[k+4] * v[cols[k+4]];
            sum += vals[k+5] * v[cols[k+5]];
        }
        // Left-over elements
        for (; k < rows[i+1]; k++) {
            sum += vals[k] * v[cols[k]];
        }
        w[i] += sum;
    }
}

// Unrolling inner loop for 8 times
void spMV_CSR8(double *vals, int *cols, int *rows, int N, double *v, double *w) {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        int k;
        for (k = rows[i]; k < rows[i+1]-8; k+=8) {
            sum += vals[k] * v[cols[k]];
            sum += vals[k+1] * v[cols[k+1]];
            sum += vals[k+2] * v[cols[k+2]];
            sum += vals[k+3] * v[cols[k+3]];
            sum += vals[k+4] * v[cols[k+4]];
            sum += vals[k+5] * v[cols[k+5]];
            sum += vals[k+6] * v[cols[k+6]];
            sum += vals[k+7] * v[cols[k+7]];
        }
        // Left-over elements
        for (; k < rows[i+1]; k++) {
            sum += vals[k] * v[cols[k]];
        }
        w[i] += sum;
    }
}

// Unrolling inner loop for 10 times
void spMV_CSR10(double *vals, int *cols, int *rows, int N, double *v, double *w) {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        int k;
        for (k = rows[i]; k < rows[i+1]-10; k+=10) {
            sum += vals[k] * v[cols[k]];
            sum += vals[k+1] * v[cols[k+1]];
            sum += vals[k+2] * v[cols[k+2]];
            sum += vals[k+3] * v[cols[k+3]];
            sum += vals[k+4] * v[cols[k+4]];
            sum += vals[k+5] * v[cols[k+5]];
            sum += vals[k+6] * v[cols[k+6]];
            sum += vals[k+7] * v[cols[k+7]];
            sum += vals[k+8] * v[cols[k+8]];
            sum += vals[k+9] * v[cols[k+9]];
        }
        // Left-over elements
        for (; k < rows[i+1]; k++) {
            sum += vals[k] * v[cols[k]];
        }
        w[i] += sum;
    }
}

// Unrolling inner loop for 12 times
void spMV_CSR12(double *vals, int *cols, int *rows, int N, double *v, double *w) {
    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        int k;
        for (k = rows[i]; k < rows[i+1]-12; k+=12) {
            sum += vals[k] * v[cols[k]];
            sum += vals[k+1] * v[cols[k+1]];
            sum += vals[k+2] * v[cols[k+2]];
            sum += vals[k+3] * v[cols[k+3]];
            sum += vals[k+4] * v[cols[k+4]];
            sum += vals[k+5] * v[cols[k+5]];
            sum += vals[k+6] * v[cols[k+6]];
            sum += vals[k+7] * v[cols[k+7]];
            sum += vals[k+8] * v[cols[k+8]];
            sum += vals[k+9] * v[cols[k+9]];
            sum += vals[k+10] * v[cols[k+10]];
            sum += vals[k+11] * v[cols[k+11]];
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
  if (argc != 3) {
    std::cout << "Usage: spMV <matrixName>\n";
    exit(1);
  }
  
  string matrixName(argv[1]);
  Matrix *csrMatrix = readMatrixFromFile(matrixName + ".mtx");

  int chosenMethod = atoi(argv[2]);
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
    if(chosenMethod == 0) {
        for (int i=0; i < ITERS; i++)
        {
            spMV_CSR(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        auto runtime = (duration / (double)ITERS);
        //  std::cout << "Runtime (usecs): " << runtime << "\n";
        std::cout << runtime << "\n";
    }
    else if(chosenMethod == 4) {
        for (int i=0; i < ITERS; i++)
        {
            spMV_CSR4(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        auto runtime = (duration / (double)ITERS);
        //  std::cout << "Runtime (usecs): " << runtime << "\n";
        std::cout << runtime << "\n";
    }
    else if(chosenMethod == 5) {
        for (int i=0; i < ITERS; i++)
        {
            spMV_CSR5(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        auto runtime = (duration / (double)ITERS);
        //  std::cout << "Runtime (usecs): " << runtime << "\n";
        std::cout << runtime << "\n";
    }
    else if(chosenMethod == 6) {
        for (int i=0; i < ITERS; i++)
        {
            spMV_CSR6(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        auto runtime = (duration / (double)ITERS);
        //  std::cout << "Runtime (usecs): " << runtime << "\n";
        std::cout << runtime << "\n";
    }
    else if(chosenMethod == 8) {
        for (int i=0; i < ITERS; i++)
        {
            spMV_CSR8(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        auto runtime = (duration / (double)ITERS);
        //  std::cout << "Runtime (usecs): " << runtime << "\n";
        std::cout << runtime << "\n";
    }
    else if(chosenMethod == 10) {
        for (int i=0; i < ITERS; i++)
        {
            spMV_CSR10(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        auto runtime = (duration / (double)ITERS);
        //  std::cout << "Runtime (usecs): " << runtime << "\n";
        std::cout << runtime << "\n";
    }
    else if(chosenMethod == 12) {
        for (int i=0; i < ITERS; i++)
        {
            spMV_CSR12(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
        }
        auto endTime = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();
        auto runtime = (duration / (double)ITERS);
        //  std::cout << "Runtime (usecs): " << runtime << "\n";
        std::cout << runtime << "\n";
    } else {
        std::cout << "Method not found" << "\n";
        exit(1);
    }
  // For debugging purposes:
  //for(int i = 0; i < N; ++i) w[i] = 0;
  //spMV_CSR4(csrMatrix->vals, csrMatrix->cols, csrMatrix->rows, csrMatrix->n, v, w);
  //for(int i = 0; i < N; ++i) printf("%g\n", w[i]);
  
  return 0;
}

