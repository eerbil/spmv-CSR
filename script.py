import os
import subprocess
import csv

methods = [0, 4, 5, 6, 8, 10, 12]

matrices = os.listdir("matrices_gpce")
print("Matrix List: ")
print(matrices)
count = 0

for matrix in matrices:
    matrices[matrices.index(matrix)] = matrix.replace(".mtx", "")
    count += 1
print(matrices)
print(count)

for k in range(0, 7):
    finalList = []
    for i in range(0, count):
        print("Running: " + matrices[i] + " with csr unrolled " + str(methods[k]) + " times")
        runtimeList = [0, 0, 0]
        for j in range(0, 3):
            result = subprocess.check_output(["/Users/eliferbil/Documents/SpMV/spmv-CSR/build/spMV",
                     "/Users/eliferbil/Documents/SpMV/matrices_gpce/" + matrices[i], str(methods[k])])
            runtimeList[j] = float(result)
        print(runtimeList)
        runtimeList.sort()
        finalResult = runtimeList[0]
        print(finalResult)
        finalList.append([matrices[i], finalResult])
    with open('spmv-csr'+str(methods[k])+'.csv', 'w') as f:
        writer = csv.writer(f)
        writer.writerows(finalList)
