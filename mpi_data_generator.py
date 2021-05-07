import subprocess
from subprocess import Popen, PIPE
import os
import csv

runs = 10
processors = [1,2,4,8,16,32]
table = []

for j in range(len(processors)):
    data = [0,0,0,0,0,0]
    print("Processors: ", processors[j])
    for i in range(runs):
        process = Popen("/apps/intel/compilers_and_libraries_2018.2.199/linux/mpi/intel64/bin/mpiexec -machinefile machines -np {p} ./a.out".format(p = processors[j]), stdout=PIPE, shell = True)
        out = process.stdout.read()
        paths = []
        for line in out.split('\n'):
            paths.append(line)
        paths.pop()
        for k in range(6):
            data[k] += float(paths[k])
    for k in range(6):
        data[k] /= runs
    table.append(data)
    print("Creating frequency table: ", data[0])
    print("Creating Huffman Tree: ", data[1])
    print("Creating hashtable of encodes: ", data[2])
    print("Compression: ", data[3])
    print("Decompression: ", data[4])
    print("Total Time: ", data[5])
    print("----------------------------------------------------")

with open("data.csv","w+") as my_csv:
    csvWriter = csv.writer(my_csv,delimiter=',')
    csvWriter.writerows(table)
    
    