Commands to compile and run codes are given at the starting of the code file itself, along with some informations about output as well.

Note: For accuracy purposes, we did the experiments 10 or 20 times and took average over it for all our time measurements.

Serial Code file name:
1. serial_code.cpp

OpenMP parallel code file names:
1. openmp_version1.cpp
2. openmp_version2.cpp
These two files are two versions of parallel codes, using OpenMP, as discussed in the report.

MPI parallel code:
1. mpi_parallel.cpp
2. mpi_data_gen_code.cpp
Use .bashrc and machines file (provided with codes) when running this code.
In MPI, to take average over many experiments we had to write a python script mpi_data_generator.py (which is also submitted).

Steps for mpi data generation:
1. compile mpi_data_gen_code.cpp
2. run mpi_data_generator.py script using "python mpi_data_generator.py"
3. Data will be printed and one file data.csv will also be generated.

Each and every code will print timings taken by all five steps in our algorithm, compressed file size and total time.
Where total time = sum of timings taken by all five steps + time taken to read the original file into buffer + time to write the decompressed file

Each code will generate one file which is the output of our code i.e. decompressed file, which is must be equal to the original file. 
We can check whether decompressed file is same as original file by using "diff command" in linux.

We have taken datasets from the website : https://klacansky.com/open-scivis-datasets/

Names of datasets we have generated results for:
Bonsai - 16 MB
Anuerism - 16 MB
Pancreas - 120 MB
Head Anuerism - 256 MB
Magnetic Reconnection Simulation - 512 MB

Some other datasets we used:
Neurons in Marmoset Visual Cortex - 314 MB

Note that the dataset size must be less than 512 MB for our code to work, since we have used (int) data type everywhere in our code.
This limitation can be improved by using long long or even (big int) type of data types.