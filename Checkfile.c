#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include <memory.h>

int main(int argc, char* argv[]){
  long int n = atol(argv[1]);
  char* inFilename = argv[2];
  char* inFilename2 = argv[3];
  
  int rank, size;

  //initial the MPI_Comm
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // open input file with MPI-IO
  MPI_File infh, infh2;
  MPI_File_open(MPI_COMM_WORLD, inFilename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infh);
  MPI_File_open(MPI_COMM_WORLD, inFilename2, MPI_MODE_RDONLY, MPI_INFO_NULL, &infh2);

  float* buffer_1 = (float*) malloc (sizeof(float)*(n));
  MPI_File_read(infh, buffer_1, n, MPI_FLOAT, &status);
  float* buffer_2 = (float*) malloc (sizeof(float)*(n));
  MPI_File_read(infh2, buffer_2, n, MPI_FLOAT, &status);

  for(int i=0; i<n; i++){
    if(buffer_1[i]!=buffer_2[i]) printf("Element %d, Test Buffer Element = %f, Sorted Buffer Element = %f\n", i, buffer_1[i], buffer_2[i]);
  }

  MPI_File_close(&infh);
  MPI_File_close(&infh2);
   
  // termination
  MPI_Finalize();
  return 0;
}           
