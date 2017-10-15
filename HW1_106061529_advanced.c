#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include <memory.h>

int compare(const void *a, const void *b){
  float c = *(float *)a;
  float d = *(float *)b;
  if(c < d) {return -1;}               
  else if (c == d) {return 0;}      
  else return 1;                         
}

void mergesort(float* unsorted_arr, int start_a, int start_b, int start_c){
  int start_loc = start_a;
  int end_a = start_b;
  int end_b = start_c;
  int total_length = end_b - start_a;
  float* slice_sort_arr = (float*) malloc (sizeof(float)*(total_length));
  int output_arr_id = 0;
  while(start_a < end_a && start_b < end_b){
    if(unsorted_arr[start_a] < unsorted_arr[start_b]){
      slice_sort_arr[output_arr_id] = unsorted_arr[start_a]; 
      start_a++;
    }
    else{
      slice_sort_arr[output_arr_id] = unsorted_arr[start_b]; 
      start_b++;
    }
    output_arr_id++;
  }
  if(start_a < end_a){
    memcpy(&slice_sort_arr[output_arr_id], &unsorted_arr[start_a], sizeof(float) * (end_a - start_a));
  }
  else{
    memcpy(&slice_sort_arr[output_arr_id], &unsorted_arr[start_b], sizeof(float) * (end_b - start_b));
  }
  memcpy(&unsorted_arr[start_loc], slice_sort_arr, sizeof(float) * total_length);
  free(slice_sort_arr);
}


int main(int argc, char* argv[]){
  long int n = atol(argv[1]);
  char* inFilename = argv[2];
  char* outFilename = argv[3];
  
  int rank, size;

  //initial the MPI_Comm
  MPI_Status status;
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // open input file with MPI-IO
  MPI_File infh, outfh;
  MPI_File_open(MPI_COMM_WORLD, inFilename, MPI_MODE_RDONLY, MPI_INFO_NULL, &infh);
  MPI_File_open(MPI_COMM_WORLD, outFilename, MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &outfh);
  
  //read file from MPI_FILE
  int item_per_prc = n/size;
  int use_size;
  if(item_per_prc < 1){
    item_per_prc = 1;
    use_size = n;
  }
  else{
    use_size = size;
  }
  if(rank < use_size){
    int start_id = rank * item_per_prc;
    int end_id, num_item;
    if(rank < use_size-1){
      num_item = item_per_prc;
      end_id = start_id + item_per_prc;
    }
    else{
      end_id = n; 
      num_item = n - start_id;
    }
    if(rank == use_size-1){
      float* recv_buffer = (float*) malloc (sizeof(float)*(n));
      MPI_File_read_at(infh, start_id*sizeof(float), &recv_buffer[start_id], num_item, MPI_FLOAT, &status);
      qsort(&recv_buffer[start_id], num_item, sizeof(float), compare);
      int* partition_arr = (int*) malloc (sizeof(int)*(use_size));
      partition_arr[use_size-1] = start_id;
      int recv_id = start_id;
      for(int remain_recv_time=rank; remain_recv_time>0; remain_recv_time--){
        recv_id -= item_per_prc;
        partition_arr[remain_recv_time-1] = recv_id;
        MPI_Recv(&recv_buffer[recv_id], item_per_prc, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
      }
      int numPartitions = use_size;
      int new_numPartitions, p;
      while(numPartitions > 1){
        p = 0;
        for(; p < numPartitions-1; p+=2){
          if(p+2 == numPartitions){
            mergesort(recv_buffer, partition_arr[p], partition_arr[p+1], n);
          }
          else{
            mergesort(recv_buffer, partition_arr[p], partition_arr[p+1], partition_arr[p+2]);
          }
          partition_arr[p/2] = partition_arr[p];
        }
        if(p == numPartitions-1){
          new_numPartitions = (numPartitions/2) + 1;
          partition_arr[new_numPartitions-1] = partition_arr[numPartitions-1];
          numPartitions = new_numPartitions;
        }
        else{
          numPartitions = numPartitions/2; 
        }
      }
      MPI_File_write(outfh, recv_buffer, n, MPI_FLOAT, &status);
    }
    else{
      float* pass_buffer = (float*) malloc (sizeof(float)*(item_per_prc));
      MPI_File_read_at(infh, start_id*sizeof(float), pass_buffer, num_item, MPI_FLOAT, &status);
      qsort(pass_buffer, num_item, sizeof(float), compare);
      if(rank < use_size-1){
        MPI_Send(pass_buffer, item_per_prc, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
        for(int remain_pass_time=rank; remain_pass_time>0; remain_pass_time--){
          MPI_Recv(pass_buffer, item_per_prc, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
          MPI_Send(pass_buffer, item_per_prc, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
        }
      }
    }
  }
  MPI_File_close(&outfh);
  MPI_File_close(&infh);
   
  // termination
  MPI_Finalize();
  return 0;
}           
