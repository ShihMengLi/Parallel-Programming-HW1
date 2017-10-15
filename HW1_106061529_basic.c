#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <mpi.h>
#include <memory.h>

bool oe_sort_inside(float* array, int start_posi, int length){
  float swap;
  bool sort_down = true;
  for(int i = start_posi; i<length-1; i+=2){
    if(array[i] > array[i+1]){
      swap = array[i];
      array[i] = array[i+1];
      array[i+1] = swap;
      sort_down = false;
    }
  }
  return sort_down;
}

bool oe_sort_outside(float external, float* array, int index, bool head_flag){
  if(head_flag && (array[index] < external)) {array[index] = external; return false;}
  if(!head_flag && (array[index] > external)) {array[index] = external; return false;}
  return true;
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

  int item_per_prc = n/size;
  int use_size;
  if(item_per_prc < 1){
    item_per_prc = 1;
    use_size = n;
  }
  else{
    use_size = size;
  }
  bool sorted = false, sorted_tmp = false;
  if(rank < use_size){
    float external;
    int start_id = rank * item_per_prc;
    int end_id, num_item, head, tail;
    bool item_is_even = false, item_is_odd = false, prc_is_even = false, prc_is_odd = false, head_prc = false, tail_prc = false;
    if (rank == 0) head_prc = true; 
    if (rank == use_size-1) tail_prc = true;
    bool head_flag;
    if(item_per_prc%2 == 0) item_is_even = true;
    else item_is_odd = true;
    if(rank%2 == 0) prc_is_even = true;
    else prc_is_odd = true;
    if(rank < use_size-1){
      num_item = item_per_prc;
      end_id = start_id + item_per_prc;
    }
    else{
      end_id = n; 
      num_item = n - start_id;
    }
    float* buffer = (float*) malloc (sizeof(float)*(num_item));
    MPI_File_read_at(infh, start_id*sizeof(float), buffer, num_item, MPI_FLOAT, &status);
    head = 0;
    tail = num_item - 1;
    while(!sorted){
      sorted = true;
      if(item_is_odd){
        //Even phase
        //Sort in-side
        if(prc_is_even) sorted = oe_sort_inside(buffer, 0, num_item) & sorted;
        if(prc_is_odd) sorted = oe_sort_inside(buffer, 1, num_item) & sorted;

        //Sort out-side   Even Up, Odd down 
        if(prc_is_even && !tail_prc){                                            //Send up
          head_flag = false;
          MPI_Send(&buffer[tail], 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
          MPI_Recv(&external, 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status);
          sorted = oe_sort_outside(external, buffer, tail, head_flag) & sorted;
        }  
        if(prc_is_odd && !head_prc){                                             //Send down
          head_flag = true;
          MPI_Recv(&external, 1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
          MPI_Send(&buffer[head], 1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);
          sorted = oe_sort_outside(external, buffer, head, head_flag) & sorted;
        }  
        //Odd Phase
        //Sort in-side
        if(prc_is_even) sorted = oe_sort_inside(buffer, 1, num_item) & sorted;
        if(prc_is_odd) sorted = oe_sort_inside(buffer, 0, num_item) & sorted;

        //Sort out-side   Even down, Odd up
        if(prc_is_even && !head_prc){                                             //Send down
          head_flag = true;
          MPI_Send(&buffer[head], 1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);
          MPI_Recv(&external, 1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
          sorted = oe_sort_outside(external, buffer, head, head_flag) & sorted;
        }  
        if(prc_is_odd && !tail_prc){                                            //Send up
          head_flag = false;
          MPI_Recv(&external, 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status);
          MPI_Send(&buffer[tail], 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
          sorted = oe_sort_outside(external, buffer, tail, head_flag) & sorted;
        }
      }
      if(item_is_even){
        //Even phase
        //Sort in-side
        sorted = oe_sort_inside(buffer, 0, num_item) & sorted;
        //Odd phase
        //Sort in-side
        sorted = oe_sort_inside(buffer, 1, num_item) & sorted;
        //Sort out-side
        //Round 1 Even Up, Odd down
        if(prc_is_even && !tail_prc){                                            //Send up
          head_flag = false;
          MPI_Send(&buffer[tail], 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
          MPI_Recv(&external, 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status);
          sorted = oe_sort_outside(external, buffer, tail, head_flag) & sorted;
        }  
        if(prc_is_odd && !head_prc){                                             //Send down
          head_flag = true;
          MPI_Recv(&external, 1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
          MPI_Send(&buffer[head], 1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);
          sorted = oe_sort_outside(external, buffer, head, head_flag) & sorted;
        } 
        //Round 2 Even down, Odd up
        if(prc_is_even && !head_prc){                                             //Send down
          head_flag = true;
          MPI_Send(&buffer[head], 1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD);
          MPI_Recv(&external, 1, MPI_FLOAT, rank-1, 0, MPI_COMM_WORLD, &status);
          sorted = oe_sort_outside(external, buffer, head, head_flag) & sorted;
        }  
        if(prc_is_odd && !tail_prc){                                            //Send up
          head_flag = false;
          MPI_Recv(&external, 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD, &status);
          MPI_Send(&buffer[tail], 1, MPI_FLOAT, rank+1, 0, MPI_COMM_WORLD);
          sorted = oe_sort_outside(external, buffer, tail, head_flag) & sorted;
        }
      }
      MPI_Allreduce(&sorted, &sorted_tmp, 1, MPI_CHAR, MPI_BAND, MPI_COMM_WORLD);
      sorted = sorted_tmp;
    }
    MPI_File_write_at(outfh, start_id*sizeof(float), buffer, num_item, MPI_FLOAT, &status);
  }
  else{
    while(!sorted){
      sorted = true;
      MPI_Allreduce(&sorted, &sorted_tmp, 1, MPI_CHAR, MPI_BAND, MPI_COMM_WORLD);
      sorted = sorted_tmp;
    }
  }
  MPI_File_close(&outfh);
  MPI_File_close(&infh);
   
  // termination
  MPI_Finalize();
  return 0;
}           
