//-------------------------------------------------------------------------
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University.
//-------------------------------------------------------------------------

// The extsort program.
// Jong Seong Lee
//
// Usage: ./extsort <input file name> <output file name>
//
//
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#define DEBUG 0
#define MINSIZE   10            // threshold for switching to bubblesort
#define INTERVAL  10
#define TAG_BUCKET 1001
#define TAG_BUCKET_SIZE 1002
#define TAG_WRITE 1003
#define PROC_MAIN 0

void debug_print(int rank, const char* msg) {
  if (DEBUG) {
    fprintf(stderr, "Proc %d %s\n", rank, msg);
  }
}

void debug_print_io(int rank, const char* msg, double time) {
  if (DEBUG) {
    fprintf(stderr, "Proc %d %s %f\n", rank, msg, time);
  }
}
// Swap two array elements
//
void swap(int *array, int i, int j) {
  if (i == j) return;
  int tmp = array[i];
  array[i] = array[j];
  array[j] = tmp;
}

// Bubble sort for the base cases
//
void bubblesort(int *array, int low, int high) {
  int i, j;

  if (low >= high)
    return;
  for (i = low; i <= high; i++)
    for (j = i+1; j <= high; j++)
      if (array[i] > array[j])
        swap(array, i, j);
}

// Pick an arbitrary element as pivot. Rearrange array
// elements into [smaller one, pivot, larger ones].
// Return pivot's index.
//
int partition(int *array, int low, int high) {
  int pivot = array[high];      // use highest element as pivot
  int middle = low;
  int i;
  for(i = low; i < high; i++)
    if(array[i] < pivot) {
      swap(array, i, middle);
      middle++;
    }
  swap(array, high, middle);
  return middle;
}

// QuickSort an array range
//
void quicksort(int *array, int low, int high) {
  if (high - low < MINSIZE) {
    bubblesort(array, low, high);
    return;
  }
  int middle = partition(array, low, high);
  if (low < middle)
    quicksort(array, low, middle-1);
  if (middle < high)
    quicksort(array, middle+1, high);
}

// Main routine for testing quicksort
//
int main(int argc, char **argv) {

  // check command line first
  if (argc < 3) {
    printf ("Usage: ./extsort <input file name> <output file name>\n");
    exit(0);
  }

  // common variables
  char* input_name = argv[1];
  char* output_name = argv[2];
  double time_all = MPI_Wtime();
  double time_io = 0;
  double time_temp = 0;
  int nprocs;
  int rank;
  MPI_File input;
  MPI_File output;
  MPI_Status st;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  if (nprocs < 2) {
    printf("Need at least 2 processes...\n");
    MPI_Finalize();
    return (1);
  }

  if (rank == PROC_MAIN) {
    debug_print(rank, "started...");

    int *buf;
    int *pivot = (int*)malloc((nprocs - 1) * sizeof(int));
    int *buckets[nprocs];
    // to count number of elements in a bucket
    int *count = (int*)malloc(nprocs * sizeof(int));
    int buf_count;
    int i;
    int j;
    MPI_Offset file_size;

    // record io time
    time_temp = MPI_Wtime();

    MPI_File_open(MPI_COMM_SELF, input_name, MPI_MODE_RDONLY, MPI_INFO_NULL, &input);

    MPI_File_get_size(input, &file_size);
    buf = (int*)malloc(file_size);
    buf_count = file_size / sizeof(int);
    MPI_File_read(input, buf, buf_count, MPI_INTEGER, &st);
    MPI_File_close(&input);

    // record io time
    time_io += MPI_Wtime() - time_temp;

    debug_print(rank, "reading input complete...");

    quicksort(buf, 0, nprocs * INTERVAL - 1);
    debug_print(rank, "quicksort complete...");

    // select elements at positions 10, 20, ... 10(P - 1) as pivots
    for (i = 0; i < nprocs - 1; i++) {
      pivot[i] = buf[INTERVAL * (i + 1)];
    }

    debug_print(rank, "pivot selected...");

    for (i = 0; i < nprocs; i++) {
      buckets[i] = (int*)malloc(buf_count * sizeof(int));
      count[i] = 0;
    }
    debug_print(rank, "bucket space allocated...");

    for (i = 0; i < buf_count; i++) {
      for (j = 0; j < nprocs; j++) {
        if (buf[i] < pivot[j]) {
          buckets[j][count[j]] = buf[i];
          count[j]++;
          break;
        }
      }
    }

    debug_print(rank, "elemetns distributed...");

    for (i = 1; i < nprocs; i++) {
      // send data to other processes
      MPI_Send(&count[i], 1, MPI_INT, i, TAG_BUCKET_SIZE, MPI_COMM_WORLD);
      MPI_Send(buckets[i], count[i], MPI_INT, i, TAG_BUCKET, MPI_COMM_WORLD);
    }

    debug_print(rank, "send complete...");

    quicksort(buckets[rank], 0, count[rank] - 1);
    debug_print(rank, "quicksort completed...");

    // record i/o time
    time_temp = MPI_Wtime();

    MPI_File_open(MPI_COMM_SELF, output_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &output);
    MPI_File_write(output, buckets[rank], count[rank], MPI_INT, &st);
    MPI_File_close(&output);

    // record i/o time
    time_io = MPI_Wtime() - time_temp;

    debug_print_io(rank, "sending i/o time:", time_io);
    MPI_Send(&time_io, 1, MPI_DOUBLE, 1, TAG_WRITE, MPI_COMM_WORLD);

    debug_print(rank, "file write completed...");

  } else {

    debug_print(rank, "started...");

    // getting count first to allocate memory
    int count;
    MPI_Recv(&count, 1, MPI_INT, 0, TAG_BUCKET_SIZE, MPI_COMM_WORLD, &st);
    int *bucket = (int*)malloc(count * sizeof(int));

    // receive bucket
    MPI_Recv(bucket, count, MPI_INT, 0, TAG_BUCKET, MPI_COMM_WORLD, &st);
    debug_print(rank, "bucket receive completed...");

    quicksort(bucket, 0, count - 1);
    debug_print(rank, "quicksort completed...");

    // receives the i/o time from former process to write in order
    debug_print_io(rank, "receiving i/o time: ", time_io);
    MPI_Recv(&time_io, 1, MPI_DOUBLE, rank - 1, TAG_WRITE, MPI_COMM_WORLD, &st);

    // record i/o time
    time_temp = MPI_Wtime();

    MPI_File_open(MPI_COMM_SELF, output_name, MPI_MODE_CREATE|MPI_MODE_WRONLY, MPI_INFO_NULL, &output);

    // append data
    MPI_File_seek(output, 0, MPI_SEEK_END);
    MPI_File_write(output, bucket, count, MPI_INT, &st);
    MPI_File_close(&output);

    // record i/o time
    time_io += MPI_Wtime() - time_temp;

    debug_print(rank, "write completed...");

    if (rank == nprocs - 1) {
      // the last one will print out timing
      double time_total = MPI_Wtime() - time_all;
      debug_print_io(rank, "has the final i/o time:", time_io);
      printf("Total elapsed time: %f\nWithout I/O: %f\n", time_total, time_total - time_io);
    } else {
      // sending the i/o to the next process for ordering
      debug_print_io(rank, "sending i/o time:", time_io);
      MPI_Send(&time_io, 1, MPI_DOUBLE, rank + 1, TAG_WRITE, MPI_COMM_WORLD);
    }
  }

  MPI_Finalize();
}
