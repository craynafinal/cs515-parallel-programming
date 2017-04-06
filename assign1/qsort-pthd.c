//-------------------------------------------------------------------------
// This is supporting software for CS415/515 Parallel Programming.
// Copyright (c) Portland State University.
//-------------------------------------------------------------------------

// qsort thread version
// Jong Seong Lee

#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <sched.h>
#include <pthread.h>

#include "task-queue.h"

#define MINSIZE   10            // threshold for switching to bubblesort

int *array= NULL,
    N = 0,
    numThreads = 1,
    complete = 0;

queue_t *queue;
pthread_mutex_t mutexQueue;
pthread_mutex_t mutexComplete;
pthread_cond_t conditionQueue;

void addTaskWrapper(int low, int high) {
  task_t *task = create_task(low, high);

  int ret;

  // shared variable is queue here
  pthread_mutex_lock(&mutexQueue);

  ret = add_task(queue, task);

  pthread_mutex_unlock(&mutexQueue);

  // wake up threads if there is a new task
  if (ret)
    pthread_cond_broadcast(&conditionQueue);
}

int isComplete() {
  int ret = 0;

  // shared variable is complete here
  pthread_mutex_lock(&mutexComplete);

  if (N == complete)
    ret = 1;

  pthread_mutex_unlock(&mutexComplete);

  return ret;
}

void increaseComplete(int count) {
  // complete is shared variable here
  pthread_mutex_lock(&mutexComplete);
  complete += count;
  pthread_mutex_unlock(&mutexComplete);
}

// Swap two array elements
//
void swap(int *array, int i, int j) {
  if (i == j) return;
  int tmp = array[i];
  array[i] = array[j];
  array[j] = tmp;
}

// Initialize array.
// - first generate [1,2,...,N]
// - then perform a random permutation
//
int *init_array(int N)  {
  int *array = (int *) malloc(sizeof(int) * N);
  int i;
  for (i = 0; i < N; i++) {
    array[i] = i + 1;
  }
  srand(time(NULL));
  for (i = 0; i < N; i++) {
    int j = (rand()*1./RAND_MAX) * (N-1);
    swap(array, i, j);
  }
  printf("Initialized array to a random permutation of [1..%d]\n", N);
  return array;
}

// Verify that array is sorted
//
void verify_array(int *array, int N) {
  int i;
  for (i = 0; i < N-1; i++)
    if (array[i]>array[i+1]) {
      printf("FAILED: array[%d]=%d, array[%d]=%d\n",
             i, array[i], i+1, array[i+1]);
      return;
    }
  printf("Result verified!\n");
}

// Bubble sort for the base cases
//
void bubblesort(int *array, int low, int high) {
  if (low >= high)
    return;
  int i, j;
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
    increaseComplete(high - low + 1);
    return;
  }
  int middle = partition(array, low, high);

  // the middle is sorted and done
  increaseComplete(1);

  // save it for later use
  addTaskWrapper(low, middle-1);

  if (middle < high)
    quicksort(array, middle+1, high);
}

void worker(long wid) {
  printf("Worker %ld started on %d\n", wid, sched_getcpu());

  while (!isComplete()) {

    // accessing queue
    pthread_mutex_lock(&mutexQueue);

    // wait here if there is no task in queue and it is not finishedd
    while (!isComplete() && queue->length == 0)
      pthread_cond_wait(&conditionQueue, &mutexQueue);

    task_t *task = remove_task(queue);

    pthread_mutex_unlock(&mutexQueue);

    // if task is legit, do quicksort
    if (task) {
      quicksort(array, task->low, task->high);
      free (task);
    }
  }
}

// Main routine for testing quicksort
//
int main(int argc, char **argv) {
  long k;

  // check command line first
  if (argc < 2) {
    // if there is no argument passed
    printf ("Usage: ./qsort-pthd <N> [<numThreads>]\n");
    exit(0);
  }

  if ((argc >= 3) && ((numThreads=atoi(argv[2])) < 1 )) {
    printf("<numThreads> must be greater than 0\n");
    exit(0);
  }
  if ((N = atoi(argv[1])) < 2) {
    printf ("<N> must be greater than 2\n");
    exit(0);
  }

  // checking the arg values
  printf("N: %d, numThreads: %d\n", N, numThreads);

  // initialize
  pthread_t thread[numThreads];
  pthread_mutex_init(&mutexQueue, NULL);
  pthread_mutex_init(&mutexComplete, NULL);
  pthread_cond_init(&conditionQueue, NULL);

  array = init_array(N);

  task_t *startTask = create_task(0, N - 1);
  queue = init_queue(0);
  add_task(queue, startTask);

  // create worker threads
  for (k = 0; k < numThreads-1; k++)
    pthread_create(&thread[k], NULL, (void*)worker, (void*)k);

  // the main thread runs worker routine too
  worker(numThreads-1);

  // join
  for (k = 0; k < numThreads-1; k++)
    pthread_join(thread[k], NULL);

  printf("Sorting started ...\n");
  quicksort(array, 0, N-1);
  printf("... completed.\n");

  verify_array(array, N);
}