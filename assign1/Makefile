#----------------------------------------------------------------------
# CS415/515 Assignment 1
#----------------------------------------------------------------------

FLAGS = -std=c99 -pthread -g -DVERBOSE
CC = gcc

%: task-queue.h task-queue.c %.c 
	$(CC) $(FLAGS) -o $@ task-queue.c $@.c
