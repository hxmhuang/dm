#include <stdio.h>
#include <stdlib.h>
#include "pthread.h"


int pthread_yield(void);

/**
 * The following names are published to Fortran source
 */
#define CRETRD cretrd_
#define YLDTRD yldtrd_
#define SLFTRD slftrd_
#define JONTRD jontrd_
#define EXITRD exitrd_
#define LOKTRD loktrd_
#define ULKTRD ulktrd_
#define GETMTX getmtx_
#define RELMTX relmtx_

/**
 * CRETRD() ­ create a new fortran thread through a subroutine.
 *
 * @param thread_func: a function pointer, pointing to a subroutine
 * @param theThread: the thread ID
 */
void CRETRD(void *(*thread_func)(void *), pthread_t *theThread) {
  pthread_create(theThread, NULL, thread_func, NULL);
} 

/**
 * YLDTRD() ­ yeild control to othre threads
 */
void YLDTRD() {
  pthread_yield();
}

/**
 * SLFTRD() ­ get the thread ID
 */
pthread_t SLFTRD() {
  return pthread_self();
}

/**
 * LOKTRD() ­ locks the execution of all threads till we have
 *            the mutex
 */
void LOKTRD(pthread_mutex_t **theMutex) {
  pthread_mutex_lock(*theMutex);
}

/**
 * ULKTRD() ­ unlocks the execution of all threads that were
 *            stopped due to this mutex
 */
void ULKTRD(pthread_mutex_t **theMutex) {
  pthread_mutex_unlock(*theMutex);
}

/**
 * GETMTX() ­ get a new mutex object
 */
void GETMTX(pthread_mutex_t **theMutex) {
  *theMutex = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t));
  pthread_mutex_init(*theMutex, NULL);
}

/**
 * RELMTX() ­ release a mutex object
 */
void RELMTX(pthread_mutex_t **theMutex) {
  pthread_mutex_destroy(*theMutex);
  free(*theMutex);
}

/**
 * JONTRD() ­ waits for thread ID to join
 */
void JONTRD(pthread_t *theThread) {
  int value = 0;
  pthread_join(*theThread, (void **)&value);
}

/**
 * EXITRD() ­ exit from a thread
 */
void EXITRD(void *status) {
  pthread_exit(status);
} /* end of function EXITRD() */
