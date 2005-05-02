
#include <stdlib.h>
#include "taucs.h"

#undef malloc
#undef calloc
#undef realloc
#undef free

void* taucs_malloc_stub (size_t size)               { return malloc(size); }
void* taucs_calloc_stub (size_t nmemb, size_t size) { return calloc(nmemb,size); }
void* taucs_realloc_stub(void* ptr, size_t size)    { return realloc(ptr,size); }
void  taucs_free_stub   (void* ptr)                 { free(ptr); }

#if !defined(TAUCS_MEMORY_TEST_yes)

double taucs_allocation_amount()   { return 0.0; }
int    taucs_allocation_count()    { return 0; }
int    taucs_allocation_attempts() { return 0; }
void   taucs_allocation_assert_clean() {}
void   taucs_allocation_mark_clean() {}
void   taucs_allocation_induce_failure(int i) {}

#else /* we do want memory testing */

#define TABLE_SIZE 100000

static int    allocation_initialized = 0;

static void*  allocation_ptr  [TABLE_SIZE];
static double allocation_size [TABLE_SIZE];
static char*  allocation_file [TABLE_SIZE];
static int    allocation_line [TABLE_SIZE];
static int    allocation_clean[TABLE_SIZE];

static int    allocation_attempts;
static int    allocation_count;
static double allocation_amount;
static int    allocation_clean_count;
static double allocation_clean_amount;

static int    allocation_induced_failure;

static void allocation_init()
{
  int i;

  allocation_initialized = 1;

  allocation_induced_failure = -1;
  allocation_attempts = 0;
  allocation_count = 0;
  allocation_amount = 0.0;

  for (i=0; i<TABLE_SIZE; i++)
    allocation_ptr[i] = NULL;
}

static void allocation_insert(void* ptr, double size, char* file, int line)
{
  int i,j,slot;
  union {
    void* p;
    int   i;
  } converter;

  slot = -1; /* none found yet */

  converter.p = ptr;
  j = converter.i % TABLE_SIZE;

  for (i=0; i<TABLE_SIZE; i++) {
    if (allocation_ptr[j] == NULL) {
      slot = j;
      break;
    }
    j++;
  }

  if (slot == -1) {
    taucs_printf("TAUCS ALLOCATION ERROR: ALLOCATION TABLE IS FULL\n");
    exit(1);
  }

  allocation_ptr  [slot] = ptr;
  allocation_size [slot] = size;
  allocation_file [slot] = file;
  allocation_line [slot] = line;
  allocation_clean[slot] = 0;
}

static double allocation_delete(void* ptr)
{
  int i,j,slot;
  union {
    void* p;
    int   i;
  } converter;
  double size;

  slot = -1; /* none found yet */

  converter.p = ptr;
  j = converter.i % TABLE_SIZE;

  for (i=0; i<TABLE_SIZE; i++) {
    if (allocation_ptr[j] == ptr) {
      slot = j;
      break;
    }
    j++;
  }

  if (slot == -1) return -1.0;

  size = allocation_size[slot];

  allocation_ptr[slot] = NULL;

  return size;
}

double taucs_allocation_amount()   { return allocation_amount; }
int    taucs_allocation_count()    { return allocation_count; }
int    taucs_allocation_attempts() { return allocation_attempts; }

void   taucs_allocation_induce_failure(int i)  { allocation_induced_failure = i; }

void taucs_allocation_assert_clean() 
{
  int i,clean = 1;
 
  for (i=0; i<TABLE_SIZE; i++) {
    if (allocation_ptr[i] != NULL && !allocation_clean[i]) {
      clean = 0;
    }
  }

  if (!clean) {
    taucs_printf("TAUCS ALLOCATION ERROR: ASSERTED CLEAN, BUT FOUND\n");

    for (i=0; i<TABLE_SIZE; i++) {
      if (allocation_ptr[i] != NULL && !allocation_clean[i]) {
	taucs_printf("\tBLOCK ALLOCATED AT %s:%d STILL ALLOCATED (%.2e BYTES)\n",
		     allocation_file[i],allocation_line[i],allocation_size[i]);
      }
    }
    exit(1);
  }
}

void taucs_allocation_mark_clean() 
{
  int i;

  allocation_attempts = 0;
  allocation_clean_count  = allocation_count;
  allocation_clean_amount = allocation_amount;

  for (i=0; i<TABLE_SIZE; i++) {
    if (allocation_ptr[i] != NULL)
      allocation_clean[i] = 1;
  }
}

void* taucs_internal_calloc(size_t nmemb, size_t size,
			    char* file, int line)
{
  void* ptr;

  if (nmemb*size == 0) {
    taucs_printf("TAUCS ALLOCATION: ZERO SIZE (%s:%d)\n",
		 file, line);
    return NULL;
  }

  if (!allocation_initialized) allocation_init();

  if (allocation_induced_failure == allocation_attempts) {
    allocation_induced_failure = -1;
    taucs_printf("TAUCS ALLOCATION: INDUCING FAILURE (%s:%d, count=%d)\n",
		 file, line, allocation_attempts);
    return NULL;
  }

  ptr = calloc(nmemb,size);
  
  if (ptr) {
    allocation_count++;
    allocation_attempts++;
    allocation_amount += (double) nmemb * (double) size;
    allocation_insert(ptr, (double) nmemb * (double) size, file, line);
  } else 
    taucs_printf("TAUCS ALLOCATION WARNING: CALLOC AT %s:%d FAILED\n",
		 file,line);

  return ptr;
}

void* taucs_internal_malloc(size_t size,
			    char* file, int line)
{
  void* ptr;

  if (size == 0) {
    taucs_printf("TAUCS ALLOCATION: ZERO SIZE (%s:%d)\n",
		 file, line);
    return NULL;
  }

  if (!allocation_initialized) allocation_init();

  if (allocation_induced_failure == allocation_attempts) {
    allocation_induced_failure = -1;
    taucs_printf("TAUCS ALLOCATION: INDUCING FAILURE (%s:%d, count=%d)\n",
		 file, line, allocation_attempts);
    return NULL;
  }

  ptr = malloc(size);

  if (ptr) {
    allocation_count++;
    allocation_attempts++;
    allocation_amount += (double) size;
    allocation_insert(ptr, (double) size, file, line);
  } else 
    taucs_printf("TAUCS ALLOCATION WARNING: CALLOC AT %s:%d FAILED\n",
		 file,line);

  return ptr;
}

void* taucs_internal_realloc(void *oldptr, size_t size,
			    char* file, int line)
     
{
  void* ptr;

  if (size == 0) {
    taucs_printf("TAUCS ALLOCATION: ZERO SIZE (%s:%d)\n",
		 file, line);
    return NULL;
  }

  if (!allocation_initialized) allocation_init();


  if (allocation_induced_failure == allocation_attempts) {
    allocation_induced_failure = -1;
    taucs_printf("TAUCS ALLOCATION: INDUCING FAILURE (%s:%d, count=%d)\n",
		 file, line, allocation_attempts);
    return NULL;
  }

  ptr= realloc(oldptr,size);

  /* if realloc returns NULL, nothing happened (memory is not freed) */

  if (ptr) {
    double oldsize;

    oldsize = allocation_delete(oldptr);
    if (oldsize == -1.0) {
      taucs_printf("TAUCS ALLOCATION ERROR: REALLOC AT %s:%d NOT ALLOCATED\n",
		   file,line);
      exit(1);
    }
    allocation_count--;
    allocation_amount -= oldsize;

    allocation_count++;
    allocation_attempts++;
    allocation_amount += (double) size;
    allocation_insert(ptr, (double) size, file, line);
  } else 
    taucs_printf("TAUCS ALLOCATION WARNING: REALLOC AT %s:%d FAILED\n",
		 file,line);

  return ptr;
}

void taucs_internal_free(void *oldptr,
			 char* file, int line)
{
  double oldsize;

  if (!allocation_initialized) allocation_init();

  if (!oldptr) return;

  oldsize = allocation_delete(oldptr);
  if (oldsize == -1.0) {
    taucs_printf("TAUCS ALLOCATION ERROR: FREE AT %s:%d NOT ALLOCATED\n",
		 file,line);
    exit(1);
  }
  allocation_count--;
  allocation_amount -= oldsize;

  free(oldptr);
}

#endif /* TAUCS_MEMORY_TEST_yes */

