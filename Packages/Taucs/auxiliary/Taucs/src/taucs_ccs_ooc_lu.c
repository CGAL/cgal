/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*                                                       */
/*********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <assert.h>

#include "taucs.h"

#define HEADER_NROWS   0
#define HEADER_NCOLS   1
#define HEADER_FLAGS   2
#define HEADER_COLPERM 3
#define HEADER_IPIVOTS 4
#define HEADER_LCLEN   5
#define HEADER_UCLEN   6
#define HEADER 7

#ifndef TAUCS_CORE_GENERAL

#ifdef OSTYPE_win32
#include <io.h>
#else
#include <unistd.h>
#include <sys/uio.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <errno.h>

#define iabs(x) ((x) > 0 ? (x) : (-(x)))

static double remaining_memory;

/*********************************************************/
/* NEW IO ROUTINES                                       */
/*********************************************************/

/*
  Layout of data structures in the file:
  m       (number of rows)
  n       (number of cols)
  flags
  colperm
  ipivots
  Lclen      (nonzeros per column in L)
  Lclen      (nonzeros per column in U)
  Lrowind[0]
  Lvalues[0]
  Urowind[0]
  Uvalues[0]
  Lrowind[1]
  Lvalues[1]
  Urowind[1]
  Uvalues[1]
  ...
*/
  
static
int Lappendcol(taucs_io_handle* LU, int j, int nnz, int* ind, taucs_datatype* re)
{
  taucs_io_append(LU,
		  0 + HEADER+j*4,
		  nnz,1,
		  TAUCS_INT,
		  ind);
  taucs_io_append(LU,
		  1 + HEADER+j*4,
		  nnz,1,
		  TAUCS_CORE_DATATYPE,
		  re);
  return 0;
}

static
int Uappendcol(taucs_io_handle* LU, int j, int nnz, int* ind, taucs_datatype* re)
{
  taucs_io_append(LU,
		  2 + HEADER+j*4,
		  nnz,1,
		  TAUCS_INT,
		  ind);
  taucs_io_append(LU,
		  3 + HEADER+j*4,
		  nnz,1,
		  TAUCS_CORE_DATATYPE,
		  re);
  return 0;
}

static
int Lreadcol(taucs_io_handle* LU, int j, int nnz, int* ind, taucs_datatype* re)
{
  taucs_io_read(LU,
		0 + HEADER+j*4,
		nnz,1,
		TAUCS_INT,
		ind);
  taucs_io_read(LU,
		1 + HEADER+j*4,
		nnz,1,
		TAUCS_CORE_DATATYPE,
		re);
  return 0;
}

static
int Ureadcol(taucs_io_handle* LU, int j, int nnz, int* ind, taucs_datatype* re)
{
  taucs_io_read(LU,
		2 + HEADER+j*4,
		nnz,1,
		TAUCS_INT,
		ind);
  taucs_io_read(LU,
		3 + HEADER+j*4,
		nnz,1,
		TAUCS_CORE_DATATYPE,
		re);
  return 0;
}

/*********************************************************/
/* ANALYSIS PHASE                                        */
/*********************************************************/

/* BUFFER SIZES AND MEMORY MANAGEMENT */

static int get_iobufsize()
{
  return 65536; /* minimal size of I/O buffers for good performance */
}

/* COLUMN-HEAP OPERATIONS */

/* The last entry in each input run must be extracted last among
   all the entries from that run, otherwise the heap might
   overflow */

static int heap_compare(int* heap, int i, int j) 
{
  /* compare columns */

  if (heap[3*i+1] < heap[3*j+1]) return -1;
  if (heap[3*i+1] > heap[3*j+1]) return  1;

  /* same column, compare runs (and order within runs) */

  if (heap[3*i+2] < heap[3*j+2]) return -1;
  if (heap[3*i+2] > heap[3*j+2]) return  1;

  return 0;
}

static void heap_heapify(int* heap, int* heapsize,int p) 
{
  int r,l,smallest;
  int i,j,k;

  l = p * 2;
  r = l + 1;

  if (l != 2*p || r != 2*p+1)
    taucs_printf("heap: left-right computation failed\n");

  if (l < *heapsize && heap_compare(heap,l,p) < 0)
    smallest = l;
  else
    smallest = p;

  if (r < *heapsize && heap_compare(heap,r,smallest) < 0)
    smallest = r;
  
  if (smallest != p) {
    i = heap[3*p]; 
    j = heap[3*p+1]; 
    k = heap[3*p+2]; 

    heap[3*p]   = heap[3*smallest];
    heap[3*p+1] = heap[3*smallest+1];
    heap[3*p+2] = heap[3*smallest+2];

    heap[3*smallest]   = i;
    heap[3*smallest+1] = j;
    heap[3*smallest+2] = k;
    
    heap_heapify(heap,heapsize,smallest);
  }
}

#if 0
static void x_heap_insert(int* heap, int* heapsize, int i, int j, int k)
{
  int l = *heapsize;

  heap[3*l]   = i;
  heap[3*l+1] = j;
  heap[3*l+2] = k;

  (*heapsize)++;
}  

static void x_heap_extract_min(int* heap, int* heapsize, int* i, int* j, int* k) 
{
  int l,m,mv;

  if (*heapsize <= 0)
    taucs_printf("heap: Trying to extract from an expty heap\n");

  mv = INT_MAX;
  for (l=0; l<(*heapsize); l++) {
    if (heap[3*l+1] < mv) {
      m = l;
      mv = heap[3*l+1];
    }
  }
  
  *i = heap[3*m];
  *j = heap[3*m+1];
  *k = heap[3*m+2];

  (*heapsize)--;
  for (l=m; l<(*heapsize); l++) {
    heap[3*l]   = heap[3*(l+1)];
    heap[3*l+1] = heap[3*(l+1) + 1];
    heap[3*l+2] = heap[3*(l+1) + 2];
  }
}
#endif

static void heap_insert(int* heap, int* heapsize, int i, int j, int k)
{
  int child, parent;

  (*heapsize)++;

  child = (*heapsize-1);
  parent = child / 2;
  while (child > 0 && 
	 ((heap[3*parent+1] > j) || 
	  ((heap[3*parent+1] == j) && (heap[3*parent+2] > k)))) {
    heap[3*child]   = heap[3*parent];
    heap[3*child+1] = heap[3*parent+1];
    heap[3*child+2] = heap[3*parent+2];
    child = parent;
    parent = child / 2;
  }

  heap[3*child]   = i;
  heap[3*child+1] = j;
  heap[3*child+2] = k;
}

static void heap_extract_min(int* heap, int* heapsize, int* i, int* j, int* k) 
{
  if (*heapsize <= 0)
    taucs_printf("heap: Trying to extract from an expty heap\n");
  
  *i = heap[0];
  *j = heap[1];
  *k = heap[2];

  heap[0] = heap[3 * ((*heapsize)-1)];
  heap[1] = heap[3 * ((*heapsize)-1) + 1];
  heap[2] = heap[3 * ((*heapsize)-1) + 2];

  (*heapsize) --;

  heap_heapify(heap,heapsize,0);
}

/* SYMMETRIC SKELETON GRAPH OPERATIONS */

static char skel_basename[256];
static int* skel_buffer;
static int  skel_buffer_size;
static int  skel_buffer_ptr;
static int  skel_outfiles;
static int  skel_infiles;
static int  skel_outfile;

static char skel_inphase;
static char skel_outphase;

static int skel_get_lastcol;
static int skel_next;

static int skel_compare(const void* e1, const void* e2) 
{
  /* we sort according to columns */
  int j1, j2;
  const int* jp1;
  const int* jp2;

  jp1 = (const int*) e1;
  jp2 = (const int*) e2;

  j1 = jp1[1];
  j2 = jp2[1];
  if (j1 < j2) return -1;
  if (j1 > j2) return 1;

  /* same column, compare rows */

  j1 = jp1[0];
  j2 = jp2[0];
  if (j1 < j2) return -1;
  if (j1 > j2) return 1;
  
  return 0;
}

static void skel_init(char* basename)
{
  sprintf(skel_basename,"%s.ssort",basename);

  /* adjust remaining memory.
     We allocate 2 io buffers for the sort phase, plus an array
     of file pointers that is guarantreed to be smaller.
     We then free them and allocate another io buffer for the
     stack phase, but we do not free the skel buffer first,
     so 3 io buffers is a conservative estimate. */

  remaining_memory -= (double) (3*get_iobufsize());

  skel_buffer_size = (int)(remaining_memory) / (2*sizeof(int));

  skel_buffer      = (int*)taucs_malloc(skel_buffer_size * 2 * sizeof(int));
  skel_buffer_ptr  = 0;
  skel_outfiles    = 0;

  skel_outphase = 'e';
  skel_inphase  = 'o';

  skel_get_lastcol = -1;
  skel_next = 0;
  skel_outfile = -1;
}

static void skel_finalize()
{
  taucs_free(skel_buffer);

  remaining_memory += (double) (3*get_iobufsize());
}

static void skel_add(int i,int j) 
{
  if (skel_buffer_ptr < skel_buffer_size) {
    skel_buffer[2*skel_buffer_ptr]   = i;
    skel_buffer[2*skel_buffer_ptr+1] = j;
    skel_buffer_ptr++;
  } else {
    int     file;
    ssize_t io_size;
    char    fname[256];

    /* SORT THIS BUFFER */

    qsort(skel_buffer, skel_buffer_ptr, 2*sizeof(int), &skel_compare);


    /* WRITE OUT */
    sprintf(fname,"%s.%c.%d",skel_basename,skel_outphase,skel_outfiles);
    taucs_printf("oocsp_colanalyze: Writing out skel sort buffer <%s> (3)\n",fname);
    file = open(fname,O_WRONLY | O_CREAT,0644);
    if (file == -1)
      taucs_printf("oocsp_colanalyze: could not create skel sort file\n");
    io_size = write(file,skel_buffer,skel_buffer_ptr * 2 * sizeof(int));
    if (io_size != skel_buffer_ptr * 2 * sizeof(int))
      taucs_printf("oocsp_colanalyze: write to skel sort file failed\n");
    close(file);

    skel_outfiles++;
    skel_buffer_ptr = 0;

    taucs_printf("oocsp_colanalyze: done (3)\n",fname);
  }
}

static void skel_sort_incore(int* postorder, int ncols, int* inv_postorder)
{
  int natural,i,col;

  natural = 1;
  for (i=0; i<ncols; i++) {
    if (postorder[i] != i) {
      natural = 0;
      break;
    }
  }

  if (natural == 0) {
    for (i=0; i<ncols; i++)
      inv_postorder[postorder[i]] = i;

    for (i=0; i<skel_buffer_ptr; i++) {
      col = skel_buffer[2*i+1];
      skel_buffer[ 2*i+1 ] = inv_postorder[col];
    }
  }

  qsort(skel_buffer, skel_buffer_ptr, 2*sizeof(int),
	&skel_compare);
}

static void skel_sort_outofcore(int* postorder, int ncols, int* inv_postorder)
{
  int  natural,i,j,k,e,f;
  int  heapsize, iobufsize, openruns, maxopenruns, runstart;
  int* infiles;
  int  outfile;

  int *inbuf;
  int *outbuf;
  int outbuf_ptr;

  int     file;
  ssize_t io_size;
  ssize_t io_count;
  char    fname[256];

  int last_col, last_extracted;

  /* FIRST, WRITE OUT THIS BUFFER */
  if (skel_buffer_ptr > 0) {
    /* SORT THIS BUFFER */

    qsort(skel_buffer, skel_buffer_ptr, 2*sizeof(int), &skel_compare);

    /* WRITE OUT */
    sprintf(fname,"%s.%c.%d",skel_basename,skel_outphase,skel_outfiles);
    taucs_printf("oocsp_colanalyze: Writing out skel sort buffer <%s> (1)\n",fname);
    file = open(fname,O_WRONLY | O_CREAT,0644);
    if (file == -1)
      taucs_printf("oocsp_colanalyze: could not create skel sort file\n");
    io_size = write(file,skel_buffer,skel_buffer_ptr * 2 * sizeof(int));
    if (io_size != skel_buffer_ptr * 2 * sizeof(int))
      taucs_printf("oocsp_colanalyze: write to skel sort file failed\n");
    close(file);

    skel_outfiles++;
    skel_buffer_ptr = 0;
  }

  /* DO WE NEED TO SORT THE RUNS AGAIN IN POSTORDER? */
  natural = 1;
  for (i=0; i<ncols; i++) {
    if (postorder[i] != i) {
      natural = 0;
      break;
    }
  }

  if (natural == 0) {
    for (i=0; i<ncols; i++)
      inv_postorder[postorder[i]] = i;

    for (f=0; f<skel_outfiles; f++) {
      sprintf(fname,"%s.%c.%d",skel_basename,skel_outphase,f);
      taucs_printf("oocsp_colanalyze: Resorting skel sort file <%s>\n",fname);
      file = open(fname,O_RDWR);
      if (file == -1)
	taucs_printf("oocsp_colanalyze: could not open skel sort file\n");
      /* read the run */
      io_size = read(file,skel_buffer,skel_buffer_size*2*sizeof(int));
      if (io_size == -1)
	taucs_printf("oocsp_colanalyze: read from skel sort file failed\n");
      io_count = io_size/(2*sizeof(int));
      close(file);

      /* sort again */

      for (i=0; i<(int)io_count; i++) {
	j = skel_buffer[2*i+1];
	skel_buffer[ 2*i+1 ] = inv_postorder[j];
      }

      qsort(skel_buffer, io_count, 2*sizeof(int),
	    &skel_compare);

      /* rewind the file and write back */
      file = open(fname,O_RDWR);
      if (file == -1)
	taucs_printf("oocsp_colanalyze: could not open skel sort file\n");
      /*      lseek(file,0,SEEK_SET);*/
      io_size = write(file,skel_buffer,io_count * 2 * sizeof(int));
      if (io_size != io_count * 2 * sizeof(int))
	taucs_printf("oocsp_colanalyze: write to skel sort file failed\n");
      close(file);
    }
  }

  /* WE NOW USE THE BUFFER AS A HEAP */
  
  heapsize = 0;

  /* in blocks of 2*sizeof(int), not bytes! */
  iobufsize   = get_iobufsize() / (2*sizeof(int)); 

  /* 
     each element in the skel_buffer is 2 ints, but 
     each element in the heap is 3 ints 
  */
  maxopenruns = (2*skel_buffer_size) / (3*iobufsize);
  if (maxopenruns < 2) {
    maxopenruns = 2;
    iobufsize   = (2*skel_buffer_size) / (3 * 2);
  }
  taucs_printf("oocsp_colanalyze: Using io buffers of %d elements (%d bytes), max runs = %d\n",
	     iobufsize,iobufsize*2,maxopenruns);

  inbuf  = (int*)taucs_malloc(iobufsize*2*sizeof(int));
  outbuf = (int*)taucs_malloc(iobufsize*2*sizeof(int));
  infiles  = (int*)taucs_malloc(skel_outfiles*sizeof(int));
  /*
  inbuf  = taucs_calloc(iobufsize,2*sizeof(int));
  outbuf = taucs_calloc(iobufsize,2*sizeof(int));
  infiles  = taucs_calloc(skel_outfiles,sizeof(int));
  */
  while (skel_outfiles > 1) {
    char phase;
    /*    int  i,j,k,runstart,openruns;*/

    taucs_printf("oocsp_colanalyze: Starting another merge phase with %d input runs\n",
	       skel_outfiles);

    skel_infiles  = skel_outfiles;
    skel_outfiles = 0;

    phase         = skel_inphase;
    skel_inphase  = skel_outphase;
    skel_outphase = phase;

    for (runstart=0; runstart<skel_infiles; runstart += maxopenruns) {

      sprintf(fname,"%s.%c.%d",
	      skel_basename,skel_outphase,runstart/maxopenruns);
      taucs_printf("oocsp_colanalyze: Opening output run <%s>\n",fname);
      outfile = open(fname,O_WRONLY | O_CREAT,0644);
      if (outfile == -1)
	taucs_printf("oocsp_colanalyze: could not open skel sort output file\n");
      skel_outfiles++;

      for (openruns=0; 
	   openruns < maxopenruns && runstart+openruns < skel_infiles;
	   openruns++) {
	
	sprintf(fname,"%s.%c.%d",skel_basename,skel_inphase,runstart+openruns);
	infiles[openruns] = open(fname,O_RDONLY);
	taucs_printf("oocsp_colanalyze: Opening input run <%s> (%d)\n",fname,infiles[openruns]);
	if (infiles[openruns] == -1)
	  taucs_printf("oocsp_colanalyze: could not open skel sort input file\n");
	io_size = read(infiles[openruns],inbuf,iobufsize*2*sizeof(int));
	if (io_size == -1)
	  taucs_printf("oocsp_colanalyze: read from skel sort file failed\n");
	io_count = io_size/(2*sizeof(int));
	if (io_count == 0) {
	  close(infiles[openruns]);
	  unlink(fname);
	}

	taucs_printf("oocsp_colanalyze: files %d %d\n",infiles[0],infiles[1]);

	taucs_printf("oocsp_colanalyze: Inserting %d elements from input run into heap\n",io_count);
	last_col = -1;
	for (e=0; e<(int)io_count; e++) {
	  i = inbuf[2*e];
	  j = inbuf[2*e+1];
	  /*printf("ij = %d\t%d\n",i,j);*/
	  if (last_col > j) {
	    taucs_printf("oocsp_colanalyze: > last = %d col = %d, (%d %d)\n",last_col,j,
		      e,io_count);
	    taucs_printf("oocsp_colanalyze: input run not sorted!\n");
	  }
	  last_col = j;
	  if (e == (int)io_count-1) /* end of inbuf marker */
	    heap_insert(skel_buffer,&heapsize,i,j,2*openruns+1);
	  else
	    heap_insert(skel_buffer,&heapsize,i,j,2*openruns);
	}
      } 

      taucs_printf("oocsp_colanalyze: files %d %d\n",infiles[0],infiles[1]);

      taucs_printf("oocsp_colanalyze: heapsize = %d\n",heapsize);

      outbuf_ptr = 0;

      last_extracted = -1;
      while (heapsize > 0) {
	int end_of_run,run;

	heap_extract_min(skel_buffer,&heapsize,&i,&j,&k);
	if (last_extracted > j)
	  taucs_printf("oocsp_colanalyze: heap order error!\n");
	last_extracted = j;
	outbuf[2*outbuf_ptr]   = i;
	outbuf[2*outbuf_ptr+1] = j;
	if (k % 2 == 1) end_of_run = 1;
	else end_of_run = 0;
	run = k / 2;
	outbuf_ptr ++;
	/*
	printf("extracted (%d,%d,%d) heapsize = %d\n",i,j,k,heapsize);
	*/


	if (end_of_run) {
	  taucs_printf("oocsp_colanalyze: Reading more from run %d\n",run);
	  io_size = read(infiles[run],inbuf,iobufsize * 2 * sizeof(int));
	  if (io_size == -1) {
	    taucs_printf("oocsp_colanalyze: errno = %d (%d)\n",errno,infiles[run]);
	    taucs_printf("oocsp_colanalyze: read from skel sort file failed\n");
	  }
	  io_count = io_size/(2*sizeof(int));
	  if (io_count == 0) {
	    sprintf(fname,"%s.%c.%d",skel_basename,skel_inphase,runstart+run);
	    taucs_printf("oocsp_colanalyze: Closing input run %d <%s>\n",run,fname);
	    close(infiles[run]);
	    unlink(fname);
	  }

	  taucs_printf("oocsp_colanalyze: Extracted %d,%d,%d\n",i,j,k);
	  taucs_printf("oocsp_colanalyze: heapsize = %d\n",heapsize);
	  /*
	  { 
	    int ii; 
	    printf("heap: \n");
	    for (ii=0; ii<heapsize; ii++) printf("(%d,%d,%d) ",
						    skel_buffer[3*ii+0],
						    skel_buffer[3*ii+1],
						    skel_buffer[3*ii+2]);
	    printf("\n");
	  }
	  */
	  taucs_printf("oocsp_colanalyze: Inserting %d elements from input run into heap\n",io_count);

	  last_col = -1;
	  for (e=0; e<(int)io_count; e++) {
	    i = inbuf[2*e];
	    j = inbuf[2*e+1];

	    if (last_col > j) 
	      taucs_printf("oocsp_colanalyze: input run not sorted!\n");
	    last_col = j;

	    if (e == (int)io_count-1) /* end of inbuf marker */
	      heap_insert(skel_buffer,&heapsize,i,j,2*run+1);
	    else
	      heap_insert(skel_buffer,&heapsize,i,j,2*run);

	    if (3*heapsize >= 2*skel_buffer_size) {
	      taucs_printf("oocsp_colanalyze: heapsize = %d, buffer_size = %d\n",
			heapsize,2*skel_buffer_size);
	      taucs_printf("oocsp_colanalyze: merge-heap overflow\n");
	    }
	  }

	  taucs_printf("oocsp_colanalyze: heapsize = %d\n",heapsize);
	}

	if (outbuf_ptr >= iobufsize) {
	  taucs_printf("oocsp_colanalyze: Writing to output run\n");
	  taucs_printf("oocsp_colanalyze: heapsize = %d\n",heapsize);
	  io_count = iobufsize;
	  io_size = write(outfile,outbuf,io_count * 2 * sizeof(int));
	  if (io_size != io_count * 2 * sizeof(int))
	    taucs_printf("oocsp_colanalyze: write to skel sort file failed\n");
	  outbuf_ptr = 0;
	}
      }

      /* write the rest of the output */
      if (outbuf_ptr > 0) {
	taucs_printf("oocsp_colanalyze: Writing to output run and closing\n");
	io_count = outbuf_ptr;
	io_size = write(outfile,outbuf,io_count * 2 * sizeof(int));
	if (io_size != io_count * 2 * sizeof(int))
	  taucs_printf("oocsp_colanalyze: write to skel sort file failed\n");
	outbuf_ptr = 0;
      }
      close(outfile);
      skel_buffer_ptr = 0;
    }
  }


  taucs_free(infiles);
  taucs_free(inbuf);
  taucs_free(outbuf);

  /*
  taucs_free(infiles);
  taucs_free(inbuf);
  taucs_free(outbuf);
  */
}

static int  stack_allocated;
static int  stack_buffer_size;
static int* stack_buffer;

static void skel_sort(int* postorder, int ncols,int* tmp)
{
  int  file;
  char fname[256];
  ssize_t io_size;
  int iobufsize;

  iobufsize   = get_iobufsize() / (2*sizeof(int)); 

  if (skel_outfiles == 0) {
    skel_sort_incore(postorder,ncols,tmp);

    if (skel_buffer_ptr <= skel_buffer_size / 2) {
      stack_buffer      = skel_buffer + (2*skel_buffer_ptr);
      stack_buffer_size = 2*(skel_buffer_size - skel_buffer_ptr);
      stack_allocated   = 0;
      taucs_printf("oocsp_colanalyze: Using remainder of skeleton buffer for stack,\n");
      taucs_printf("oocsp_colanalyze: size = %d ints\n",stack_buffer_size);
    } else {
      sprintf(fname,"%s.%c.%d",skel_basename,skel_outphase,0);
      taucs_printf("oocsp_colanalyze: Writing out skel sort buffer <%s> (2)\n",fname);
      file = open(fname,O_WRONLY | O_CREAT,0644);
      if (file == -1)
	taucs_printf("oocsp_colanalyze: could not create skel sort file\n");
      io_size = write(file,skel_buffer,skel_buffer_ptr * 2 * sizeof(int));
      if (io_size != skel_buffer_ptr * 2 * sizeof(int))
	taucs_printf("oocsp_colanalyze: write to skel sort file failed\n");
      close(file);

      skel_outfiles++;
      skel_buffer_ptr = 0;

      stack_buffer      = skel_buffer;
      stack_buffer_size = 2*skel_buffer_size;
      stack_allocated   = 1; /* we need to free it */

      skel_buffer       = (int*)taucs_malloc(iobufsize*2*sizeof(int));
      skel_buffer_size  = iobufsize;
      taucs_printf("oocsp_colanalyze: Using skeleton buffer for stack, allocating \n");
      taucs_printf("oocsp_colanalyze: new skeleton buffer\n\n");
    }
  }
  else {
    skel_sort_outofcore(postorder,ncols,tmp);

    stack_buffer      = skel_buffer;
    stack_buffer_size = 2*skel_buffer_size;
    skel_buffer       = (int*)taucs_malloc(iobufsize*2*sizeof(int));
    skel_buffer_size  = iobufsize;
    taucs_printf("oocsp_colanalyze: Using skeleton buffer for stack, allocating \n");
    taucs_printf("oocsp_colanalyze: new skeleton buffer\n");
  }
}

static int skel_get_next(int j)
{
  int row, col;
  char fname[256];
  ssize_t io_size;

  if (skel_next >= skel_buffer_ptr) {
    if (skel_outfiles > 0) {
      if (skel_outfile == -1) {
	sprintf(fname,"%s.%c.%d",skel_basename,skel_outphase,0);
	taucs_printf("oocsp_colanalyze: Opening skel sort buffer <%s> (2)\n",fname);
	skel_outfile = open(fname,O_RDONLY);
	if (skel_outfile == -1)
	  taucs_printf("oocsp_colanalyze: could not open skel sort file\n");
      }
      io_size = read(skel_outfile,
		     skel_buffer,skel_buffer_size * 2 * sizeof(int));
      if (io_size == -1) 
	taucs_printf("oocsp_colanalyze: I/O error while trying to read skel sort file\n");
      if (io_size == 0) { /* end of file */
	taucs_printf("oocsp_colanalyze: Closing and removing skel file, col=%d\n",j);
	sprintf(fname,"%s.%c.%d",skel_basename,skel_outphase,0);
	close(skel_outfile);
	unlink(fname);
	skel_outfiles   = 0;
	skel_next       = 0;
	skel_buffer_ptr = 0;
	return -1;
      } else {
	skel_next       = 0;
	skel_buffer_ptr = io_size / (2*sizeof(int));
	taucs_printf("oocsp_colanalyze: read %d elements from skel sort file\n",skel_buffer_ptr);
      }
    } else
      return -1;
  }

  if (skel_next >= skel_buffer_ptr)
    return -1;

  row = skel_buffer[2*skel_next];
  col = skel_buffer[2*skel_next+1];
  if (col == j) {
    skel_next++;
    return row;
  }
  else
    return -1;
}
    

static void skel_get_postordercol(int* found, int flag,
				  int j,
				  int* nnz, int* rowind)
{
  int row;

  *nnz = 0;
  while ((row = skel_get_next(j)) != -1) {
    if (found[row] < flag) {
      found[row] = flag;
      rowind[ *nnz ] = row;
      (*nnz)++;
    }
  }
}


/*
static void skel_get_postordercol(int* found, int flag,
				  int j,
				  int* nnz, int* rowind)
{
  int row;
  int i;
  int next;

  next = skel_get_next;
  while (skel_buffer[2*next + 1] < j && next < skel_buffer_ptr) {
    next++;
  }

  if (skel_buffer[2*next + 1] == j && next < skel_buffer_ptr) {
    if (next != skel_get_lastcol+1)
      taucs_printf("oocsp_colanalyze: internal error in get_postordercol\n");
  }

  *nnz = 0;
  for (i = next;
       skel_buffer[2*i + 1] == j && i < skel_buffer_ptr;
       i++) {
    row = skel_buffer[2*i];
    if (found[row] < flag) {
      found[row] = flag;
      rowind[ *nnz ] = row;
      (*nnz)++;
    }
  }

  skel_get_next = i;
  skel_get_lastcol = i - 1;
}
*/

/* UNION FIND ROUTINES */

static int uf_makeset(int* uf, int i)        { uf[i] = i; return i; }
static int uf_union  (int* uf, int s, int t) { uf[s] = t; return t; }
static int uf_find   (int* uf, int i)        { if (uf[i] != i) 
                                                 uf[i] = uf_find(uf,uf[i]); 
                                               return uf[i]; }

/* FILL STACK ROUTINES */

static int  stack_files;
static char stack_basename[256];

/*
static int  stack_buffer_size;
static int* stack_buffer;
*/
static int  stack_buffer_ptr;

static int  stack_top;

static double stack_size;
static double stack_max_size;

static void stack_init(char* basename,
		       int* colptr, int* colstack, 
		       int ncols)
{
  int j;

  sprintf(stack_basename,"%s.fstack",basename);
  stack_files      = 0;
  stack_buffer_ptr = 0;
  stack_top        = -1;

  for (j=0; j<ncols; j++) colptr[j] = -1;

  /*
  stack_buffer_size = 1048576;
  stack_buffer = mxCalloc(stack_buffer_size,sizeof(int));
  */

  stack_size     = 0.0;
  stack_max_size = 0.0;
}

static void stack_finalize()
{
  if (stack_files != 0 || stack_buffer_ptr != 0)
    taucs_printf("oocsp_colanalyze: fill stack did not get empty\n");

  if (stack_allocated)
    taucs_free(stack_buffer);

  taucs_printf("oocsp_colanalyze: max stack size = %.0lf\n",stack_max_size);
}

static void stack_push(int* colptr, int* colstack, int i, int j)
{
  if (stack_top < 0 || colstack[stack_top] != j) {
    if (colptr[j] != -1) {
      taucs_printf("oocsp_colanalyze: fill stack internal error (push)\n");
    }
    stack_top++;
    colstack[stack_top] = j;
    colptr[j] = (stack_buffer_size*stack_files) + stack_buffer_ptr;
    if (colptr[j] > (INT_MAX/2))
      taucs_printf("oocsp_colanalyze: Warning! Pointers to fill stack may overflow\n");
  }
   
  stack_buffer[stack_buffer_ptr] = i;
  stack_buffer_ptr++;

  if (stack_buffer_ptr >= stack_buffer_size) {
    int     file;
    ssize_t io_size;
    char    fname[256];

    sprintf(fname,"%s.%d",stack_basename,stack_files);
    taucs_printf("oocsp_colanalyze: Writing out fill stack buffer <%s>\n",fname);
    file = open(fname,O_WRONLY | O_CREAT,0644);
    if (file == -1)
      taucs_printf("oocsp_colanalyze: could not create stack file\n");
    io_size = write(file,stack_buffer,stack_buffer_size * sizeof(int));
    if (io_size != stack_buffer_size * sizeof(int))
      taucs_printf("oocsp_colanalyze: write to stack file failed\n");
    close(file);

    stack_files++;
    stack_buffer_ptr = 0;
  }

  /*
  stack_size++;
  if (stack_size > stack_max_size) stack_max_size = stack_size;
  */
}

static void stack_pop(int* colptr, int* colstack, 
		      int* found, int flag,
		      int j, int* nnz, int* rowind)
{
  int row;
  int   i;

  if (stack_top < 0 || colstack[stack_top] != j) { /* empty fill column */
    for (i=0; i<=stack_top; i++)
      if (colstack[i] == j)
	taucs_printf("oocsp_colanalyze: fill stack internal error (pop)\n");
    
    *nnz    = 0;
  } else {
    *nnz = 0;
    while (colptr[j] < (stack_buffer_size*stack_files) +stack_buffer_ptr) {

      if (stack_buffer_ptr == 0) {
	int     file;
	ssize_t io_size;
	char    fname[256];
	
	stack_files--;
	stack_buffer_ptr = stack_buffer_size;
	sprintf(fname,"%s.%d",stack_basename,stack_files);
	taucs_printf("oocsp_colanalyze: Reading a fill stack buffer <%s>\n",fname);
	file = open(fname,O_RDONLY);
	if (file == -1)
	  taucs_printf("oocsp_colanalyze: could not open stack file\n");
	io_size = read(file,stack_buffer,stack_buffer_size * sizeof(int));
	if (io_size != stack_buffer_size * sizeof(int))
	  taucs_printf("oocsp_colanalyze: read from stack file failed\n");
	close(file);
	unlink(fname);
      }

      stack_buffer_ptr --;
      row = stack_buffer[ stack_buffer_ptr ];

      if (found[row] < flag) {
	found[row] = flag;
	rowind[ *nnz ] = row;
	(*nnz)++;
      }
    }

    stack_top--;
    colptr[j] = -1;
  }
  /*
  stack_size -= (double) (*nnz);
  */
}

/* MAIN ROUTINE */

static
void oocsp_colanalyze(taucs_ccs_matrix* matrix,
		      char* basename,
		      int*  colperm,
		      int** ptrparent,
		      int** ptrpostorder,
		      int** ptrlcolcount,
		      int** ptrucolcount)
{
  int i,j,ip,p,jp;
  int nnz,cset,rset,rroot,fcol;
  
  int  postnum, depth;
  int* first_kid;
  int* next_kid;
  int* stack_vertex;
  int* stack_child;

  /*int* colptr;*/
  int* rowind;
  int* firstcol;
  int* uf;
  int* found;
  int* root;
  int* stack_colptr;
  int* stack_colstk;
  int* tmp_col;
  
  int* parent;
  int* postorder;
  int* lcolcount;
  int* ucolcount;  
  
  int*    nrows;
  int*    ncols;

  taucs_printf("oocsp_colanalyze: In colanalyze\n");
  taucs_printf("oocsp_colanalyze: using %.0lf MBytes of memory\n",(remaining_memory)/1048576.0);
  
  nrows     = &matrix->m;
  ncols     = &matrix->n;

  /* START THE ANALYSYS */

  skel_init(basename);

  (remaining_memory) -= (double) ( 4 * (*ncols) * 4 /* sizeof(int32) */);

  /* +1 for stack_vertex, stack_child */
  parent          = (int*)taucs_malloc((*ncols+1)*sizeof(int));
  *ptrparent = parent;
  lcolcount       = (int*)taucs_malloc((*ncols+1)*sizeof(int));
  *ptrlcolcount = lcolcount;
  ucolcount       = (int*)taucs_malloc((*ncols+1)*sizeof(int));
  *ptrucolcount = ucolcount;
  postorder       = (int*)taucs_malloc((*ncols+1)*sizeof(int));
  *ptrpostorder = postorder;

  (remaining_memory) -= (double) ( 2 * ((*ncols)+1) * sizeof(int));
  (remaining_memory) -= (double) ( 2 * (*nrows) * sizeof(int));

  uf        = (int*)taucs_malloc((*ncols+1)*sizeof(int));
  root      = (int*)taucs_malloc((*ncols+1)*sizeof(int));
  firstcol  = (int*)taucs_malloc((*nrows)*sizeof(int));

  tmp_col   = (int*)taucs_malloc((*nrows)*sizeof(int));

  /* we can reuse the same space */
  first_kid = uf;
  next_kid  = root;

  found        = firstcol;
  stack_colptr = uf;
  stack_colstk = root;

  /* we use the output arrays before they are used. */
  stack_vertex = lcolcount; 
  stack_child  = ucolcount;

  for (i=0; i < (*nrows); i++) {
    firstcol[i] = (*ncols);
  }

  for (j=0; j < (*ncols); j++) {

#ifdef OLD
    {
      int   file;
      off_t offset;
      ssize_t io_size;
	
      nnz = clen[j];
      /* taucs_printf("oocsp_colanalyze: Debug_1: j= %d nnz= %d\n",j,nnz);*/
      file   = dfile_fid[ ind_fid[j] ];
      offset = ind_off[j];
      if (nnz != 0 && file != -1 && offset != -1) {
	if (lseek(file,offset,SEEK_SET) == -1) taucs_printf("oocsp_colanalyze: lseek failed\n");
	io_size = read(file, tmp_col, nnz * sizeof(int));
	if (io_size != nnz * sizeof(int)) taucs_printf("oocsp_colanalyze: Error reading data.\n");
      }
      rowind = tmp_col;

    }
#else
    /*
    nnz = clen[j];
    oocsp_readcol_structure(matrix,j,tmp_col);
    rowind = tmp_col;
    */

    /* new code: Sivan 28 Feb 2002 */
    nnz = (matrix->colptr)[colperm[j]+1] - (matrix->colptr)[colperm[j]];
    rowind = (matrix->rowind) + (matrix->colptr)[colperm[j]];
#endif

    cset       = uf_makeset(uf, j);
    root[cset] = j;
    parent[j]  = (*ncols);
    
    for (ip=0; ip<nnz; ip++) {
      
      i = rowind[ip];
      fcol = firstcol[i];
      if (fcol >= j) {
	firstcol[i] = j;
	fcol        = j;
      } else {
	rset = uf_find(uf,fcol);
	rroot = root[rset]; 
	if (rroot != j) {
	  parent[rroot] = j;
	  cset          = uf_union(uf,cset,rset);
	  root[cset]    = j;
	}
      }
      /* ADD (j,fcol) TO SKELETON */
      skel_add(j,fcol);
    }
    /*
    mxDestroyArray(output_args[0]);
    */
  }

  /* COMPUTE POSTORDER OF ETREE */

  /* create linked lists of children */
  for (j=0; j <= *ncols; j++) first_kid[j] = -1;
  for (j = (*ncols)-1; j >= 0; j--) {
    p              = parent[j];
    next_kid[j]    = first_kid[p];
    first_kid[p] = j;
  }

  
  /* do dfs in a loop */
  postnum = 0;
  depth = 0;
  stack_vertex[depth] = *ncols; /* root */
  stack_child [depth] = first_kid[ stack_vertex[depth] ];
  while (depth >= 0) {
    if ( stack_child[depth] != -1 ) {
      stack_vertex[depth+1] = stack_child[depth];
      stack_child [depth+1] = first_kid[  stack_vertex[depth+1] ];
      depth++;
    } else {
      if ( stack_vertex[depth] != (*ncols) ) { /* not root */
	if (stack_vertex[depth] >= *ncols) 
	  taucs_printf("oocsp_colanalyze: internal error in dfs (0)\n");
	postorder[ postnum ] = stack_vertex[depth];
	postnum++;
      }
      depth--;
      if (depth >= 0) /* sivan June 30, to avoid crash, seems to be right */
	stack_child[depth] = next_kid[  stack_child[depth] ];
    }
    if (depth > *ncols) {
      int i;
      taucs_printf("oocsp_colanalyze: depth=%d, ncols=%d\n",depth,*ncols);
      for (i=0; i<(*ncols); i++) {
	taucs_printf("oocsp_colanalyze: %d: [%d %d]\n",i,first_kid[i],next_kid[i]);
      }
      taucs_printf("oocsp_colanalyze: internal error in dfs (1)\n");
    }
  }

  
  if (postnum != *ncols) taucs_printf("oocsp_colanalyze: internal error in dfs (2)\n");
    
  /* SORT THE SKELETON MATRIX */

  skel_sort(postorder, *ncols, 
	    found /* temporary */);


  
  /* SECOND PHASE, COMPUTE COLCOUNTS */

  /* we reuse the space of uf and root */
  stack_init(basename,stack_colptr,stack_colstk,*ncols); 

  for (i=0; i < (*nrows); i++) {
    found[i] = -1;
  }

  for (j=0; j<(*ncols); j++) {
    lcolcount[j] = 1;
    ucolcount[j] = 1;
  }

  for (jp=0; jp<(*ncols); jp++) {
    j = postorder[jp];
    found[j] = jp;
    p = parent[j];

    if (p < (*ncols)) {
      found[p] = jp;
      lcolcount[j]++;
      ucolcount[p]++;
    }

    stack_pop(stack_colptr,stack_colstk,
	      found, jp,
	      j,
	      &nnz,tmp_col);
    rowind = tmp_col;
    for (ip=0; ip<nnz; ip++) {
      i = rowind[ip];

      lcolcount[j]++;
      ucolcount[i]++;
      if (p < *ncols) stack_push(stack_colptr,stack_colstk,i,p);
    }

    skel_get_postordercol(found,jp,
			  jp, /* use postorder column index */
			  &nnz,tmp_col);
    rowind = tmp_col;
    for (ip=0; ip<nnz; ip++) {
      i = rowind[ip];

      lcolcount[j]++;
      ucolcount[i]++;

      if (lcolcount[j] > *(ncols))
	taucs_printf("oocsp_colanalyze: Internal error while producing lcolcounts\n");
      if (ucolcount[i] > *(ncols))
	taucs_printf("oocsp_colanalyze: Internal error while producing ucolcounts\n");

      if (p < *ncols) stack_push(stack_colptr,stack_colstk,i,p);
    }
  }

  stack_finalize();

  /*
  mxDestroyArray(tmp1_array);
  mxDestroyArray(tmp2_array);
  mxDestroyArray(tmp3_array);
  mxDestroyArray(tmp4_array);
  */

  taucs_free(uf);
  taucs_free(firstcol);
  taucs_free(root);
  taucs_free(tmp_col);

  (remaining_memory) += (double) ( 2 * ((*ncols)+1) * sizeof(int));
  (remaining_memory) += (double) ( 2 * (*nrows) * sizeof(int));

  skel_finalize();

  /* MAKE AND POSTORDER PARENT 1-BASED AND MARK ROOTS WITH A ZERO */

  /*for (j=0; j < (*ncols); j++) {
    postorder[j]++;
    if (parent[j] == (*ncols))
      parent[j] = 0;
    else 
      parent[j]++;
      }  */
  /* fix up by Vladi */

  for (j=0; j < (*ncols); j++) {
    if (parent[j] == (*ncols))
      parent[j] = -1;
    /*lcolcount[j]--;
      ucolcount[i]--;*/
  }

  taucs_printf("oocsp_colanalyze: done\n");
}

/*********************************************************/
/* PANELIZATION                                          */
/*********************************************************/


/* There seems to be a confusion here between spawidth and remaining memory; sivan */

static
void oocsp_panelize_simple(
			   int  nrows,             /* input  */  
			   int  ncols,             /* input  */
			   int* postorder,         /* input  */
			   int* l_colcounts,       /* input  */
			   int* u_colcounts,       /* input  */
			   int* parents,           /* input  */

                           int* spawidth,          /* input  */  
                           int* maxsn,             /* input  */  
			   int** ptrpanels,        /* output */
			   int** ptrschedstart,    /* output */
			   int** ptrschedend,      /* output */
			   int** ptrfetchnext,     /* output */
			   int** ptrejectnext      /* output */
			   )

{
  int i,j,c;
  int panelsize, panelcols, panelnumber;
  int stop,eject,newpanel;
  int colcount;
  int* panels;
  int* schedstart;
  int* schedend;
  int* fetchnext;
  int* ejectnext;  
  double maxpanelsize;
  double memuse, width_multiplier;
  int  maxcolcount;

  maxcolcount = 0;
  for (j=0; j<ncols; j++)
    maxcolcount = max( maxcolcount, l_colcounts[j]+u_colcounts[j] );

  taucs_printf("oocsp_panelize: max col count = %d, nrows = %d\n",
	     maxcolcount, nrows);

  /*
    memory usage in numfact, width > 0, snodes:
           8*nrows*I + 2*nrows*D + 2*nrows*C + 3*nrows*P  +
           2*ncols*I + 
           1*w*nrows*I + 1*w*nrows*D + 1*w*I +
           5*w*maxc *I + 2*w*maxc *D +
           1*s*I +
           1*s*maxc *I + 2*s*maxc *D 

    we also account for the memory required for the in-core representation of A, L, U,
    which is 2*ncols*sizeof(char)+3*ncols*sizeof(int)
  */

  *maxsn = 8 - 4; /* 8 is the minimum value we use here */
  do {
    *maxsn += 4;
    memuse = 
      8.0*nrows*sizeof(int) + 2.0*nrows*sizeof(taucs_datatype) 
      + 2.0*nrows*sizeof(char) + 3.0*nrows*sizeof(void*)
      + 2.0*ncols*sizeof(int)
      + 1.0*(*maxsn)*sizeof(int)
      + 1.0*(*maxsn)*maxcolcount*sizeof(int) + 2.0*(*maxsn)*maxcolcount*sizeof(taucs_datatype);

    memuse +=
      3 * ( 3.0*ncols*sizeof(int) + 2.0*ncols*sizeof(char) );

    width_multiplier = 
      1.0*nrows*sizeof(int) + 1.0*nrows*sizeof(taucs_datatype) + 1.0*sizeof(int)
      + 5.0*maxcolcount*sizeof(int) + 2.0*maxcolcount*sizeof(taucs_datatype);
    *spawidth = (int) floor( (remaining_memory - memuse) / width_multiplier );
  } while (*spawidth > 4*(*maxsn));

  if (*spawidth < 8) *spawidth = 8; /* it might go over the limit */

  maxpanelsize = nrows * *spawidth; 

  taucs_printf("oocsp_panelize: spawidth = %d, max supernode = %d\n",*spawidth,*maxsn);

  panels         = (int*)taucs_malloc(ncols*sizeof(int)); (*ptrpanels)     = panels;
  schedstart     = (int*)taucs_malloc(ncols*sizeof(int)); (*ptrschedstart) = schedstart;
  schedend       = (int*)taucs_malloc(ncols*sizeof(int)); (*ptrschedend)   = schedend;
  fetchnext      = (int*)taucs_malloc(ncols*sizeof(int)); (*ptrfetchnext)  = fetchnext;
  ejectnext      = (int*)taucs_malloc(ncols*sizeof(int)); (*ptrejectnext)  = ejectnext;

  panelnumber = 1;
  j           = 0;
  eject = 0;
  
  while (j<ncols)
  {
    panelsize = 0;
    panelcols = 0;
    stop = 0;

    while (stop == 0) 
      {
	c = postorder[j];
	colcount = l_colcounts[c] + u_colcounts[c]; /* was only L, don't remember why, sivan */
	newpanel = 0;
	if (*spawidth > 0)
	  {
	  if (panelcols >= *spawidth) 
	    newpanel=1; 
	  else
	    if (panelsize + colcount > maxpanelsize) newpanel=1; 
	  }

	if (newpanel == 0)
	  {
	    /* add c to panel */
	    panelsize = panelsize + colcount;
	    panelcols = panelcols + 1;
	    panels[c] = panelnumber;
	    ejectnext[j] = c;
	    schedstart[c] = panelnumber;
	    schedend[c]   = panelnumber;
	    j = j+1;
	    if (j >= ncols) stop = 1; 
	  } 
	else
	  {
	    panelnumber = panelnumber + 1;
	    stop = 1;
	  }
      }
    
    /* now reverse panel, compute fetchnext */
#ifdef VLADIMIR
    if(j<ncols)
      for (i=0;i<j-eject;i++)
	fetchnext[eject+i] = ejectnext[j-i];
    else {
      printf("j >= ncols???\n");
      exit(1);
      for (i=0;i<j-eject;i++)
	fetchnext[eject+i] = ejectnext[j-i-1];
    }
#else
    for (i=0;i<j-eject;i++)
      fetchnext[eject+i] = ejectnext[j-i-1];
#endif

    eject = j;
     
  }
}


/*********************************************************/
/* NUMERICAL PHASE                                       */
/*********************************************************/

#define SNODES
#define SNODE_THRESHOLD 4
#define SNODE_BLOCK 8
#define SIMPLE_COL_COL_no
#define SPA_ONEARRAY
#define USE_BLAS
#define BLAS_THRESHOLD 10
#define BLOCK 16

#ifdef DETAILED_TIMING
static double flops_extra = 0.0;
static double flops_dense = 0.0;
#endif

/*
  Out-of-core sparse LU

  Numerical factorization.

  Memory management
    ipivots              nrows*I                            ! not freed, returned to caller
    rowlists             nrows*I + 3*width*maxlucols*I      ! heads,colind,next,prev
    heap                 nrows*I                            !                        

    spa                  nrows*I + nrows*D                  ! if spawidth < 0
    spa (width>0)        width*nrows*I + width*nrows*D      !

    lindices             nrows*C                            ! L or U 
    nnzmap               nrows*C                            ! ?

    panel                2*ncols*I                          ! id, nnz
                         + 3*ncols*P                        ! ind, inrowlist, re

    panel_compressed     sum_j lucols_j*(D+2*I), j=1:width  ! re, ind, inrowlist
                         <= width*maxlucols*(D+2*I)
                         
    snodes               2*ncols*I                          ! pivots, snode_index
                         +   MAX_SNODE*I                    ! pivrows
			 + 2*nrows*I                        ! ind, map
			 +   nrows*D                        ! re ???
			 +   MAX_SNODE*maxcolcount*I        ! lu_ind
			 +   MAX_SNODE*maxcolcount*D        ! lu_re

    snodes_dense         SNODE_MAX*maxcolcont*D             ! S
                         + width*maxcolcount*D              ! P
                         + width*I                          ! updcols

    no_snodes            nrows*I + nrows*D                  ! lu_ind, lu_re

    total: width > 0, snodes
           8*nrows*I + 2*nrows*D + 2*nrows*C + 3*nrows*P  +
           2*ncols*I + 
           1*w*nrows*I + 1*w*nrows*D + 1*w*I +
           5*w*maxc *I + 2*w*maxc *D +
           1*s*I +
           1*s*maxc *I + 2*s*maxc *D 

    total: width < 0, snodes
           8*nrows*I + 2*nrows*D + 2*nrows*C + 3*nrows*P  +
           2*ncols*I + 
           1*w*I +
           5*w*maxc *I + 2*w*maxc *D +
           1*s*I +
           1*s*maxc *I + 2*s*maxc *D 

    total: width > 0, no snodes
           7*nrows*I + 2*nrows*D + 2*nrows*C + 3*nrows*P  +
           1*w*nrows*I + 1*w*nrows*D +
           5*w*maxc *I + 1*w*maxc *D +

    total: width < 0, no snodes
           7*nrows*I + 2*nrows*D + 2*nrows*C + 3*nrows*P  +
           5*w*maxc *I + 1*w*maxc *D +

    where:
      I=sizeof(int) D=sizeof(taucs_datatype) C=sizeof(char) P=sizeof(void*)
      maxc=maxcolcount w=spawidth (=active panelwidth) s=max_snode
*/



static double time_total;

#ifdef DETAILED_TIMING
static double time_colcol;
static double time_colcol_1;
static double time_colcol_2;
static double time_factor;
static double time_scatter;
static double time_gather;
static double time_append;
static double time_read;
#ifdef SNODES
static double time_snode_tmp;
static double time_snode_1 = 0;
static double time_snode_2 = 0;
static double time_snode_21 = 0;
static double time_snode_3 = 0;
static double time_snode_4 = 0;
static double time_snode_detect;
static double time_snode_prepare;
static double time_snode_dense;
#endif

static double bytes_read;
static double bytes_appended;
static double col_ooc_updates;
static double col_read;

static double  flops;
static double  scatters;
static double  gathers;
static double  rowlist_ops;
static double  num_heap_ops;
#endif /* DETAILED_TIMING */

/****************************************************/
/*                                                  */
/* Heap operations                                  */
/*                                                  */
/****************************************************/

/* HEAP OPERATIONS */


static void num_heap_heapify(int* heap, int* heapsize, 
			     int* ipivots, int p) 
{
  int r,l,smallest;
  int temp;

#ifdef DETAILED_TIMING
  num_heap_ops += 1.0;
#endif

  r = (p+1) * 2;
  l = r - 1;

#ifdef INTERNAL_CHECKS
  if ((l-1)/2 != p || (r-1)/2 != p)
    taucs_printf("oocsp_numfact: left-right computation in heap failed\n");
#endif

  if (l < *heapsize && ipivots[heap[l]] < ipivots[heap[p]])
    smallest = l;
  else
    smallest = p;

  if (r < *heapsize && ipivots[heap[r]] < ipivots[heap[smallest]])
    smallest = r;
  
  if (smallest != p) {
    temp           = heap[p]; 
    heap[p]        = heap[smallest];
    heap[smallest] = temp;
    
    num_heap_heapify(heap, heapsize, ipivots, smallest);
  }
}

static void num_heap_insert(int* heap, int* heapsize, int* ipivots, int i)
{
  int child, parent;

  (*heapsize)++;

#ifdef DETAILED_TIMING
  num_heap_ops += 1.0;
#endif

  child = (*heapsize-1);
  parent = (child-1) / 2;
  while (child > 0 && (ipivots[heap[parent]] > ipivots[i])) {
    heap[child]   = heap[parent];
    child = parent;
    parent = (child-1) / 2;

#ifdef DETAILED_TIMING
    num_heap_ops += 1.0;
#endif
  }

  heap[child]   = i;
}

static int num_heap_extractmin(int* heap, int* heapsize, int* ipivots) 
{
  int m; 

#ifdef DETAILED_TIMING
  num_heap_ops += 1.0;
#endif

  if (*heapsize <= 0) return -1;
  
  m = heap[0];
  
  heap[0] = heap[(*heapsize)-1];

  (*heapsize)--;

  num_heap_heapify(heap,heapsize,ipivots,0);

  return m;
}

/****************************************************/
/*                                                  */
/* Row lists                                        */
/*                                                  */
/****************************************************/

static int* rowlists_head;   /* one head per row */
static int* rowlists_colind; 
static int* rowlists_next;
static int* rowlists_prev;
static int  rowlists_size;
static int  rowlists_freehead;
/*static int  rowlists_freenext;*/

static void rowlists_finalize()
{
  taucs_free(rowlists_head);
  taucs_free(rowlists_colind);
  taucs_free(rowlists_next);
  taucs_free(rowlists_prev);
}

static void rowlists_init(int size, int nrows)
{
  int i;

  rowlists_size = size;

  rowlists_head   = (int*)taucs_malloc(nrows*sizeof(int));
  rowlists_colind = (int*)taucs_malloc(rowlists_size*sizeof(int));
  rowlists_next   = (int*)taucs_malloc(rowlists_size*sizeof(int));
  rowlists_prev   = (int*)taucs_malloc(rowlists_size*sizeof(int));
  assert(rowlists_head && rowlists_colind && rowlists_next && rowlists_prev); 

  for (i=0; i<nrows; i++) rowlists_head[i] = -1;

  /* link the entire rowlist as one freelist */

  rowlists_freehead = 0;
  for (i=0; i<rowlists_size; i++) {
    rowlists_next[i] = i+1;
    /* freelist does not need prev pointers */ 
    /* rowlists_prev[i] = i-1; */ 
  }
  rowlists_next[ rowlists_size - 1 ] = -1;
}

static int rowlists_insert(int row, int panelcol)
{
  int new;

#ifdef DETAILED_TIMING
  rowlist_ops += 1.0;
#endif /* DETAILED_TIMING */

  /* get memory from the freelist */

  if ((new = rowlists_freehead) == -1) {
    taucs_printf("oocsp_numfact: Out of rowlist memory\n");
    exit(1);
  }

  /* remove this memory from the freelist; freelist does now use prev */

  rowlists_freehead = rowlists_next[ new ];

  /* link to row list */

  rowlists_next[ new ] = rowlists_head[ row ];
  rowlists_prev[ new ] = -1;
  rowlists_colind[ new ] = panelcol;

  if (rowlists_next[new] != -1)
    rowlists_prev[ rowlists_next[new] ] = new;

  rowlists_head[ row ] = new;

  return new;
}

static void rowlists_delete(int row, int index)
{
#ifdef DETAILED_TIMING
  rowlist_ops += 1.0;
#endif /* DETAILED_TIMING */

  if (rowlists_head[ row ] == index)
    rowlists_head[ row ] = rowlists_next[ index ];

  if (rowlists_next[ index ] != -1)
    rowlists_prev[ rowlists_next[index] ] = rowlists_prev[ index ];

  if (rowlists_prev[ index ] != -1)
    rowlists_next[ rowlists_prev[index] ] = rowlists_next[ index ];

  rowlists_next[ index ] = rowlists_freehead;
  rowlists_freehead = index;
}

static int rowlists_isempty()
{
  int i,count;

  i = rowlists_freehead;
  count = 0;
  while (i != -1) {
    count++;
    i = rowlists_next[i];
  }

  if (count == rowlists_size) return 1;
  else return 0;
}
  
/****************************************************/
/*                                                  */
/* Spa routines                                     */
/*                                                  */
/****************************************************/

#if 0
static
int intcmp(const void* v1, const void* v2)
{
  const int* i1 = (const int*) v1;
  const int* i2 = (const int*) v2;

  if      (*i1 < *i2) { return -1; }
  else if (*i1 > *i2) { return  1; }
  else                  return  0;
}
#endif

static taucs_datatype*  spa;
static char*            spamap;

static void spa_finalize()
{
  taucs_free(spa);
  taucs_free(spamap);
}

static void spa_init(int nrows)
{
  int i;

  spa    = (taucs_datatype*) taucs_malloc(nrows*sizeof(taucs_datatype));
  spamap = (char*)  taucs_malloc(nrows*sizeof(char));
  assert(spa && spamap);

  for (i=0; i<nrows; i++) {spa[i] = taucs_zero; spamap[i] = 0;}
}

static void
gather(int             a_nnz,
       taucs_datatype* a_re,
       int*            a_ind,
       taucs_datatype* spa,
       char*           spamap)

{
  int i,ip;

#ifdef DETAILED_TIMING
  double time_tmp;

  gathers += ((double) a_nnz);
  time_tmp = taucs_wtime();
#endif

  for (ip=0; ip<a_nnz; ip++) {
    i = a_ind[ip];
    a_re[ip] = spa[i];
    spamap[i] = 0;

    spa[i] = taucs_zero;
  }

#ifdef DETAILED_TIMING
  time_gather += (taucs_wtime() - time_tmp);
#endif
}

static void 
scatter(int             a_nnz,
	taucs_datatype* a_re,
	int*            a_ind,
	taucs_datatype* spa,
	char*           spamap)
     
{
  int i,ip;

#ifdef DETAILED_TIMING
  double time_tmp;

  scatters += ((double) a_nnz);
  time_tmp = taucs_wtime();
#endif

  for (ip=0; ip<a_nnz; ip++) {
    i = a_ind[ip];
    spa[i]    = a_re[ip];
    spamap[i] = 1;
  }

#ifdef DETAILED_TIMING
  time_scatter += (taucs_wtime() - time_tmp);
#endif
}

/****************************************************/
/*                                                  */
/* Column updates                                   */
/*                                                  */
/****************************************************/

static
void spcol_spa_update(int pivotindex,
		      taucs_datatype* l_re,
		      int*    l_ind,
		      int     l_nnz,
		      int*    panelcols,
 		      int     panelcols_n,
		      int     nrows,
		      int**   a_inrowlist,
		      taucs_datatype* spa,
		      char*   spamap,
		      int**   a_ind,
		      int*    a_nnz)

{
  int i,ip,j,q;
  taucs_datatype v;

  for(j = 0; j < panelcols_n; j++) {
    q = panelcols[j];
    if (   taucs_re(spa[q*nrows + pivotindex]) == 0.0 
	&& taucs_im(spa[q*nrows + pivotindex]) == 0.0) continue;
    for (ip=0; ip<l_nnz; ip++) {
      i = l_ind[ip];
      v = l_re[ip];
      if (spamap[q*nrows + i] == 0) {
	spamap[q*nrows + i] = 1;
	spa   [q*nrows + i] = taucs_zero;
	a_ind      [q][ a_nnz[q] ] = i;
	a_inrowlist[q][ a_nnz[q] ] = rowlists_insert(i,q);
	(a_nnz[q])++;
      }
      /*spa[q*nrows + i] -= (spa[q*nrows + pivotindex] * v);*/
      spa[q*nrows + i] = taucs_sub(spa[q*nrows + i],
				   taucs_mul(spa[q*nrows + pivotindex] , v));
    }
  }
  return;
}

#ifdef SIMPLE_COL_COL
static void
spcol_spcol_update(int pivotindex,
		   taucs_datatype* l_re,
		   int*    l_ind,
		   int     l_nnz,
		   int     panelcol,
		   int*    a_inrowlist,
		   taucs_datatype* spa,
		   char*  spamap,
		   int*    a_ind,
		   int*    a_nnz,
		   int     lu_colcount)
{
  int i,ip;
  taucs_datatype pv;

  pv = spa[pivotindex];
  
  if (taucs_iszero(pv)) return;

#ifdef DETAILED_TIMING
  flops += 2.0 * ((double) l_nnz);
#endif /* DETAILED_TIMING */

  for (ip=0; ip<l_nnz; ip++) {
    i = l_ind[ip];
    if (spamap[i] == 0) {
      spamap[i] = 1;
      spa[i] = taucs_zero;
      a_ind[ *a_nnz ] = i;
      a_inrowlist[ *a_nnz ] = rowlists_insert(i,panelcol);
      (*a_nnz)++;
    }
    /*spa[i] -= (spa[pivotindex] * l_re[ip]);*/
    spa[i] = taucs_sub(spa[i],taucs_mul(spa[pivotindex] , l_re[ip]));
  }

  /*
  if ((*a_nnz) > lu_colcount) {
    taucs_printf("oocsp_numfact: prediction=%d, size now=%d\n",lu_colcount,*a_nnz);
    taucs_printf("oocsp_numfact: Column expands beyond prediction\n");
  }
  */

}
#else /* simple col col */

static int oocsp_spcol_n1 = 0;
static int oocsp_spcol_n2 = 0;

static void
spcol_spcol_update(int pivotindex,
			  taucs_datatype* l_re,
			  int*    l_ind,
			  int     l_nnz,
			  int     panelcol,
			  int*    a_inrowlist,
			  taucs_datatype* spa,
			  char*  spamap,
			  int*    a_ind,
			  int*    a_nnz,
			  int     lu_colcount)
{
  register int i;
  register int ip;
  register int ip_block;
  register int loop_bound;
  register char flag;
  register taucs_datatype pv;

  pv = spa[pivotindex];
  
  if (taucs_iszero(pv)) return;

#ifdef DETAILED_TIMING
  flops += 2.0 * ((double) l_nnz);
#endif /* DETAILED_TIMING */

  for (ip_block=0; ip_block<l_nnz; ip_block += BLOCK) {

    loop_bound = min(ip_block + BLOCK,l_nnz);

    flag = 1;
    oocsp_spcol_n1++;
    for (ip=ip_block; ip<loop_bound; ip++) {
      i = l_ind[ip];
      flag &= spamap[i];
      /* spa[i] -= (pv * l_re[ip]); */
      spa[i] = taucs_sub(spa[i], taucs_mul(pv , l_re[ip]));
    }

    if (!flag) {
      oocsp_spcol_n2++;

      for (ip=ip_block; ip<loop_bound; ip++) {
	i = l_ind[ip];
	if (spamap[i] == 0) {
	  spamap[i] = 1;
	  a_ind[ *a_nnz ] = i;
	  a_inrowlist[ *a_nnz ] = rowlists_insert(i,panelcol);
	  (*a_nnz)++;
	  /* we essentially zero and update */
	  /*spa[i] = - (pv * l_re[ip]);*/
	  spa[i] = taucs_neg(taucs_mul(pv , l_re[ip]));
	}
      }
    }
  }
  return;
}
#endif

#if 0
static void spcol_panel_update(int pivotindex,
			       taucs_datatype* l_re,
			       int*    l_ind,
			       int     l_nnz,
			       int*    subpanel, int subpanel_size, 
			       taucs_datatype* subpanel_tmp,
			       int**   a_inrowlist,
			       taucs_datatype* spa,
			       char*   spamap,
			       int**   a_ind,
			       int*    a_nnz,
			       int     nrows)
{
  int i,ii,ip,q,j;
  taucs_datatype x;

  assert(0);

#ifdef DETAILED_TIMING
  flops += ((double) subpanel_size) * 2.0 * ((double) l_nnz);
#endif /* DETAILED_TIMING */

  for (j=0; j<subpanel_size; j++) {
    q = subpanel[j];
    subpanel_tmp[q] = spa[q*nrows+pivotindex];
  }

  for (ip=0; ip<l_nnz; ip++) {
    i = l_ind[ip];
    x = l_re[ip];
    for (j=0; j<subpanel_size; j++) {
      q = subpanel[j];
      ii = (q*nrows) + i;
      if (spamap[ii] == 0) {
	spamap[ii] = 1;
	spa[ii] = taucs_zero;
	a_ind[q][ a_nnz[q] ] = i;
	a_inrowlist[q][ a_nnz[q] ] = rowlists_insert(i,q);
	(a_nnz[q])++;
      }
      /*spa[ii] -= (subpanel_tmp[q] * x);*/
      spa[ii] = taucs_sub(spa[ii], taucs_mul(subpanel_tmp[q] , x));
    }
  }
}
#endif

/****************************************************/
/*                                                  */
/* OLD STUFF                                        */
/*                                                  */
/****************************************************/

#if 0
void x_heap_insert( int* heap, int* heapsize,  int* ipivots, int value) 
{
  heap[ *heapsize ] = value;
  (*heapsize)++;

  taucs_printf("oocsp_numfact: heap insert %d\n",value);
}

int x_heap_extractmin( int* heap, int* heapsize, int* ipivots) 
{
  int i,m,mindex;

  if (*heapsize == 0) return -1;

  m = INT_MAX;

  for (i=0; i<(*heapsize); i++) 
    if (ipivots[heap[i]] < ipivots[m]) {m = heap[i]; mindex=i;}

  if (m == INT_MAX) return -1;

  for (i=mindex; i<(*heapsize)-1; i++) 
    heap[ i ] = heap[i+1];

  (*heapsize)--;
 
  taucs_printf("oocsp_numfact: heap extractmin %d (%d)\n",m,ipivots[m]);

  return m;
}
#endif

/****************************************************/
/*                                                  */
/* NUMERICAL FACTORIZATION MAIN ROUTINE             */
/*                                                  */
/****************************************************/

static
void oocsp_numfact (taucs_ccs_matrix* A, int* colperm,
		    taucs_io_handle* LU,
		    int* panels,
		    int* schedstart,
		    int* schedend,
		    int* fetchnext,
		    int* ejectnext,
		    int* lcolcount,
		    int* ucolcount,
		    int spawidth,
		    int maxsn
		    )
{
  int i,j,k,ip,p,q,qp,ii,ks;/* ip_next,jp_next,jp omer*/

  int nrows, ncols;

  
  /*char    fname[256]; 
  int     file; 
  int     mode; 
  mode_t  perm; omer*/

  /*int     len,status; omer*/

  /*
  double* pr;
  */

  /* NEW VARS (Sivan, for this function) */

  int fn;
  int en;
  int nsteps;
  int step;

  /*int nnz; omer*/
  /*off_t offset; omer*/
  /*ssize_t io_size; omer*/

  int  heapsize;
  int* heap;

  char* nnzmap;   /* bit vector */
  char* lindices; /* bit vector */
  /* char* uindices; */ /* bit vector */

  int* panel_id;
  int* panel_nnz;
  int**    panel_ind;
  int**    panel_inrowlist;
  taucs_datatype** panel_re;

  int* Lclen;
  int* Uclen;

  /*
  taucs_datatype* update_tmp;
  int* update_vec;
  int  update_vec_next;
  */

#ifdef SPA_ONEARRAY
  taucs_datatype*  panel_spa = NULL; /* warning */
  char*    panel_spamap = NULL; /* warning */
#else
  taucs_datatype**  panel_spa;
  char**    panel_spamap;
#endif

  int     unext,lnext;
  int*    lu_ind;
  taucs_datatype* lu_re;

  int* ipivots;

  int    maxind;
  double maxval, absval;
  int    pivotindex;

#ifdef SNODES
  int*    snode_ind;
  taucs_datatype* snode_re;
  /*char*   snode_map;*/ /* sivan changed to int to support m2 */
  int*    snode_map;
  int*    snode_pivrows;
  int*    snode_index;
  int*    pivots;

  int     snode_hash=0; /* warning */
  int     snode_size, snode_flag, snode_nnz, snode_lastcol;
  int     snode_last;
  int     hash;
  int     snode_id;

  taucs_datatype* S;
  taucs_datatype* P;
  int*    spa_updcols;
  int     spa_n;
  int*    m2;

  int*    srows;       /* indices of rows in the supernode */
  int     srows_n;     /* number of rows in the supernode  */
  int     srow_next;

  int     dense_flag;
  int     tmp;
#endif

  int maxcolcount;
  
  /*double time_tmp; omer*/

  /* READ GLOBALS */

  taucs_printf("oocsp_numfact: Using %.0lf MBytes of memory\n",
	     remaining_memory/1048576.0);

    
  /* START THE FACTORIZATION */
  /*
  nrows = A->nrows;
  ncols = A->ncols;
  */
  nrows = A->m;
  ncols = A->n;

#ifdef DETAILED_TIMING
  flops       = 0.0;
  scatters    = 0.0;
  gathers     = 0.0;
  num_heap_ops    = 0.0;
  rowlist_ops = 0.0;
  time_append = 0.0;
  time_read   = 0.0;
  time_colcol = 0.0;
  time_colcol_1 = 0.0;
  time_colcol_2 = 0.0;
  time_factor = 0.0;
  time_scatter= 0.0;
  time_gather = 0.0;
#ifdef SNODES
  time_snode_detect  = 0.0;
  time_snode_prepare = 0.0;
  time_snode_dense   = 0.0;
#endif
  bytes_read  = 0.0;
  bytes_appended = 0.0;
  col_read    = 0.0;
  col_ooc_updates = 0.0;
#endif /* DETAILED_TIMING */

  time_total  = taucs_wtime();
  
  maxcolcount = 0;
  for (j=0; j<ncols; j++)
    maxcolcount = max( maxcolcount, lcolcount[j]+ucolcount[j] );
  taucs_printf("oocsp_numfact: maxcolcount = %d, nrows = %d, spawidth = %d, maxsn = %d\n",
	     maxcolcount, nrows, spawidth, maxsn);

  Lclen            = (int*) taucs_calloc(ncols,sizeof(int));
  Uclen            = (int*) taucs_calloc(ncols,sizeof(int));
  assert(Uclen && Lclen);

  ipivots          = (int*)taucs_malloc(ncols*sizeof(int));
  assert(ipivots);
  for (i=0; i<nrows; i++) ipivots[i] = INT_MAX;

  /* create row lists */

  rowlists_init(maxcolcount * iabs(spawidth),nrows);
  
  fn = 0;
  en = 0;
  
  lindices = (char*)taucs_malloc(nrows*sizeof(char));
  /* uindices = (char*)taucs_malloc(nrows*sizeof(char)); */

  for (i=0; i<nrows; i++) {
    lindices[i] = 1;
    /*uindices[i] = 0;*/
  }

  nnzmap   = (char*)taucs_malloc(nrows*sizeof(char));
  heap     = (int*) taucs_malloc(nrows*sizeof(int));

#ifdef SNODES
  /* the supernodes consists of contiguous columns and we need to know */
  /* the corresponding pivot rows */
  pivots           = (int*)taucs_malloc(ncols*sizeof(int));
  assert(pivots);
  for (i=0; i<ncols; i++) pivots[i] = INT_MAX;

  snode_lastcol = snode_nnz = snode_size = 0;
  snode_pivrows = (int*)   taucs_malloc(maxsn*sizeof(int)); /* size was ncols */
  snode_ind     = (int*)   taucs_malloc(nrows*sizeof(int));
  snode_re      = (taucs_datatype*)taucs_malloc(nrows*sizeof(taucs_datatype));
  snode_map     = (int*)   taucs_malloc(nrows*sizeof(int));

  snode_id      = 0;
  snode_last    = 0;
  snode_index   = (int*)taucs_malloc(ncols*sizeof(int));
  for (i=0; i<ncols; i++) snode_index[i] = -1;
  for (i=0; i<nrows; i++) snode_map[i]   = -1;

  S = (taucs_datatype*) taucs_malloc( maxcolcount * maxsn  * sizeof(taucs_datatype) );
  P = (taucs_datatype*) taucs_malloc( maxcolcount * spawidth   * sizeof(taucs_datatype) );
  srows = (int*) taucs_malloc( maxcolcount * sizeof(int) );
  spa_updcols = (int*) taucs_malloc( spawidth * sizeof(int) );
  assert(spa_updcols);
  assert(srows);
  assert(S);
  assert(P);

  /*
  lu_re    = (taucs_datatype*)taucs_malloc(maxsn * nrows * sizeof(taucs_datatype) );
  lu_ind   = (int*)   taucs_malloc(maxsn * nrows * sizeof(int)    );
  */
  lu_re    = (taucs_datatype*)taucs_malloc(maxsn * maxcolcount * sizeof(taucs_datatype) );
  lu_ind   = (int*)   taucs_malloc(maxsn * maxcolcount * sizeof(int)    );
  taucs_printf("lu_re  = %08x -> %08x\n",lu_re,lu_re+(maxsn*maxcolcount));
  taucs_printf("lu_ind = %08x -> %08x\n",lu_ind,lu_ind+(maxsn*maxcolcount));
#else
  lu_re    = (taucs_datatype*)taucs_malloc(maxcolcount * sizeof(taucs_datatype));
  lu_ind   = (int*)   taucs_malloc(maxcolcount * sizeof(int));
  /*
  lu_re    = (taucs_datatype*)taucs_malloc(nrows*sizeof(taucs_datatype));
  lu_ind   = (int*)taucs_malloc(nrows*sizeof(int));
  */
#endif

  /* These two can be smaller */

  /*
  update_vec = (int*)taucs_malloc(ncols*sizeof(int));
  update_tmp = (taucs_datatype*)taucs_malloc(ncols*sizeof(taucs_datatype));
  */

  panel_id  = (int*)taucs_malloc(ncols*sizeof(int));
  panel_nnz = (int*)taucs_malloc(ncols*sizeof(int));
  panel_ind = (int**)taucs_malloc(ncols*sizeof(int*));
  panel_inrowlist = (int**)taucs_malloc(ncols*sizeof(int*));
  panel_re  = (taucs_datatype**)taucs_malloc(ncols*sizeof(taucs_datatype*));
  for (i=0; i<ncols; i++) panel_id[i] = -1;

  if (spawidth > 0) {
#ifdef SPA_ONEARRAY
    panel_spa = (taucs_datatype*)taucs_malloc(spawidth*nrows*sizeof(taucs_datatype));
    panel_spamap = (char*)taucs_malloc(spawidth*nrows*sizeof(char) );
    assert(panel_spa && panel_spamap);
#else
    panel_spa = (taucs_datatype**) taucs_malloc(spawidth*sizeof(taucs_datatype*));
    panel_spamap = (char**) taucs_malloc(spawidth*sizeof(char*) );
    assert(panel_spa && panel_spamap);
#endif
  } 
    
  spa_init( nrows );

  nsteps=0; 
  for (i=0; i<ncols; i++) {
    if ((int) (schedstart[i]) > nsteps) nsteps = (int) schedstart[i];
  }


  /* INITIALIZE HEAP */

  heapsize = 0;
  for (i=0; i<nrows; i++) nnzmap[i] = 0;

  taucs_printf("oocsp_numfact: Starting numerical factorization (%d steps)\n",nsteps);

  p = 0;

  for (step=1; step<=nsteps; step++) {

    taucs_printf("oocsp_numfact: Starting step %d/%d\r",step,nsteps);

    if (p==0) {
      /*taucs_printf("oocsp_numfact: (new panel)\n");*/
#ifdef SNODES
      snode_last = -1;
#endif

      if (!rowlists_isempty())
	taucs_printf("oocsp_numfact: Internal Error (row lists not empty)\n");
      if (heapsize) 
	taucs_printf("oocsp_numfact: Internal Error (heap not empty; 1)\n");
      for (i=0; i<nrows; i++) 
	if (nnzmap[i]) taucs_printf("oocsp_numfact: Internal Error (heap not empty; 2)\n");
    } else {
      taucs_printf("oocsp_numfact: (same panel)\n");
    }
    
    /* LOAD A PARTIAL PANEL */
    
    while (fn < ncols && 
            (schedstart[fetchnext[fn]]) == step) {
      j = fetchnext[fn];

      /*
      if (A->clen[j] > lcolcount[j]+ucolcount[j]) {
	taucs_printf("oocsp_numfact: Column %d in matrix larger than L+U estimate\n",j);
	taucs_printf("oocsp_numfact: Size %d > %d+%d\n",A->clen[j],lcolcount[j],ucolcount[j]); 
	taucs_printf("oocsp_numfact: Aborting\n");
      }

      panel_nnz[p] = A->clen[j];
      */

      /* New code: Sivan 28 Feb 2002 */

      panel_nnz[p] = (A->colptr)[colperm[j]+1] - (A->colptr)[colperm[j]];
      panel_id[p]  = j;
      if (lcolcount[j] < 1 
	  || lcolcount[j] > nrows
	  || ucolcount[j] < 1 
	  || ucolcount[j] > nrows) {
	taucs_printf("oocsp_numfact: Column %d: l,u colcounts = (%d+%d)\n",
		  j,lcolcount[j],ucolcount[j]);
	taucs_printf("oocsp_numfact: Invalid column nonzero count\n\n");
      }	

      panel_inrowlist[p] = (int*)    taucs_malloc((lcolcount[j] + ucolcount[j]) * sizeof(int));
      panel_ind[p]       = (int*)    taucs_malloc((lcolcount[j] + ucolcount[j]) * sizeof(int));
      panel_re[p]        = (taucs_datatype*) taucs_malloc((lcolcount[j] + ucolcount[j]) * sizeof(taucs_datatype));
      if (!panel_inrowlist[p] || 
	  !panel_ind[p] || 
	  !panel_re[p]) {
	fprintf(stderr,"out of memory for panel compressed vector\n");
	fprintf(stderr,"j = %d lcolcount = %d ucolcount = %d\n",
		j,lcolcount[j],ucolcount[j]);
	exit(1);
      }
      assert (panel_inrowlist[p] && panel_ind[p] && panel_re[p]);

      if (spawidth > 0) {
#ifdef SPA_ONEARRAY
#else
	panel_spa[p]     = (taucs_datatype*) taucs_malloc(nrows * sizeof(taucs_datatype));
	panel_spamap[p]  = (char*)   taucs_malloc(nrows * sizeof(char));
	assert(panel_spa[p] && panel_spamap[p]);
#endif
      }

      /*
      time_tmp = taucs_wtime();
      oocsp_readcol(A,j,panel_ind[p],panel_re[p]);
      time_read += (taucs_wtime() - time_tmp);
      bytes_read += (double) (panel_nnz[p] * (sizeof(taucs_datatype)+sizeof(int)));
      */

      /* new code: Sivan 28 Feb 2002 */
      {
	int jp,ip,i;
	jp = colperm[j];
	for (i=0, ip=(A->colptr)[jp]; ip<(A->colptr)[jp+1]; i++, ip++) {
	  (panel_ind[p])[i] = (A->rowind)[ip];
	  (panel_re [p])[i] = (A->taucs_values)[ip];
	  /*
	  printf(">>> reading column number %d (index %d), row %d value %.2e\n",
		 j,jp,(panel_ind[p])[i],(panel_re[p])[i]);
	  */
	}
      }

      if (spawidth > 0) {
	if (p >= spawidth) {
	  taucs_printf("oocsp_numfact: p=%d spawidth=%lg\n",p,spawidth);
	  taucs_printf("oocsp_numfact: Panel wider than spawidth\n");
	} 
	scatter(panel_nnz[p],panel_re[p],panel_ind[p],
#ifdef SPA_ONEARRAY
		panel_spa+(p * nrows),
		panel_spamap+(p * nrows)
#else
		panel_spa[p],
		panel_spamap[p]
#endif
		);
      }

      for (i=0; i<panel_nnz[p]; i++) {
	panel_inrowlist[p][i] = rowlists_insert(panel_ind[p][i],
						p);

	/*if (nnzmap[ panel_ind[p][i] ] == 0 && uindices[ panel_ind[p][i] ]) {*/
	if (nnzmap[ panel_ind[p][i] ] == 0 && !lindices[ panel_ind[p][i] ]) {
	  nnzmap[ panel_ind[p][i] ] = 1;
	  num_heap_insert( heap, &heapsize, ipivots, panel_ind[p][i] );
	}
      }

      p++;
      fn++;
    }

    /*taucs_printf("oocsp_numfact: Updating (step %d)\n",step);*/

#ifndef SNODES
    /* we need to prevent panel cols from being loaded! */

    while ((i = num_heap_extractmin( heap, &heapsize, ipivots)) != -1) {
      if (i == INT_MAX) continue;

      nnzmap[ i ] = 0;
      
      k = ipivots[i];
      if (k == INT_MAX) continue;
      
      /*taucs_printf("oocsp_numfact: extracted row %d (col %d) from heap. %d\n",i,k,uindices[i]);*/

      if (panels[k] == panels[ panel_id[0] ]) {
	/*taucs_printf("oocsp_numfact: skipping column from this panel (col %d, panel %lg)\n",
	       k,panels[k]);*/
	continue;
      }
      
      /*if (lindices[i]) continue;*/   /*if (!uindices[i]) continue;*/

#ifdef DETAILED_TIMING
      time_tmp = taucs_wtime();
#endif
      /* Sivan: replaced 2 March 2002 */
      /*oocsp_readcol(L,k,lu_ind,lu_re);*/
      Lreadcol(LU,k,Lclen[k],lu_ind,lu_re);

#ifdef DETAILED_TIMING
      bytes_read += (double) (Lclen[k] * (sizeof(taucs_datatype)+sizeof(int)));
      col_read += 1.0;
#endif /* DETAILED_TIMING */

      /*
      for (ks = k+1; ks<ncols && snode_index[ks]==snode_index[k]; ks++) {
	oocsp_readcol(L,ks,lu_ind+((ks-k)*nrows),lu_re+((ks-k)*nrows));
	bytes_read += (double) ((L.clen[ks]) * (sizeof(taucs_datatype)+sizeof(int)));
#ifdef DETAILED_TIMING
	col_read += 1.0;
#endif
      }
      taucs_printf("oocsp_numfact: Read supernode, %d cols\n",ks-k+1);
      */

#ifdef DETAILED_TIMING
      time_read += (taucs_wtime() - time_tmp);
#endif


      /*
      taucs_printf("oocsp_numfact: col %d (panel %lg, now %lg; read from disk) for row %d, length %d\n",
	     k,panels[k],panels[panel_id[0]],i,(L.clen)[k]);
	     */

      if (Lclen[k] < 0) {
	taucs_printf("oocsp_numfact: k=%d len=%d\n",k,Lclen[k]);
	taucs_printf("oocsp_numfact: Internal Error\n");
      }

      for (ii=0; ii<Lclen[k]; ii++) {
	/*if (nnzmap[ lu_ind[ii] ] == 0 && uindices[ lu_ind[ii] ]*/
	if (nnzmap[ lu_ind[ii] ] == 0 && !lindices[ lu_ind[ii] ]
	    && !taucs_iszero(lu_re[ii])) {
	  nnzmap[ lu_ind[ii] ] = 1;
	  num_heap_insert( heap, &heapsize, ipivots, lu_ind[ii] );
	  /*taucs_printf("oocsp_numfact: inserting row %d into heap\n",lu_ind[ii]);*/
	}
      }

      for(qp = rowlists_head[i]; qp != -1; qp = rowlists_next[qp]) {
	q = rowlists_colind[qp];
#ifdef DETAILED_TIMING
	time_tmp = taucs_wtime();
#endif
	if (spawidth <= 0) {
	  scatter(panel_nnz[q],panel_re[q],panel_ind[q],
		  spa,spamap);
	  spcol_spcol_update(i,
			     lu_re,lu_ind,Lclen[k],
			     q,panel_inrowlist[q],
			     spa,
			     spamap,
			     panel_ind[q],&(panel_nnz[q]),
			     lcolcount[panel_id[q]]+ucolcount[panel_id[q]]);
	  gather(panel_nnz[q],panel_re[q],panel_ind[q],
		 spa,spamap);
	} else {
	  spcol_spcol_update(i,
			     lu_re,lu_ind,Lclen[k],
			     q,panel_inrowlist[q],
#ifdef SPA_ONEARRAY
			     panel_spa+(q * nrows),
			     panel_spamap+(q * nrows),
#else
			     panel_spa[q],
			     panel_spamap[q],
#endif
			     panel_ind[q],&(panel_nnz[q]),
			     lcolcount[panel_id[q]]+ucolcount[panel_id[q]]);
	}
#ifdef DETAILED_TIMING
	time_colcol += (taucs_wtime() - time_tmp);
	col_ooc_updates += 1.0;
#endif
      }

      /*
      update_vec_next = 0;
      for(qp = rowlists_head[i]; qp != -1; qp = rowlists_next[qp]) {
	q = rowlists_colind[qp];
	col_ooc_updates += 1.0;
	update_vec[update_vec_next] = q;
	update_vec_next++;
      }

#ifdef DETAILED_TIMING
      time_tmp = taucs_wtime();
#endif

#ifdef SPA_ONEARRAY
      spcol_panel_update(i,
			 lu_re,lu_ind,Lclen[k],
			 update_vec,update_vec_next,update_tmp,
			 panel_inrowlist,
			 panel_spa,
			 panel_spamap,
			 panel_ind,panel_nnz,
			 nrows);
#endif

#ifdef DETAILED_TIMING
      time_colcol += (taucs_wtime() - time_tmp);
#endif
      */
    }
#else  /* with SNODES */
    /* we need to prevent panel cols from being loaded! */

    while ((i = num_heap_extractmin( heap, &heapsize, ipivots)) != -1) {
      if (i == INT_MAX) continue;

      nnzmap[ i ] = 0;
      
      k = ipivots[i];
      if (k == INT_MAX) continue;
      
      /*printf("last supernode %d this column is %d its supernode %d\n",snode_last,k,snode_index[k]);*/
      if (snode_index[k] == snode_last) continue; /* skip rest of supernode */
      snode_last = snode_index[k]; /* mark for next time */
      
      /*taucs_printf("oocsp_numfact: extracted row %d (col %d) from heap. %d\n",i,k,uindices[i]);*/

      if (panels[k] == panels[ panel_id[0] ]) {
	/*taucs_printf("oocsp_numfact: skipping column from this panel (col %d, panel %lg)\n",
	  k,panels[k]);*/
	continue;
      }
      
      /*if (lindices[i]) continue;*/ /*if (!uindices[i]) continue;*/

#ifdef DETAILED_TIMING
      time_tmp = taucs_wtime();
#endif

      for (ks = k; ks<ncols && snode_index[ks]==snode_index[k]; ks++) {
	/*printf(">>> %d %d\n",ks,snode_index[ks]);*/
	/*oocsp_readcol(L,ks,lu_ind+((ks-k)*maxcolcount),lu_re+((ks-k)*maxcolcount));*/
	Lreadcol(LU,ks,Lclen[ks],lu_ind+((ks-k)*maxcolcount),lu_re+((ks-k)*maxcolcount));
#ifdef DETAILED_TIMING
	bytes_read += (double) (Lclen[ks] * (sizeof(taucs_datatype)+sizeof(int)));
	col_read += 1.0;
#endif /* DETAILED_TIMING */
      }

#ifdef DETAILED_TIMING
      time_read += (taucs_wtime() - time_tmp);
#endif
      /*taucs_printf("oocsp_numfact: Read supernode, %d cols %d:%d\n",ks-k,k,ks-1);*/
      
      if (Lclen[k] < 0) {
	taucs_printf("oocsp_numfact: k=%d len=%d\n",k,Lclen[k]);
	taucs_printf("oocsp_numfact: Internal Error\n");
      }

      for (ii=0; ii<Lclen[k]; ii++) {
	/*if (nnzmap[ lu_ind[ii] ] == 0 && uindices[ lu_ind[ii] ]*/
	if (nnzmap[ lu_ind[ii] ] == 0 && !lindices[ lu_ind[ii] ]
	    && !taucs_iszero(lu_re[ii])) {
	  nnzmap[ lu_ind[ii] ] = 1;
	  num_heap_insert( heap, &heapsize, ipivots, lu_ind[ii] );
	  /*taucs_printf("oocsp_numfact: inserting row %d into heap\n",lu_ind[ii]);*/
	}
      }

      /* determine which panel cols need to be updates */
      dense_flag = 1;
      for (tmp = 0; tmp<1; tmp++) {
	int jj,ii,kk,iip,jjp;

	/*	printf("snode %d:%d\n",k,ks-1);*/

	if (ks-k < SNODE_THRESHOLD) {dense_flag=0; break;}

	assert(spawidth > 0); /* the code for the other case is not implemented yet; Sivan */

#ifdef DETAILED_TIMING
	time_tmp = taucs_wtime();
#endif

	srows_n = Lclen[k] + 1; /* the diagonal element in column k is not */
                                /* represented explicitely in L, it's 1    */

	if (srows_n < SNODE_THRESHOLD) {dense_flag=0; break;}

	m2 = snode_map; /* reuse the vector */

	/* create an array of row indices, with pivot rows first, in order */
	/* so that the supernode array will be a trapezoidal matrix        */
	/* we should keep this and restore the -1 invariant                */

#ifdef DETAILED_TIMING
	time_snode_tmp = taucs_wtime();
#endif

	for (jj=k; jj<ks; jj++) {
	  ii = pivots[jj];
	  srows[jj-k] = ii;
	  m2[ii] = 1;
	}
	/* now the rest of the rows */
	srow_next = ks-k;
	for (iip=0; iip<Lclen[k]; iip++) {
	  ii = lu_ind[iip];
	  if (ii == pivots[k]) srows_n--; /* I have no idea why this row shows up; always with 0.0 */
	  if (m2[ii] == -1) {
	    srows[srow_next] = ii;
	    srow_next++;
	  } else {
	    m2[ii] = -1; /* restore the invariant */
	  }
	}
	assert(srow_next == srows_n);
	
	/* we begin by figuring out which columns of the panel are updated    */
	/* by this supernode.                                                 */

#ifdef DETAILED_TIMING
	time_snode_tmp = taucs_wtime();
#endif

	spa_n = 0;
	for (jj=k; jj<ks; jj++) {
	  ii = pivots[jj];
	  for(qp = rowlists_head[ii]; qp != -1; qp = rowlists_next[qp]) {
	    int skip; /* don't add a column twice to spa_updcols */
	    q = rowlists_colind[qp];
	    for (skip=0, jjp=0; jjp<spa_n; jjp++) {
	      if (spa_updcols[jjp] == q) {
		skip = 1;
		break;
	      }
	    }
	    if (!skip) {
	      /*if (jj-k > 4) printf("*** jj-k %d ks-k %d\n",jj-k,ks-k);*/

#ifdef DETAILED_TIMING
	      flops += 2.0 * ( (ks-jj) * (srows_n - (jj-k)) - 0.5*(ks-jj)*(ks-jj) );
	      flops_extra += 2.0 * ( (jj-k) * (srows_n) - 0.5*(jj-k)*(jj-k) );
#endif /* DETAILED_TIMING */

	      spa_updcols[spa_n] = q;
	      /*spa_updptrs[spa_n] = jj-k;*/
	      spa_n++;
	    }
	  }
	}
#ifdef DETAILED_TIMING
	time_snode_1 += (taucs_wtime()-time_snode_tmp);
#endif

	if (spa_n < SNODE_THRESHOLD) {
	  for (jj=k; jj<ks; jj++) {
	    ii = pivots[jj];
	    m2[ii] = -1;
	  }
	  for (iip=0; iip<Lclen[k]; iip++) {
	    ii = lu_ind[iip];
	    m2[ii] = -1; /* restore the invariant */
	  }
	  dense_flag=0; 
	  break;
	}

	/*printf("snode %d:%d (dense)\n",k,ks-1);*/

	/* copy snode columns into a dense array */

	for (jj=k; jj<ks; jj++) {
	  for (iip=0; iip<Lclen[jj]; iip++)
	    /* mark in m2 where each row is */
	    m2[ lu_ind[(jj-k)*maxcolcount + iip] ] = iip;
	  
	  for (iip=0; iip<srows_n; iip++) {
	    if (m2[ srows[iip] ] != -1) {
	      S[ (jj-k)*srows_n + iip ] = lu_re[(jj-k)*maxcolcount + m2[ srows[iip] ]];
	      m2[ srows[iip] ] = -1; /* restore invariant for next column */
	    } else
	      S[ (jj-k)*srows_n + iip ] = taucs_zero;
	  }
	}
#ifdef DETAILED_TIMING
	time_snode_3 += (taucs_wtime()-time_snode_tmp);
#endif

	/* now the snode is stored in a dense array S, with row indices srows */
	/* and with the diagonal block of L on top.                           */

	/* next we copy the columns of the panels that need to be updated     */
	/* into the dense array P.                                            */

	/* we then copy these columns into the dense array P, and if          */
	/* fill occurs, we update the nonzero bitmap and row lists.           */

#ifdef DETAILED_TIMING
	time_snode_tmp = taucs_wtime();
#endif
#define OLD_1
#ifdef OLD_1
	for (jjp=0; jjp<spa_n; jjp++) {
	  jj = spa_updcols[jjp];

	  for (iip=0; iip<srows_n; iip++) {
	  /*for (iip=spa_updptrs[jjp]; iip<srows_n; iip++) {*/
	    if (panel_spamap[jj*nrows + srows[iip]] == 0) {
	      /* fill will occur here */
	      panel_spa[jj*nrows + srows[iip]] = taucs_zero;
	      /*P[jjp*srows_n + iip] = taucs_zero;*/
	      
	      panel_spamap   [jj*nrows + srows[iip]] = 1;
	      panel_ind      [jj][panel_nnz[jj] ] = srows[iip];
	      panel_inrowlist[jj][panel_nnz[jj] ] = rowlists_insert(srows[iip],jj);
	      (panel_nnz[jj])++;
	    }
#ifdef JUNK
	    } else { 
	      /* no fill, just copy the element */
  	      P[jjp*srows_n + iip] = panel_spa[jj*nrows + srows[iip]];
	    }
#endif
	  }
	} 

	{ 
#ifdef DETAILED_TIMING
	  double x = taucs_wtime();
#endif
	  for (jjp=0; jjp<spa_n; jjp++) {
	    jj = spa_updcols[jjp];
	    for (iip=0; iip<srows_n; iip++) {
	      P[jjp*srows_n + iip] = panel_spa[jj*nrows + srows[iip]];
	    }
	  }
#ifdef DETAILED_TIMING
	  time_snode_21 += (taucs_wtime()-x);
#endif
        }
#else
	for (jjp=0; jjp<spa_n; jjp++) {
	  int iip_block, loop_bound, flag;
	  
	  jj = spa_updcols[jjp];

	  for (iip_block=0; iip_block<srows_n; iip_block += SNODE_BLOCK) {

	    loop_bound = min(iip_block + SNODE_BLOCK,srows_n);
	    
	    flag = 1;
	    for (iip=iip_block; iip<loop_bound; iip++)
	      flag &= panel_spamap[jj*nrows + srows[iip]];
	     
	    if (!flag) {
	      for (iip=iip_block; iip<loop_bound; iip++) {
		if (panel_spamap[jj*nrows + srows[iip]] == 0) {
		  panel_spamap   [jj*nrows + srows[iip]] = 1;
		  panel_spa      [jj*nrows + srows[iip]] = taucs_zero;
		  panel_ind      [jj][panel_nnz[jj] ] = srows[iip];
		  panel_inrowlist[jj][panel_nnz[jj] ] = rowlists_insert(srows[iip],jj);
		  (panel_nnz[jj])++;
		}
	      }
	    }
	  }
	  
	  { 
#ifdef DETAILED_TIMING
	    double x = taucs_wtime();
#endif
	    for (iip=0; iip<srows_n; iip++)
	      P[jjp*srows_n + iip] = panel_spa[jj*nrows + srows[iip]];
#ifdef DETAILED_TIMING
	    time_snode_21 += (taucs_wtime()-x);
#endif
	  }
	}
#endif

#ifdef DETAILED_TIMING
	time_snode_2 += (taucs_wtime()-time_snode_tmp);
#endif

	/*printf("supernode update: col %d pivotrow %d updates col %d\n",jj,ii,panel_id[q]);*/

#ifdef DETAILED_TIMING
	time_snode_prepare += (taucs_wtime() - time_tmp);
#endif

#ifdef DETAILED_TIMING
	time_tmp = taucs_wtime();
#endif

	/*flops += (2.0 * spa_n * srows_n * (ks-k) - 2.0); */ /* over estimate; sivan. */
	/* we can subract triangle in estimate, skip zero pivots in flops & code */
#ifdef DETAILED_TIMING
	flops_dense += (2.0 * spa_n * srows_n * (ks-k) -
			1.0 * spa_n * (ks-k) * (ks-k)); 
#endif /* DETAILED_TIMING */

#ifdef JUNK
	printf("TRSM's RHS:\n");
	for (kk=k; kk<ks; kk++) {
	  for (jj=0; jj<spa_n; jj++) {
	    printf("%c", P[jj*srows_n + (kk-k)] ? '*' : 'o');
	  }
	  printf("\n");
	}
#endif

#ifdef USE_BLAS
	{ 
	  int m = ks-k;
	  int n = spa_n;
	  int M = srows_n - m;
	  int N = n;
	  int K = m;
	  
	  if (m > BLAS_THRESHOLD && n > BLAS_THRESHOLD) {
	    taucs_trsm("Left",
		       "Lower",
		       "No transpose",
		       "Unit",
		       &m,&n,
		       &taucs_one_const,
		       S,&srows_n,
		       P,&srows_n
		       );
	  } else {
	    /* TRSM */
	    for (jj=0; jj<spa_n; jj++) {
	      for (kk=k; kk<ks; kk++) {
		for (ii=kk-k+1; ii<(ks-k); ii++) {
		  /*P[jj*srows_n + ii] -= (P[jj*srows_n + (kk-k)] * S[(kk-k)*srows_n + ii]);*/
		  P[jj*srows_n + ii] = 
		  taucs_sub(P[jj*srows_n + ii],
			    taucs_mul(P[jj*srows_n + (kk-k)] , S[(kk-k)*srows_n + ii]));
		}
	      }
	    }
	  }

	  if (M > BLAS_THRESHOLD && N > BLAS_THRESHOLD && K > BLAS_THRESHOLD) {
	    taucs_gemm("No transpose",
		       "No transpose",
		       &M,&N,&K,
		       &taucs_minusone_const,
		       S+m, &srows_n,
		       P  , &srows_n,
		       &taucs_one_const,
		       P+m, &srows_n);
	  } else {
	    /* GEMM */
	    for (jj=0; jj<spa_n; jj++) {
	      for (kk=k; kk<ks; kk++) {
		for (ii=(ks-k); ii<srows_n; ii++) {
		  /*P[jj*srows_n + ii] -= (P[jj*srows_n + (kk-k)] * S[(kk-k)*srows_n + ii]);*/
		  P[jj*srows_n + ii] =
		    taucs_sub(P[jj*srows_n + ii],
			      taucs_mul(P[jj*srows_n + (kk-k)] , S[(kk-k)*srows_n + ii]));
		}
	      }
	    }
	  }
	}
#else	   
	/* TRSM */
	for (jj=0; jj<spa_n; jj++) {
	  for (kk=k; kk<ks; kk++) {
	    for (ii=kk-k+1; ii<(ks-k); ii++) {
	      P[jj*srows_n + ii] -= (P[jj*srows_n + (kk-k)] * S[(kk-k)*srows_n + ii]);
	    }
	  }
	}

	/* GEMM */
	for (jj=0; jj<spa_n; jj++) {
	  for (kk=k; kk<ks; kk++) {
	    for (ii=(ks-k); ii<srows_n; ii++) {
	      P[jj*srows_n + ii] -= (P[jj*srows_n + (kk-k)] * S[(kk-k)*srows_n + ii]);
	    }
	  }
	}
#endif
#ifdef DETAILED_TIMING
	time_snode_dense += (taucs_wtime() - time_tmp);
	col_ooc_updates += 1.0;
#endif

#ifdef DETAILED_TIMING
	time_tmp = taucs_wtime();
#endif

	/* now copy panel columns out of the dense P */

#ifdef DETAILED_TIMING
	time_snode_tmp = taucs_wtime();
#endif
	for (jjp=0; jjp<spa_n; jjp++) {
	  jj = spa_updcols[jjp];
	  for (iip=0; iip<srows_n; iip++) {
	    /*for (iip=spa_updptrs[jjp]; iip<srows_n; iip++) {*/
	    panel_spa[jj*nrows + srows[iip] ] = P[jjp*srows_n + iip];
	  }
	}
#ifdef DETAILED_TIMING
	time_snode_4 += (taucs_wtime()-time_snode_tmp);
        time_snode_prepare += (taucs_wtime() - time_tmp);
#endif
      }
      if (!dense_flag) { /* we didn't do it using the blas since m,n, or k were too small */
	/* not worth copying into dense arrays etc */
	/* determine which panel cols need to be updates */
	int jj,ii;
	for (jj=k; jj<ks; jj++) {
	  ii = pivots[jj];
	  if (ii == INT_MAX) {
	    taucs_printf("oocsp_numfact: internal error (supernode update)\n");
	    exit(1);
	  }
	  for(qp = rowlists_head[ii]; qp != -1; qp = rowlists_next[qp]) {
	    q = rowlists_colind[qp];
	    /*printf("supernode update: col %d pivotrow %d updates col %d\n",jj,ii,panel_id[q]);*/

#ifdef DETAILED_TIMING
	    time_tmp = taucs_wtime();
#endif
	    if (spawidth <= 0) {
	      taucs_printf("oocsp_numfact: internal error (supernode without a spa)\n");
	      exit(1);
	    } else {
	      spcol_spcol_update(ii,
				 lu_re+((jj-k)*maxcolcount),lu_ind+((jj-k)*maxcolcount),Lclen[jj],
				 q,panel_inrowlist[q],
#ifdef SPA_ONEARRAY
				 panel_spa+(q * nrows),
				 panel_spamap+(q * nrows),
#else
				 panel_spa[q],
				 panel_spamap[q],
#endif
				 panel_ind[q],&(panel_nnz[q]),
				 lcolcount[panel_id[q]]+ucolcount[panel_id[q]]);
	    }
#ifdef DETAILED_TIMING
	    time_colcol   += (taucs_wtime() - time_tmp);
	    time_colcol_1 += (taucs_wtime() - time_tmp);
	    col_ooc_updates += 1.0;
#endif
	  }
	}
      }
    }

#endif /* else SNODES */

    
    /*taucs_printf("oocsp_numfact: Factoring and Writing (step %d)\n",step);*/

    while (en < ncols && 
            (schedend[ejectnext[en]]) == step) {
      p--;
      j = ejectnext[en];

#ifdef DETAILED_TIMING
      time_tmp = taucs_wtime();
#endif
      
      if (panel_id[p] < 0)  taucs_printf("oocsp_numfact: internal error (panel stack)\n");
      if (panel_id[p] != j) {
	taucs_printf("oocsp_numfact: en=%d p=%d, j=%d, panel_id[p]=%d\n",
	       en,p,j,panel_id[p]);
	taucs_printf("oocsp_numfact: internal error (panel id's)\n");
      }

      /* gather */

      if (spawidth > 0) {
	gather(panel_nnz[p],panel_re[p],panel_ind[p],
#ifdef SPA_ONEARRAY
	       panel_spa + (p * nrows),
	       panel_spamap + (p * nrows));
#else
	       panel_spa[p],
	       panel_spamap[p]);
#endif
      }

      /* Find pivot */
      
      maxval = 0;
      maxind = -1;
      for (ip = 0; ip<panel_nnz[p]; ip++) {
	absval = (double) taucs_abs(panel_re[p][ip]);
	if (!lindices[panel_ind[p][ip]]) continue;
	if (absval > maxval) {
	  maxval = absval;
	  maxind = ip;
	}
      }
      
      if ( maxind == -1 ) {
	taucs_printf("oocsp_numfact: Zero Column!\n");
      }
      
      if ( taucs_iszero(panel_re[p][maxind]) ) {
	taucs_printf("oocsp_numfact: Zero Pivot!\n");
      }
      
      pivotindex = panel_ind[p][maxind];

      /*
      taucs_printf("oocsp_numfact: pivot for column %d is %d (%lg)\n",j,pivotindex,
	     panel_re[p][maxind]);
	     */

      if ( ipivots[ pivotindex ] != INT_MAX ) 
	taucs_printf("oocsp_numfact: Pivoting twice on the same row\n");
      
      ipivots[ pivotindex ] = j;
#ifdef SNODES
      pivots[ j ] = pivotindex;
#endif
      /*printf("### pivot row for column %d is %d\n",j,pivotindex);*/

      lindices[ pivotindex ] = 0;
      /*uindices[ pivotindex ] = 1;*/
      
      /* copy to L, U */
      
      /*taucs_printf("oocsp_numfact: copying to l,u\n");*/
      
      lnext = 0;
      /* unext = nrows-1; */
      unext = maxcolcount-1;
      for (ip = 0; ip<panel_nnz[p]; ip++) {
	i = panel_ind[p][ip];
	/*if (uindices[i]) {*/
	if (!lindices[i]) {
	  lu_re[unext]  = panel_re[p][ip];
	  lu_ind[unext] = i;
	  /*assert( unext > 0 && unext < maxsn*maxcolcount);*/
	  unext--;
	} else {
	  /*lu_re[lnext]  = (panel_re[p][ip])/(panel_re[p][maxind]);*/
	  lu_re[lnext]  = taucs_div(panel_re[p][ip] , panel_re[p][maxind]);
	  lu_ind[lnext] = i;
	  lnext++;
	}	  
      }

#ifdef DETAILED_TIMING
      flops += (double) lnext;
#endif /* DETAILED_TIMING */
      
#ifdef DETAILED_TIMING
      time_factor += (taucs_wtime() - time_tmp);
#endif

      /* Write out column of L, U */
      /*
      taucs_printf("oocsp_numfact: writing column %d of l,u (sizes %d %d)\n",
	     j,lnext,nrows - 1 - unext);
	     */

#ifdef DETAILED_TIMING
      time_tmp = taucs_wtime();
#endif
      /*
      oocsp_appendcol(U,j,nrows - 1 - unext,
		      lu_ind + (unext+1),
		      lu_re  + (unext+1));
      */

      /* Sivan: replaced this code on 2 March 2002 */
      /*
      oocsp_appendcol(U,j,maxcolcount - 1 - unext,
		      lu_ind + (unext+1),
		      lu_re  + (unext+1));
      */
      assert(Uclen[j] == 0);
      Uclen[j] = maxcolcount - 1 - unext;
      Uappendcol(LU,j,maxcolcount - 1 - unext,
		 lu_ind + (unext+1),
		 lu_re  + (unext+1));
#ifdef DETAILED_TIMING
      time_append += (taucs_wtime() - time_tmp);
#endif

#ifdef DETAILED_TIMING
      bytes_appended += (double) ((maxcolcount - 1 - unext) 
				  * (sizeof(taucs_datatype)+sizeof(int)));
#endif /* DETAILED_TIMING */

#ifdef SNODES
      /* detect supernodes */

#ifdef DETAILED_TIMING
      time_tmp = taucs_wtime();
#endif
      
      snode_flag = 1; /* we assume so for now */

      /* a quick check using a hash function to quickly rule out columns */

      if (j>0 && panels[j] != panels[j-1]) snode_flag = 0;

      hash = pivotindex;
      for (ii=0; ii<snode_size; ii++)
	hash ^= snode_pivrows[ii];
      for (ii=0; ii<lnext; ii++)
	hash ^= lu_ind[ii];

      if (hash == snode_hash) {
	/* The hash is identical, but is it in the same supernode? */

	/* first, mark this column's structure in the snode_map    */
	/* then make sure all the supernode's nonzeros are marked  */

	/* for (ii=0; ii<nrows; ii++) assert( snode_map[ii] == -1 ); */

	                                snode_map[pivotindex       ] = 1;
	for (ii=0; ii<lnext; ii++)      snode_map[lu_ind[ii]       ] = 1;
	for (ii=0; ii<snode_size; ii++) snode_map[snode_pivrows[ii]] = 1;
	for (ii=0; ii<snode_nnz; ii++)
	  if (snode_map[snode_ind[ii]] != 1) snode_flag = 0;

	/* next, zero this column's structure in the snode_map     */
	/* then make sure all the supernode's nonzeros are marked  */


    	                                snode_map[pivotindex]        = -1;
	for (ii=0; ii<lnext; ii++)	snode_map[lu_ind[ii]       ] = -1;
	for (ii=0; ii<snode_size; ii++) snode_map[snode_pivrows[ii]] = -1;
	
	/* mark this supernodes nonzeros in the map                */
	
	for (ii=0; ii<snode_nnz; ii++) snode_map[snode_ind[ii]] = 1;

	/* check that all the column's nonzeros are in the snode   */

	for (ii=0; ii<lnext; ii++) if (snode_map[lu_ind[ii]] != 1) snode_flag = 0;

	/* zero the bitmap for next time                           */

	for (ii=0; ii<snode_nnz; ii++) snode_map[snode_ind[ii]] = -1;

	/* for (ii=0; ii<nrows; ii++) assert( snode_map[ii] == -1 ); */

      } else
	snode_flag = 0;

      if (snode_size >= maxsn) snode_flag = 0;

      if (j != snode_lastcol+1) snode_flag = 0;

      if (snode_flag) {

	/* THIS COLUMN BELONGS TO AN EXISTING SUPERNODE */
	
	snode_pivrows[snode_size] = pivotindex;
	snode_size++;

	for (ii=0; ii<snode_size; ii++)
	  spa[snode_pivrows[ii]] = taucs_zero;
	for (ii=0; ii<lnext; ii++)
	  spa[lu_ind[ii]] = lu_re[ii];

	for (ii=0; ii<snode_nnz; ii++)
	  snode_re[ii] = spa[snode_ind[ii]];

#ifdef DETAILED_TIMING
	time_snode_detect += (taucs_wtime() - time_tmp);
	time_tmp = taucs_wtime();
#endif

	/* Sivan: replaced 2 March 2002 */
	/* oocsp_appendcol(L,j,snode_nnz,snode_ind,snode_re);*/
	assert(Lclen[j] == 0);
	Lclen[j] = snode_nnz;
	Lappendcol(LU,j,snode_nnz,snode_ind,snode_re);

#ifdef DETAILED_TIMING
	time_append += (taucs_wtime() - time_tmp);
	bytes_appended += (double) (snode_nnz
				    * (sizeof(taucs_datatype)+sizeof(int)));
#endif

	/*taucs_printf("oocsp_numfact: supernode size = %d (column %d)\n",snode_size,j);*/
      } else {
	/* THIS COLUMN BELONGS TO A NEW SUPERNODE */

	snode_id ++;

	/*taucs_printf("oocsp_numfact: new supernode, column %d row %d\n",j,pivotindex);*/

#ifdef DETAILED_TIMING
	time_snode_detect += (taucs_wtime() - time_tmp);
	time_tmp = taucs_wtime();
#endif
      
	/*oocsp_appendcol(L,j,lnext,lu_ind,lu_re);*/
	assert(Lclen[j] == 0);
	Lclen[j] = lnext;
	Lappendcol(LU,j,lnext,lu_ind,lu_re);

#ifdef DETAILED_TIMING
	time_append += (taucs_wtime() - time_tmp);
	bytes_appended += (double) (lnext 
				    * (sizeof(taucs_datatype)+sizeof(int)));
	time_tmp = taucs_wtime();
#endif

	snode_hash = pivotindex;
	for (ii=0; ii<lnext; ii++) {
	  snode_hash ^= lu_ind[ii];
	  snode_ind[ii] = lu_ind[ii];
	}
	snode_pivrows[0] = pivotindex;
	snode_nnz = lnext;
	snode_size = 1;

#ifdef DETAILED_TIMING
	time_snode_detect += (taucs_wtime() - time_tmp);
#endif
      }

      snode_index[j] = snode_id;

      snode_lastcol = j;

#else /* SNODES */

#ifdef DETAILED_TIMING
      time_tmp = taucs_wtime();
#endif
      /*oocsp_appendcol(L,j,lnext,lu_ind,lu_re);*/
      assert(Lclen[j] == 0);
      Lclen[j] = lnext;
      Lappendcol(LU,j,lnext,lu_ind,lu_re);
#ifdef DETAILED_TIMING
      time_append += (taucs_wtime() - time_tmp);
      bytes_appended += (double) (lnext 
				  * (sizeof(taucs_datatype)+sizeof(int)));
#endif /* DETAILED_TIMING */

#endif /* SNODES */
     
      /* Update rest of panel */

      /*taucs_printf("oocsp_numfact: updating\n");*/

#ifdef DETAILED_TIMING
      time_tmp = taucs_wtime();
#endif
      if (spawidth <= 0) {
	for(qp = rowlists_head[pivotindex]; qp != -1; qp = rowlists_next[qp]) {
	  q = rowlists_colind[qp];
	  scatter(panel_nnz[q],panel_re[q],panel_ind[q],
		  spa,spamap);
	  spcol_spcol_update(pivotindex,
			     lu_re,lu_ind,lnext,
			     q,panel_inrowlist[q],
			     spa,
			     spamap,
			     panel_ind[q],&(panel_nnz[q]),
			     lcolcount[panel_id[q]]+ucolcount[panel_id[q]]);
	  gather(panel_nnz[q],panel_re[q],panel_ind[q],
		 spa,spamap);
	}
      } else {
	int spa_n = 0;
	for(qp = rowlists_head[pivotindex]; qp != -1; qp = rowlists_next[qp]) {
	  q = rowlists_colind[qp];
	  spa_updcols[spa_n] = q;
	  spa_n++;
	}
	/*
	for(qp = 0; qp < spa_n; qp++) {
	  q = spa_updcols[qp];
	  spcol_spcol_update(pivotindex,
			     lu_re,lu_ind,lnext,
			     q,
			     panel_inrowlist[q],
			     panel_spa    + (q*nrows),
			     panel_spamap + (q*nrows),
			     panel_ind[q],&(panel_nnz[q]),
			     0);
	}
	*/

	spcol_spa_update(pivotindex,
			 lu_re,lu_ind,lnext,
			 spa_updcols,spa_n,nrows,
			 panel_inrowlist,
			 panel_spa,
			 panel_spamap,
			 panel_ind,panel_nnz);
      }
#ifdef DETAILED_TIMING
      time_colcol   += (taucs_wtime() - time_tmp);
      time_colcol_2 += (taucs_wtime() - time_tmp);
#endif

      /*taucs_printf("oocsp_numfact: done updating\n");*/

      /* Eject from panel */

      for (ip = 0; ip<panel_nnz[p]; ip++) {
	i = panel_ind[p][ip];
	rowlists_delete(i,panel_inrowlist[p][ip]);
      }

      taucs_free(panel_inrowlist[p]); panel_inrowlist[p] = NULL;
      taucs_free(panel_ind[p]); panel_ind[p] = NULL;
      taucs_free(panel_re[p]);  panel_re[p]  = NULL;
      if (spawidth > 0) {
#ifndef SPA_ONEARRAY
      taucs_free(panel_spa[p]); panel_spa[p] = NULL;
      taucs_free(panel_spamap[p]);  panel_spamap[p]  = NULL;
#endif
      }
      panel_id[p] = -1;

      /*taucs_printf("oocsp_numfact: done with col\n");*/

      en++;
    }    
  }

#ifdef SNODESHIST
#define HIST_SIZE 32
#define HIST_INC  4
  {
    int i,j,k,s,last;
    int hist[HIST_SIZE];
    for (i=0; i<HIST_SIZE; i++) hist[i] = 0;

    last = 0;
    s = 0;
    for (j=0; j<ncols; j++) {
      if (snode_index[j] == last) {
	s++;
      } else {
	i = s / HIST_INC;
	if (i > HIST_SIZE-1) i=HIST_SIZE-1;
	(hist[i])++;
	/*printf("last snode index %d size %d hist %d\n",last,s,i);*/
	last = snode_index[j];
	s = 1;
      }
    }
    /*printf("last snode index %d size %d\n",last,s);*/
    i = s / HIST_INC;
    if (i > HIST_SIZE-1) i=HIST_SIZE-1;
    (hist[i])++;
    taucs_printf("oocsp_numfact: snode histogram:\n");
    for (i=0; i<HIST_SIZE-1; i++)
      taucs_printf("oocsp_numfact:   %02d-%02d: %d snodes\n",i*HIST_INC+1,(i+1)*HIST_INC,hist[i]);
    taucs_printf("oocsp_numfact:   %02d&up: %d snodes:\n",
	       (HIST_SIZE-1)*HIST_INC+1,hist[HIST_SIZE-1]);
  }
#endif

  taucs_io_append(LU, HEADER_NROWS  , 1,      1, TAUCS_INT, &(A->m) );
  taucs_io_append(LU, HEADER_NCOLS  , 1,      1, TAUCS_INT, &(A->n) );
  taucs_io_append(LU, HEADER_FLAGS  , 1,      1, TAUCS_INT, &(A->flags));
  taucs_io_append(LU, HEADER_COLPERM, (A->n), 1, TAUCS_INT, colperm    );
  taucs_io_append(LU, HEADER_IPIVOTS, (A->n), 1, TAUCS_INT, ipivots    );
  taucs_io_append(LU, HEADER_LCLEN  , (A->n), 1, TAUCS_INT, Lclen      );
  taucs_io_append(LU, HEADER_UCLEN  , (A->n), 1, TAUCS_INT, Uclen      );

  (remaining_memory) += (double) ( 2 * (ncols+1) * sizeof(int));

  {
    double nnzL = 0;
    double nnzU = 0;
    for (j=0; j < (A->n); j++) nnzL += (double) (Lclen[j]);
    for (j=0; j < (A->n); j++) nnzU += (double) (Uclen[j]);
    taucs_printf("oocsp_numfact: nnz(L) = %.2le nnz(U) = %.2le nnz(L+U) = %.2leM\n",
	       nnzL,nnzU,(nnzL+nnzU)/1e6);
  }
  
  if (spawidth > 0) {
    taucs_free(panel_spa);
    taucs_free(panel_spamap);
  }

  spa_finalize();

  rowlists_finalize();

  taucs_free(ipivots);

  taucs_free(Uclen);
  taucs_free(Lclen);

  taucs_free(panel_re);
  taucs_free(panel_ind);
  taucs_free(panel_inrowlist);
  taucs_free(panel_nnz);
  taucs_free(panel_id);

/*
  taucs_free(update_vec);
  taucs_free(update_tmp);
*/

#ifdef SNODES
  taucs_free(snode_index);
  taucs_free(snode_pivrows);
  taucs_free(snode_ind);
  taucs_free(snode_re);
  taucs_free(snode_map);

  taucs_free(pivots);

  taucs_free(spa_updcols);
  taucs_free(P);
  taucs_free(S);
  taucs_free(srows);
#endif

  taucs_free(heap);
  taucs_free(nnzmap);
  taucs_free(lindices);
  /*taucs_free(uindices);*/
  taucs_free(lu_re);
  taucs_free(lu_ind);

  /*taucs_free(pivots);
  taucs_free(ipivots);
  */

  /* OLD */

  time_total = taucs_wtime() - time_total;
  taucs_printf("oocsp_numfact: %lg sec total\n",time_total);

#ifdef DETAILED_TIMING
  taucs_printf("oocsp_numfact: %lg extra flops, %2.0lf %\n",flops_extra,100.0*flops_extra/flops);
  taucs_printf("oocsp_numfact: %lg dense flops, %2.0lf %\n",flops_dense,100.0*flops_dense/flops);
  taucs_printf("oocsp_numfact: %lg Mflop dense/s\n",(flops_dense*1e-6)/(time_snode_dense));
 
  taucs_printf("oocsp_numfact: %lg flops\n",flops);
  taucs_printf("oocsp_numfact: %lg scatter ops\n",scatters);
  taucs_printf("oocsp_numfact: %lg gather  ops\n",gathers);
  taucs_printf("oocsp_numfact: %lg heap    ops\n",num_heap_ops);
  taucs_printf("oocsp_numfact: %lg rowlist ops\n",rowlist_ops);
  taucs_printf("oocsp_numfact: %lg sec col/col ops\n",time_colcol);
  taucs_printf("oocsp_numfact: %lg sec col/col ops (%lg+%lg=%lg)\n",time_colcol,
	     time_colcol_1,time_colcol_2,(time_colcol_1+time_colcol_2));
  taucs_printf("oocsp_numfact: %lg sec column factor ops\n",time_factor);
  taucs_printf("oocsp_numfact: %lg sec scatter ops\n",time_scatter);
  taucs_printf("oocsp_numfact: %lg sec gather  ops\n",time_gather);
  taucs_printf("oocsp_numfact: %lg sec io read\n",time_read);
  taucs_printf("oocsp_numfact: %lg sec io write\n",time_append);
#ifdef SNODES
  taucs_printf("oocsp_numfact: %lg sec snode preparations for dense ops\n",time_snode_prepare);
  taucs_printf("oocsp_numfact: %lg sec snode dense ops\n",time_snode_dense);
  taucs_printf("oocsp_numfact: %lg sec snode detection\n",time_snode_detect);
  taucs_printf("oocsp_numfact: %lg sec snode 1\n",time_snode_1);
  taucs_printf("oocsp_numfact: %lg sec snode 2 (%lg)\n",time_snode_2,
		                                      time_snode_21);
  taucs_printf("oocsp_numfact: %lg sec snode 3\n",time_snode_3);
  taucs_printf("oocsp_numfact: %lg sec snode 4\n",time_snode_4);
#endif
  taucs_printf("oocsp_numfact: \n");
  taucs_printf("oocsp_numfact: %lg Mflop/s\n",(flops*1e-6)/(time_total));
  taucs_printf("oocsp_numfact: %lg MB/s IO read\n",(bytes_read*1e-6)/(time_read));
  taucs_printf("oocsp_numfact: %lg MB/s IO write\n",(bytes_appended*1e-6)/(time_append));
  taucs_printf("oocsp_numfact: %lg col reuse\n",col_ooc_updates/col_read);
  taucs_printf("oocsp_numfact: %lg percent IO\n",(time_read+time_append)/time_total);
  taucs_printf("oocsp_numfact: %lg percent col/col\n",(time_colcol)/time_total);
#endif


#ifndef SIMPLE_COL_COL
  {
    extern int oocsp_spcol_n1,oocsp_spcol_n2;
    taucs_printf("oocsp_numfact: spcol counts %d %d\n",oocsp_spcol_n1,oocsp_spcol_n2);
  }
#endif
}


/*********************************************************/
/* FACTOR                                                */
/*********************************************************/

static 
int
oocsp_factor(taucs_ccs_matrix* A_in,
	     taucs_io_handle* LU,
             int*  colperm)
{
  /*int i,j; omer*/
  
  int* etree=NULL;
  int* postorder=NULL;
  int* l_colcounts=NULL;
  int* u_colcounts=NULL;
  int* panels=NULL;
  int* schedstart=NULL;
  int* schedend=NULL;
  int* fetchnext=NULL;
  int* ejectnext=NULL;
  int  spawidth;
  int  maxsn;

  taucs_printf("taucs_ooc_lu: starting\n");

  taucs_printf("taucs_ooc_lu: calling colanalyze\n");
  oocsp_colanalyze(A_in,
		taucs_io_get_basename(LU),
		colperm,&etree,&postorder,&l_colcounts,&u_colcounts);
  
  taucs_printf("taucs_ooc_lu: calling panelize\n");
  oocsp_panelize_simple(A_in->m,A_in->n,
		        postorder,
			l_colcounts,u_colcounts,etree,
			&spawidth,&maxsn,
			&panels,&schedstart,&schedend,&fetchnext,&ejectnext);

  taucs_printf("taucs_ooc_lu: calling numfact\n");
  oocsp_numfact(A_in,colperm,
		/*L,U,*/
		LU,
		panels,
		schedstart,
		schedend,
		fetchnext,
		ejectnext,
		l_colcounts,
		u_colcounts,
		spawidth,maxsn);

  taucs_printf("taucs_ooc_lu: done, returning\n");

  return 0;
}

void
taucs_dtl(ooc_factor_lu)(taucs_ccs_matrix* A_in,
		         int    colperm[],
                         taucs_io_handle* LU,
		         double memory)
{
  remaining_memory = memory;
  taucs_printf("taucs_ooc_factor_lu: using %.0lf MBytes of in-core memory\n",
	     (remaining_memory)/1048576.0);
  oocsp_factor(A_in,LU,colperm);
}

/*********************************************************/
/* SOLVE                                                 */
/*********************************************************/

static int
oocsp_solve(taucs_io_handle* LU,
	    taucs_datatype* X,
	    taucs_datatype* B)
{
  int i,ip,j,n;
  taucs_datatype* Y;
  taucs_datatype  Aij;
  taucs_datatype* values;
  int*    indices;
  int*    irowperm;
  int*    Lclen;
  int*    Uclen;
  int*    colperm;
  int*    ipivots;

  int     found;
  
  double  time_solve = taucs_wtime();
  double  bytes_read = 0;

  taucs_printf("oocsp_solve: starting\n");

  /* n = L->nrows; */
  taucs_io_read(LU, HEADER_NROWS, 1, 1, TAUCS_INT, &n);

  Y       = (taucs_datatype*) taucs_malloc(n * sizeof(taucs_datatype));
  values  = (taucs_datatype*) taucs_malloc(n * sizeof(taucs_datatype));
  indices = (int*)    taucs_malloc(n * sizeof(int));
  irowperm= (int*)    taucs_malloc(n * sizeof(int));

  Lclen   = (int*)    taucs_malloc(n * sizeof(int));
  Uclen   = (int*)    taucs_malloc(n * sizeof(int));

  colperm = (int*)    taucs_malloc(n * sizeof(int));
  ipivots = (int*)    taucs_malloc(n * sizeof(int));
  assert(Y && values && indices && irowperm && Lclen && Uclen && colperm && ipivots);

  taucs_io_read(LU, HEADER_LCLEN, n, 1, TAUCS_INT, Lclen);
  taucs_io_read(LU, HEADER_UCLEN, n, 1, TAUCS_INT, Uclen);

  taucs_io_read(LU, HEADER_COLPERM, n, 1, TAUCS_INT, colperm);
  taucs_io_read(LU, HEADER_IPIVOTS, n, 1, TAUCS_INT, ipivots);

  for(i=0; i<n; i++)
     irowperm[ipivots[i]]=i;

  /*
  taucs_printf("colperm = [");
  for(i=0; i<n; i++)
    taucs_printf("%d, ",colperm[i]);
  taucs_printf("\b\b]\n");

  taucs_printf("ipivots = [");
  for(i=0; i<n; i++)
    taucs_printf("%d, ",ipivots[i]);
  taucs_printf("\b\b]\n");
  */
  

  /* start by permuting B, X=PB */

  /*
  for (i=0; i<n; i++)
    PB[i] = B[ ipivots[i] ];
  */

  for (i=0; i<n; i++)
    Y[i] = B[i];

  for (j=0; j<n; j++) {
    /*oocsp_readcol(L,j,indices,values);*/
    Lreadcol(LU,j,Lclen[j],indices,values);
    bytes_read += Lclen[j] * (sizeof(int) + sizeof(taucs_datatype));
    for (ip=0; ip < Lclen[j]; ip++) {
      i = indices[ip];
      Aij = values[ip];
      /*Y[i] = Y[i] - Aij*Y[irowperm[j]];*/
      Y[i] = taucs_sub(Y[i] , taucs_mul( Aij , Y[irowperm[j]] ));
    }
  }

  for (i=0; i<n; i++) X[i] = Y[i];

  for (j=n-1; j>=0; j--) {
    /*oocsp_readcol(U,j,indices,values);*/
    Ureadcol(LU,j,Uclen[j],indices,values);
    bytes_read += Uclen[j] * (sizeof(int) + sizeof(taucs_datatype));

    found = 0;
    for (ip=0; ip < Uclen[j]; ip++) {
      i = indices[ip];
      if (i == irowperm[j]) {
	found = 1;
	Aij = values[ip];
	/*X[i] = X[i] / Aij;*/
	X[i] = taucs_div( X[i] , Aij );
	values[ip] = taucs_zero; /* so we don't multiply in the daxpy below */
      }
    }
    assert( found );

    for (ip=0; ip < Uclen[j]; ip++) {
      i = indices[ip];
      Aij = values[ip];
      /*X[i] = X[i] - Aij*X[irowperm[j]];*/
      X[i] = taucs_sub( X[i] , taucs_mul( Aij , X[irowperm[j]] ));
    }
  }


  for (i=0; i<n; i++) Y[i] = X[i];
  for (i=0; i<n; i++)
    X[ ipivots[i] ] = Y[ i ];


  for (i=0; i<n; i++) Y[i] = X[i];
  for (i=0; i<n; i++)
    X[ colperm[i] ] = Y[ i ];

  /*
    X[ colperm[i] ] = Y[ i ];

  for (i=0; i<n; i++) Y[i] = X[i];
  for (i=0; i<n; i++)
    X[ ipivots[i] ] = Y[ i ];

  */

  /*
  for (i=0; i<n; i++) if (colperm[i] == 0) j=i;
  printf("rowperm[0]=%d irowperm[0]=%d colperm[0]=%d icp[0]=%d\n",ipivots[0], irowperm[0], colperm[0],j);
  for (i=0; i<n; i++) if (colperm[i] == 1) j=i;
  printf("rowperm[1]=%d irowperm[1]=%d colperm[1]=%d icp[1]=%d\n",ipivots[1], irowperm[1], colperm[1],j);
  */

  taucs_free(Y);
  taucs_free(values);
  taucs_free(indices);
  taucs_free(irowperm);
  taucs_free(Uclen);
  taucs_free(Lclen);
  taucs_free(ipivots);
  taucs_free(colperm);

  time_solve = (taucs_wtime() - time_solve);
  taucs_printf("oocsp_solve: done in %.0lf seconds, read %.0lf bytes (%.0lf MBytes)\n",
	     time_solve,bytes_read,bytes_read/1048576.0);

  return 0; /* success */
  
} /* main */

int taucs_dtl(ooc_solve_lu)(taucs_io_handle*   LU,
			    taucs_datatype* x, taucs_datatype* b)
{
  /*oocsp_solve(LU,x,b); omer*/
	return oocsp_solve(LU,x,b);
}

#endif /*#ifndef TAUCS_CORE_GENERAL*/

/*********************************************************/
/* USER CALLABLE ROUTINES                                */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL

int taucs_ooc_factor_lu(taucs_ccs_matrix* A,
			int*              colperm,
                        taucs_io_handle*  LU,                       
			double            memory)
{
#ifdef TAUCS_CONFIG_DREAL
  if (A->flags & TAUCS_DOUBLE) {
    taucs_dooc_factor_lu(A,colperm,LU,memory);
    return 0;
  }
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (A->flags & TAUCS_DCOMPLEX) {
    taucs_zooc_factor_lu(A,colperm,LU,memory);
    return 0;
  }
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (A->flags & TAUCS_SINGLE) {
    taucs_sooc_factor_lu(A,colperm,LU,memory);
    return 0;
  }
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (A->flags & TAUCS_SCOMPLEX) {
    taucs_cooc_factor_lu(A,colperm,LU,memory);
    return 0;
  }
#endif

  assert(0);
  return -1;
}

int taucs_ooc_solve_lu (taucs_io_handle* LU,
			void* x, void* b)
{
  int flags;

  taucs_io_read(LU, HEADER_FLAGS, 1, 1, TAUCS_INT, &flags);

  printf("taucs_ooc_solve_lu: starting, DZSC=%d%d%d%d\n",
	 (flags & TAUCS_DOUBLE  ) != 0,
	 (flags & TAUCS_DCOMPLEX) != 0,
	 (flags & TAUCS_SINGLE  ) != 0,
	 (flags & TAUCS_SCOMPLEX) != 0);
  
#ifdef TAUCS_CONFIG_DREAL
  if (flags & TAUCS_DOUBLE) {
    taucs_dooc_solve_lu(LU,x,b);
    return 0;
  }
#endif

#ifdef TAUCS_CONFIG_DCOMPLEX
  if (flags & TAUCS_DCOMPLEX) {
    taucs_zooc_solve_lu(LU,x,b);
    return 0;
  }
#endif

#ifdef TAUCS_CONFIG_SREAL
  if (flags & TAUCS_SINGLE) {
    taucs_sooc_solve_lu(LU,x,b);
    return 0;
  }
#endif

#ifdef TAUCS_CONFIG_SCOMPLEX
  if (flags & TAUCS_SCOMPLEX) {
    taucs_cooc_solve_lu(LU,x,b);
    return 0;
  }
#endif

  assert(0);
  return -1;
}

#endif /*#ifdef TAUCS_CORE_GENERAL*/

/*********************************************************/
/* END OF FILE                                           */
/*********************************************************/


