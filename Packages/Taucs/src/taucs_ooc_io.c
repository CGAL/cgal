/*********************************************************/
/* TAUCS                                                 */
/* Author: Vladimir Rotking and Sivan Toledo             */
/*                                                       */
/* Out-of-Core Sparse Matrix I/O Subroutines             */
/*********************************************************/

/*************************************************************/
/*                                                           */
/*************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "taucs.h"

#include <assert.h>
#include <math.h>

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

/*************************************************************/
/* io routines                                               */
/*************************************************************/

/*
#define TAUCS_PACKED  256
#define TAUCS_BYROW   512
*/
/* #define IO_TYPE_COMPLEX       2 */

/* out-of-core matrix file types */

#define IO_TYPE_SINGLEFILE    1
#define IO_TYPE_MULTIFILE     0

/* maximum file size in GB */

#define IO_FILE_RESTRICTION   1024

/* in taucs.h:
typedef struct {
  int   type;
  int   nmatrices;
  void* type_specific;
} taucs_io_handle;
*/

typedef struct {
  int   m;
  int   n;
  int   flags;
  off_t offset;
} taucs_io_matrix_singlefile;

typedef struct {
  int   f;
  off_t last_offset;
  taucs_io_matrix_singlefile* matrices;
} taucs_io_handle_singlefile;

typedef struct {
  int    m;
  int    n;
  int    flags;
  double offset;
} taucs_io_matrix_multifile;

typedef struct {
  int    f[1024];
  double last_offset;
  int    last_created_file;
  char   basename[256];
  taucs_io_matrix_multifile* matrices;
} taucs_io_handle_multifile;


#define TAUCS_FILE_SIGNATURE "taucs"

/*
double disc_read = 0.0;
double bytes_read =0.0;
double time_read = 0.0;
double disc_write = 0.0;
double bytes_write = 0.0;
double time_write = 0.0;
*/

/*************************************************************/
/*                                                           */
/*************************************************************/

static int element_size(int flags)
{
  if (flags & TAUCS_SINGLE)   return sizeof(taucs_single);
  if (flags & TAUCS_DOUBLE)   return sizeof(taucs_double);
  if (flags & TAUCS_SCOMPLEX) return sizeof(taucs_scomplex);
  if (flags & TAUCS_DCOMPLEX) return sizeof(taucs_dcomplex);
  if (flags & TAUCS_INT)      return sizeof(int);
  assert(0);
  return -1;
}

/*************************************************************/
/*                                                           */
/*************************************************************/

taucs_io_handle* taucs_io_create_multifile(char* basename)
{
  int f;
  ssize_t nbytes;
  int     nmatrices;
  double  offset;
  taucs_io_handle* h;
  mode_t mode;
  mode_t perm;
  char filename[256];

  sprintf(filename,"%s.%d",basename,0);

#ifdef OSTYPE_win32
  mode = _O_RDWR | _O_CREAT | _O_BINARY;
  perm = _S_IREAD | _S_IWRITE | _S_IEXEC;
#else
  mode = O_RDWR | O_CREAT;
  perm = 0644;
#endif

  f = open(filename,mode,perm);

  if (f == -1) {
    taucs_printf("taucs_create: Could not create metadata file %s\n",filename);
    return NULL;
  }

  nbytes = write(f,
		 TAUCS_FILE_SIGNATURE,
		 strlen(TAUCS_FILE_SIGNATURE));
 if (nbytes != strlen(TAUCS_FILE_SIGNATURE)) { 
    taucs_printf("taucs_create: Error writing metadata.\n");
    return NULL;
  }

  nmatrices = 0;
  offset = (double)(strlen(TAUCS_FILE_SIGNATURE)+sizeof(int)+sizeof(double));

  nbytes = write(f,&nmatrices,sizeof(int));
  if (nbytes != sizeof(int)) { 
    taucs_printf("taucs_create: Error writing metadata (2).\n");
    return NULL;
  }

  nbytes = write(f,&offset   ,sizeof(double));
  if (nbytes != sizeof(double)) { 
    taucs_printf("taucs_create: Error writing metadata (3).\n");
    return NULL;
  }

  h = (taucs_io_handle*) taucs_malloc(sizeof(taucs_io_handle));
  if (!h) {
    taucs_printf("taucs_create: out of memory (4)\n");
    return NULL;
  }
  h->type      = IO_TYPE_MULTIFILE;
  h->nmatrices = 0;
  h->type_specific = (taucs_io_handle_multifile*) taucs_malloc(sizeof(taucs_io_handle_multifile));
  if (!h->type_specific) {
    taucs_printf("taucs_create: out of memory (5)\n");
    taucs_free(h);
    return NULL;
  }
  ((taucs_io_handle_multifile*)h->type_specific)->f[0] = f;
  ((taucs_io_handle_multifile*)h->type_specific)->matrices = NULL;
  ((taucs_io_handle_multifile*)h->type_specific)->last_offset = offset;
  ((taucs_io_handle_multifile*)h->type_specific)->last_created_file = 0;
  strcpy(((taucs_io_handle_multifile*)h->type_specific)->basename,basename);

  h->nreads = h->nwrites = h->bytes_read =
    h->bytes_written = h->read_time =h->write_time = 0.0;

  return h;
}

taucs_io_handle* taucs_io_create_singlefile(char* filename)
{
  int f;
  ssize_t nbytes;
  int     nmatrices;
  off_t   offset;
  taucs_io_handle* h;
  mode_t mode;
  mode_t perm;

#ifdef OSTYPE_win32
  mode = _O_RDWR | _O_CREAT | _O_BINARY;
  perm = _S_IREAD | _S_IWRITE | _S_IEXEC;
#else
  mode = O_RDWR | O_CREAT;
  perm = 0644;
#endif

  f = open(filename,mode,perm);

  if (f == -1) {
    taucs_printf("taucs_create: Could not create metadata file %s\n",filename);
    return NULL;
  }

  nbytes = write(f,
		 TAUCS_FILE_SIGNATURE,
		 strlen(TAUCS_FILE_SIGNATURE));
 if (nbytes != strlen(TAUCS_FILE_SIGNATURE)) { 
    taucs_printf("taucs_create: Error writing metadata.\n");
    return NULL;
  }

  nmatrices = 0;
  offset = strlen(TAUCS_FILE_SIGNATURE) 
    + sizeof(int)
    + sizeof(off_t);

  nbytes = write(f,&nmatrices,sizeof(int  ));
  if (nbytes != sizeof(int)) { 
    taucs_printf("taucs_create: Error writing metadata (2).\n");
    return NULL;
  }
  nbytes = write(f,&offset   ,sizeof(off_t));
  if (nbytes != sizeof(off_t)) { 
    taucs_printf("taucs_create: Error writing metadata (3).\n");
    return NULL;
  }

  h = (taucs_io_handle*) taucs_malloc(sizeof(taucs_io_handle));
  if (!h) {
    taucs_printf("taucs_create: out of memory (4)\n");
    return NULL;
  }
  h->type      = IO_TYPE_SINGLEFILE;
  h->nmatrices = 0;
  h->type_specific = (taucs_io_handle_singlefile*) taucs_malloc(sizeof(taucs_io_handle_singlefile));
  if (!h->type_specific) {
    taucs_printf("taucs_create: out of memory (5)\n");
    taucs_free(h);
    return NULL;
  }
  ((taucs_io_handle_singlefile*)h->type_specific)->f = f;
  ((taucs_io_handle_singlefile*)h->type_specific)->matrices = NULL;
  ((taucs_io_handle_singlefile*)h->type_specific)->last_offset = offset;

  h->nreads = h->nwrites = h->bytes_read =
    h->bytes_written = h->read_time =h->write_time = 0.0;

  return h;
}

int taucs_io_append(taucs_io_handle* f,
		    int   index,
		    int   m,int n,
		    int   flags,
		    void* data
		    )
{
  int this_size = 0; /* warning */
  int next_size;
  int first_size = 0;
  int added_files = 0;
  int written_bytes = 0;
  off_t this_offset;
  double this_multi_offset,new_last_offset;
  ssize_t nbytes;
  int i;
  mode_t mode;
  mode_t perm;
  char filename[256];
  int file_id;
  double wtime;
 
  wtime = taucs_wtime();
 
  if (f->type == IO_TYPE_SINGLEFILE) {
    taucs_io_handle_singlefile* h = ((taucs_io_handle_singlefile*) f->type_specific);
    taucs_io_matrix_singlefile* matrices;
   
    if (index >= f->nmatrices){    
      ((taucs_io_handle_singlefile*)f->type_specific)->matrices = 
	(taucs_io_matrix_singlefile*) taucs_realloc(h->matrices,
					      (index + 1) * 
					      sizeof(taucs_io_matrix_singlefile));
      for(i=f->nmatrices;i<index;i++){
	h->matrices[i].m = -1;
	h->matrices[i].n = -1;
	h->matrices[i].flags = -1;
	h->matrices[i].offset = -1;
      }
      f->nmatrices = index+1;
    }
    else
      if(h->matrices[index].m!=-1||h->matrices[index].n!=-1){
	taucs_printf("taucs_append: try append more than once for index=%d \n",index);
	return -1;
      }
    
    if (!((taucs_io_handle_singlefile*)f->type_specific)->matrices) {
      taucs_printf("taucs_append: out of memory \n");
      return -1;
    }

    matrices = h->matrices;
    this_offset = h->last_offset;
    matrices[index].m = m;
    matrices[index].n = n;
    matrices[index].flags = flags;
    matrices[index].offset = this_offset;
    /*this_size = m * n * ((flags & TAUCS_INT) ? sizeof(int) : sizeof(double));*/
    this_size = m * n * element_size(flags);
    h->last_offset +=this_size; 
    

    /*taucs_printf("debug1: index = %d offset = %d\n ",index,this_offset);*/
    if (lseek(h->f, this_offset, SEEK_SET) == -1) {
      taucs_printf("taucs_append: lseek failed\n");
      return -1;
    }
    
    nbytes = write(h->f, data, this_size);
    /*    taucs_printf("debug for nbytes = %d this_size = %d \n ",nbytes,this_size);*/
    /*if (nbytes != this_size) { omer*/
		if ((int)nbytes != this_size) { 
      taucs_printf("taucs_append: Error writing data (%s:%d).\n",__FILE__,__LINE__);
      return -1;
    }
  }

  if (f->type == IO_TYPE_MULTIFILE) {
    taucs_io_handle_multifile* h = ((taucs_io_handle_multifile*) f->type_specific);
    taucs_io_matrix_multifile* matrices;
   
    if (index >= f->nmatrices){    
      ((taucs_io_handle_multifile*)f->type_specific)->matrices = 
	(taucs_io_matrix_multifile*) taucs_realloc(h->matrices,
					      (index + 1) * 
					      sizeof(taucs_io_matrix_multifile));
      for(i=f->nmatrices;i<index;i++){
	h->matrices[i].m = -1;
	h->matrices[i].n = -1;
	h->matrices[i].flags = -1;
	h->matrices[i].offset = -1.0;
      }
      f->nmatrices = index+1;
    }
    else
      if(h->matrices[index].m!=-1||h->matrices[index].n!=-1){
	taucs_printf("taucs_append: try append more than once for index=%d \n",index);
	return -1;
      }
    
    if (!((taucs_io_handle_multifile*)f->type_specific)->matrices) {
      taucs_printf("taucs_append: out of memory \n");
      return -1;
    }

    matrices = h->matrices;
    matrices[index].m = m;
    matrices[index].n = n;
    matrices[index].flags = flags;
    matrices[index].offset =  h->last_offset;
    /*this_size = m * n * ((flags & TAUCS_INT) ? sizeof(int) : sizeof(double));*/
    this_size = m * n * element_size(flags);
    new_last_offset = h->last_offset + ((double)this_size);
    /*    taucs_printf("debug1: index = %d offset = %lf\n ",index,h->last_offset);*/
    
    if(new_last_offset < ((h->last_created_file+1)*IO_FILE_RESTRICTION*1024.0*1024.0)){    
      this_multi_offset = h->last_offset - ((h->last_created_file)*IO_FILE_RESTRICTION*1024.0*1024.0);   
      if (lseek(h->f[h->last_created_file],(off_t)this_multi_offset, SEEK_SET) == -1) {
	taucs_printf("taucs_append: lseek failed\n");
	return -1;
      }
					    
      nbytes = write(h->f[h->last_created_file], data, this_size);
      /*if (nbytes != this_size) { omer*/
			if ((int)nbytes != this_size) { 
	taucs_printf("taucs_append: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	taucs_printf("taucs_append: index %d n %d m %d\n",index,n,m);
	taucs_printf("taucs_append: trying to write %d bytes from %08x, wrote %d\n",
		     this_size,data,nbytes);
	if (nbytes==-1) perror("taucs_append");
	return -1;
      }
    }
    else{
      if(h->last_offset < ((h->last_created_file+1)*IO_FILE_RESTRICTION*1024.0*1024.0)){
	this_multi_offset = h->last_offset - 
	  ((h->last_created_file)*IO_FILE_RESTRICTION*1024.0*1024.0);
	if (lseek(h->f[h->last_created_file],(off_t)this_multi_offset, SEEK_SET) == -1) {
	  taucs_printf("taucs_append: lseek failed\n");
	  return -1;
	}
	first_size = (int)((IO_FILE_RESTRICTION*1024.0*1024.0) - this_multi_offset);
	nbytes = write(h->f[h->last_created_file], data, first_size);
	/*if (nbytes != first_size) { omer*/
	if ((int)nbytes != first_size) { 
	  taucs_printf("taucs_append: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	  return -1;
	}
      }

      this_multi_offset = 0.0;
      next_size = this_size - first_size;
      written_bytes = first_size;
      while(next_size>0){
	if(next_size>IO_FILE_RESTRICTION*1024*1024)
	  next_size = IO_FILE_RESTRICTION*1024*1024;
	sprintf(filename,"%s.%d",h->basename,(h->last_created_file+1));

#ifdef OSTYPE_win32
	mode = _O_RDWR | _O_CREAT | _O_BINARY;
	perm = _S_IREAD | _S_IWRITE | _S_IEXEC;
#else
	mode = O_RDWR | O_CREAT;
	perm = 0644;
#endif

	file_id = open(filename,mode,perm);
      
	if (file_id == -1) {
	  taucs_printf("taucs_append: Could not create metadata file %s\n",filename);
	  return -1;
	}
	added_files++;
	h->last_created_file++;
	h->f[h->last_created_file] = file_id;
	nbytes = write(h->f[h->last_created_file],((char*)data)+written_bytes,next_size);
	/*if (nbytes != next_size) { omer*/
	if ((int)nbytes != next_size) { 
	  taucs_printf("taucs_append: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	  return -1;
	}
	written_bytes += next_size;
	next_size = this_size - written_bytes;
      }
    }
    h->last_offset = new_last_offset; 
 }
  
  wtime = taucs_wtime()-wtime;

  f->nwrites       += 1.0;
  f->bytes_written += (double) this_size;
  f->write_time    += wtime;

  /*disc_write += 1.0;*/
  /*bytes_write += (double)this_size;*/
  /*time_write += wtime;*/

  return 0;
}

int   taucs_io_write(taucs_io_handle* f,
		     int   index,
		     int   m,int n,
		     int   flags,
		     void* data
		     )
{
  int this_size;
  off_t this_offset;
  ssize_t nbytes;
  double curr_file_offset;
  int first_size;
  int next_size,start_file_index;
  int write_bytes;

  if (f->type == IO_TYPE_SINGLEFILE) {
    taucs_io_handle_singlefile* h = ((taucs_io_handle_singlefile*) f->type_specific);
    taucs_io_matrix_singlefile* matrices;
    
    if (index>=f->nmatrices) return -1;
    matrices = h->matrices;
    /*this_size = m * n * ((flags & TAUCS_INT) ? sizeof(int) : sizeof(double));*/
    this_size = m * n * element_size(flags);
    this_offset = matrices[index].offset;
        
    if (lseek(h->f, this_offset, SEEK_SET) == -1) {
      taucs_printf("taucs_write: lseek failed\n");
      return -1;
    }
    nbytes = write(h->f, data, this_size);
    /*if (nbytes != this_size) { omer*/
		if ((int)nbytes != this_size) { 
      taucs_printf("taucs_write: Error writing data (%s:%d).\n",__FILE__,__LINE__);
      return -1;
    }
  }

if (f->type == IO_TYPE_MULTIFILE) {
    taucs_io_handle_multifile* h = ((taucs_io_handle_multifile*) f->type_specific);
    taucs_io_matrix_multifile* matrices;
    
    if (index>=f->nmatrices) return -1;
    matrices = h->matrices;
    /*this_size = m * n * ((flags & TAUCS_INT) ? sizeof(int) : sizeof(double));*/
    this_size = m * n * element_size(flags);
    start_file_index = (int)floor((matrices[index].offset/(IO_FILE_RESTRICTION*1024*1024)));
    curr_file_offset = matrices[index].offset - start_file_index*(IO_FILE_RESTRICTION*1024.0*1024.0);

    if (lseek(h->f[start_file_index],(off_t) curr_file_offset, SEEK_SET) == -1) {
      taucs_printf("taucs_write: lseek failed\n");
      return -1;
    }
    first_size = (int)((IO_FILE_RESTRICTION*1024.0*1024.0) - curr_file_offset);
    if(this_size<first_size) first_size = this_size;

    nbytes = write(h->f[start_file_index], data, first_size);
    /*if (nbytes != first_size) { omer*/
		if ((int)nbytes != first_size) { 
      taucs_printf("taucs_write: Error writing data (%s:%d).\n",__FILE__,__LINE__);
      return -1;
    }
    next_size = this_size - first_size;
    write_bytes = first_size;

    while(next_size>0){
      if(next_size>IO_FILE_RESTRICTION*1024*1024)
	next_size = IO_FILE_RESTRICTION*1024*1024;
	start_file_index++;
	if (lseek(h->f[start_file_index],0, SEEK_SET) == -1) {
	  taucs_printf("taucs_write: lseek failed\n");
	  return -1;
	}
	nbytes = write(h->f[start_file_index],((char*)data)+write_bytes,next_size);
	/*if (nbytes != next_size) { omer*/
	if ((int)nbytes != next_size) { 
	  taucs_printf("taucs_write: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	  return -1;
	}
	write_bytes += next_size;
	next_size = this_size - write_bytes;
      }
  }

  return 0;
}

int   taucs_io_read(taucs_io_handle* f,
		    int    index,
		    int    m,int n,
		    int    flags,
		    void* data
		    )
{
  int this_size = 0; /* warning */
  off_t this_offset;
  ssize_t nbytes;
  double curr_file_offset;
  int first_size;
  int next_size,start_file_index;
  int read_bytes;
  double wtime;

  wtime = taucs_wtime();

  if (f->type == IO_TYPE_SINGLEFILE) {
    taucs_io_handle_singlefile* h = ((taucs_io_handle_singlefile*) f->type_specific);
    taucs_io_matrix_singlefile* matrices;
    
    if (index>=f->nmatrices) return -1;
    matrices = h->matrices;
    /*this_size = m * n * ((flags & TAUCS_INT) ? sizeof(int) : sizeof(double));*/
    this_size = m * n * element_size(flags);
    this_offset = matrices[index].offset;
        
    if (lseek(h->f, this_offset, SEEK_SET) == -1) {
      taucs_printf("taucs_read: lseek failed\n");
      return -1;
    }
    nbytes = read(h->f, data, this_size);
    /*if (nbytes != this_size) { omer*/
		if ((int)nbytes != this_size) { 
      taucs_printf("taucs_read: Error reading data .\n");
      return -1;
    }
  }

  if (f->type == IO_TYPE_MULTIFILE) {
    taucs_io_handle_multifile* h = ((taucs_io_handle_multifile*) f->type_specific);
    taucs_io_matrix_multifile* matrices;
    
    if (index>=f->nmatrices) return -1;
    matrices = h->matrices;
    /*this_size = m * n * ((flags & TAUCS_INT) ? sizeof(int) : sizeof(double));*/
    this_size = m * n * element_size(flags);
    start_file_index = (int)floor((matrices[index].offset/(IO_FILE_RESTRICTION*1024*1024)));
    curr_file_offset = matrices[index].offset - start_file_index*(IO_FILE_RESTRICTION*1024.0*1024.0);

    /* for find overflow */
    assert(curr_file_offset < IO_FILE_RESTRICTION*1024.0*1024.0);

    if (lseek(h->f[start_file_index],(off_t) curr_file_offset, SEEK_SET) == -1) {
      taucs_printf("taucs_read: lseek failed\n");
      return -1;
    }
    first_size = (int)(IO_FILE_RESTRICTION*1024.0*1024.0 - curr_file_offset);
    if(this_size<first_size) first_size = this_size;

    nbytes = read(h->f[start_file_index], data, first_size);
    /*if (nbytes != first_size) { omer*/
		if ((int)nbytes != first_size) { 
      taucs_printf("taucs_read: Error reading data .\n");
      return -1;
    }
    next_size = this_size - first_size;
    read_bytes = first_size;

    while(next_size>0){
      if(next_size>IO_FILE_RESTRICTION*1024*1024)
	next_size = IO_FILE_RESTRICTION*1024*1024;
	start_file_index++;
	if (lseek(h->f[start_file_index],0, SEEK_SET) == -1) {
	  taucs_printf("taucs_read: lseek failed\n");
	  return -1;
	}
	nbytes = read(h->f[start_file_index],((char*)data)+read_bytes,next_size);
	/*if (nbytes != next_size) { omer*/
	if ((int)nbytes != next_size) { 
	  taucs_printf("taucs_read: Error reading data .\n");
	  return -1;
	}
	read_bytes += next_size;
	next_size = this_size - read_bytes;
      }
  }
  
  wtime = taucs_wtime()-wtime;

  f->nreads     += 1.0;
  f->read_time  += wtime;
  f->bytes_read += (double) this_size;

  /*disc_read += 1.0;*/
  /*bytes_read += (double)this_size;*/
  /*time_read += wtime;*/
  
  return 0;
}

int   taucs_io_close(taucs_io_handle* f)
{
  int i;
  /*int this_size,last_size; omer*/
  /*off_t offset; omer*/
  ssize_t nbytes;
  double curr_file_offset;
  mode_t mode;
  mode_t perm;
  char filename[256];
  int file_id;
  int first_size;

  if (f->type == IO_TYPE_SINGLEFILE) {
    taucs_io_handle_singlefile* h = ((taucs_io_handle_singlefile*) f->type_specific);
    taucs_io_matrix_singlefile* matrices;
    matrices = h->matrices;

    /* this no longer works! must deal with flags correctly for all data types */
    /*last_size = matrices[f->nmatrices-1].m * matrices[f->nmatrices-1].n * ((matrices[f->nmatrices-1].flags & TAUCS_INT) ? sizeof(int) : sizeof(double));
      offset = matrices[f->nmatrices-1].offset + last_size;*/

    if (lseek(h->f, strlen(TAUCS_FILE_SIGNATURE), SEEK_SET) == -1) {
      taucs_printf("taucs_close: lseek failed\n");
      return -1;
    }
    /* writing number of matrices */
    nbytes = write(h->f,&f->nmatrices, sizeof(int));
    if (nbytes != sizeof(int)) { 
      taucs_printf("taucs_close: Error writing metadata.\n");
      return -1;
    }
    /* writing start offset of metadata */
    nbytes = write(h->f,&(h->last_offset), sizeof(int));
    if (nbytes != sizeof(int)) { 
      taucs_printf("taucs_close: Error writing metadata.\n");
      return -1;
    }
    if (lseek(h->f, h->last_offset, SEEK_SET) == -1) {
	taucs_printf("taucs_close: lseek failed\n");
	return -1;
      }
    /* writing metadata  for every matrix */
    for(i=0; i<f->nmatrices; i++){
      nbytes = write(h->f,&matrices[i].m, sizeof(int));
      if (nbytes != sizeof(int)) { 
	taucs_printf("taucs_close: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	return -1;
	}
      nbytes = write(h->f,&matrices[i].n, sizeof(int));
      if (nbytes != sizeof(int)) { 
	taucs_printf("taucs_close: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	return -1;
	}
      nbytes = write(h->f,&matrices[i].flags, sizeof(int));
      if (nbytes != sizeof(int)) { 
	taucs_printf("taucs_close: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	return -1;
	}
      nbytes = write(h->f,&matrices[i].offset, sizeof(off_t));
      if (nbytes != sizeof(off_t)) { 
	taucs_printf("taucs_close: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	return -1;
      }
    }
    taucs_free(matrices);
  }

  if (f->type == IO_TYPE_MULTIFILE) {
    taucs_io_handle_multifile* h = ((taucs_io_handle_multifile*) f->type_specific);
    taucs_io_matrix_multifile* matrices;
    
    matrices = h->matrices;

    if (lseek(h->f[0], strlen(TAUCS_FILE_SIGNATURE), SEEK_SET) == -1) {
      taucs_printf("taucs_close: lseek failed\n");
      return -1;
    }
    /* writing number of matrices */
    nbytes = write(h->f[0],&(f->nmatrices), sizeof(int));
    if (nbytes != sizeof(int)) { 
      taucs_printf("taucs_close: Error writing metadata.\n");
      return -1;
    }
    /* writing start offset of metadata */
    nbytes = write(h->f[0],&(h->last_offset), sizeof(double));
    if (nbytes != sizeof(double)) { 
	taucs_printf("taucs_close: Error writing metadata.\n");
	return -1;
    }

    curr_file_offset = h->last_offset - (h->last_created_file)*(IO_FILE_RESTRICTION*1024.0*1024.0);
    if(!((int)curr_file_offset)){
      sprintf(filename,"%s.%d",h->basename,(h->last_created_file+1));
#ifdef OSTYPE_win32
      mode = _O_RDWR | _O_CREAT | _O_BINARY;
      perm = _S_IREAD | _S_IWRITE | _S_IEXEC;
#else
      mode = O_RDWR | O_CREAT;
      perm = 0644;
#endif

      file_id = open(filename,mode,perm);
      
      if (file_id == -1) {
	taucs_printf("taucs_close: Could not create metadata file %s\n",filename);
	return -1;
      }
      h->last_created_file++;
	h->f[h->last_created_file] = file_id;
    }
    else
      if (lseek(h->f[h->last_created_file],(off_t) curr_file_offset, SEEK_SET) == -1) {
	taucs_printf("taucs_close: lseek failed\n");
	return -1;
      }
    
    /* writing metadata  for every matrix */
    for(i=0; i<f->nmatrices; i++){
      /* write m */
      if(curr_file_offset+(double)sizeof(int)<IO_FILE_RESTRICTION*1024.0*1024.0){
	nbytes = write(h->f[h->last_created_file],&(matrices[i].m),sizeof(int));
	if (nbytes != sizeof(int)){ 
	  taucs_printf("taucs_close: Error writing data (%s:%d).\n",__FILE__,__LINE__);
	  return -1;
	}
	curr_file_offset += (double)sizeof(int);
      }
      else{
	first_size = IO_FILE_RESTRICTION*1024*1024 - (int)curr_file_offset;
	nbytes = write(h->f[h->last_created_file],&(matrices[i].m),first_size);
	/*if (nbytes != first_size) { omer*/
	if ((int)nbytes != first_size) { 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	sprintf(filename,"%s.%d",h->basename,(h->last_created_file+1));

#ifdef OSTYPE_win32
	mode = _O_RDWR | _O_CREAT | _O_BINARY;
	perm = _S_IREAD | _S_IWRITE | _S_IEXEC;
#else
	mode = O_RDWR | O_CREAT;
	perm = 0644;
#endif

	file_id = open(filename,mode,perm);
	
	if (file_id == -1) {
	  taucs_printf("taucs_close: Could not create metadata file %s\n",filename);
	  return -1;
	}
	h->last_created_file++;
	h->f[h->last_created_file] = file_id;
	nbytes = write(h->f[h->last_created_file],(char*)&(matrices[i].m)+first_size,sizeof(int)-first_size);
	if (nbytes != sizeof(int)-first_size){ 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	curr_file_offset = (double)(sizeof(int)-first_size);
      }
      /* write n */
      if(curr_file_offset+(double)sizeof(int)<IO_FILE_RESTRICTION*1024.0*1024.0){
	nbytes = write(h->f[h->last_created_file],&(matrices[i].n),sizeof(int));
	if (nbytes != sizeof(int)){ 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	curr_file_offset += (double)sizeof(int);
      }
      else{
	first_size = IO_FILE_RESTRICTION*1024*1024 - (int)curr_file_offset;
	nbytes = write(h->f[h->last_created_file],&(matrices[i].n),first_size);
	/*if (nbytes != first_size) { omer*/
	if ((int)nbytes != first_size) { 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	sprintf(filename,"%s.%d",h->basename,(h->last_created_file+1));
#ifdef OSTYPE_win32
	mode = _O_RDWR | _O_CREAT | _O_BINARY;
	perm = _S_IREAD | _S_IWRITE | _S_IEXEC;
#else
	mode = O_RDWR | O_CREAT;
	perm = 0644;
#endif
	file_id = open(filename,mode,perm);
	
	if (file_id == -1) {
	  taucs_printf("taucs_close: Could not create metadata file %s\n",filename);
	  return -1;
	}
	h->last_created_file++;
	h->f[h->last_created_file] = file_id;
	nbytes = write(h->f[h->last_created_file],(char*)&(matrices[i].n)+first_size,sizeof(int)-first_size);
	if (nbytes != sizeof(int)-first_size){ 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	curr_file_offset = (double)(sizeof(int)-first_size);
      }
      
      /* write flags */
      if(curr_file_offset+(double)sizeof(int)<IO_FILE_RESTRICTION*1024.0*1024.0){
	nbytes = write(h->f[h->last_created_file],&(matrices[i].flags),sizeof(int));
	if (nbytes != sizeof(int)){ 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	curr_file_offset += (double)sizeof(int);
      }
      else{
	first_size = IO_FILE_RESTRICTION*1024*1024 - (int)curr_file_offset;
	nbytes = write(h->f[h->last_created_file],&(matrices[i].flags),first_size);
	/*if (nbytes != first_size) { omer*/
	if ((int)nbytes != first_size) { 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	sprintf(filename,"%s.%d",h->basename,(h->last_created_file+1));
#ifdef OSTYPE_win32
	mode = _O_RDWR | _O_CREAT | _O_BINARY;
	perm = _S_IREAD | _S_IWRITE | _S_IEXEC;
#else
	mode = O_RDWR | O_CREAT;
	perm = 0644;
#endif
	file_id = open(filename,mode,perm);
	
	if (file_id == -1) {
	  taucs_printf("taucs_close: Could not create metadata file %s\n",filename);
	  return -1;
	}
	h->last_created_file++;
	h->f[h->last_created_file] = file_id;
	nbytes = write(h->f[h->last_created_file],(char*)&(matrices[i].flags)+first_size,sizeof(int)-first_size);
	if (nbytes != sizeof(int)-first_size){ 
	    taucs_printf("taucs_close: Error writing data .\n");
	    return -1;
	}
	curr_file_offset = (double)(sizeof(int)-first_size);
      }
      /* write offset */
      if(curr_file_offset+(double)sizeof(double)<IO_FILE_RESTRICTION*1024.0*1024.0){
	nbytes = write(h->f[h->last_created_file],&(matrices[i].offset),sizeof(double));
	if (nbytes != sizeof(double)){ 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	curr_file_offset += (double)sizeof(double);
      }
      else{
	first_size = IO_FILE_RESTRICTION*1024*1024 - (int)curr_file_offset;
	nbytes = write(h->f[h->last_created_file],&(matrices[i].offset),first_size);
	/*if (nbytes != first_size) { omer*/
	if ((int)nbytes != first_size) { 
	  taucs_printf("taucs_close: Error writing data .\n");
	  return -1;
	}
	sprintf(filename,"%s.%d",h->basename,(h->last_created_file+1));
#ifdef OSTYPE_win32
	mode = _O_RDWR | _O_CREAT | _O_BINARY;
	perm = _S_IREAD | _S_IWRITE | _S_IEXEC;
#else
	mode = O_RDWR | O_CREAT;
	perm = 0644;
#endif

	file_id = open(filename,mode,perm);
	if (file_id == -1) {
	  taucs_printf("taucs_close: Could not create metadata file %s\n",filename);
	    return -1;
	}
	h->last_created_file++;
	h->f[h->last_created_file] = file_id;
	nbytes = write(h->f[h->last_created_file],(char*)&(matrices[i].offset)+first_size,sizeof(double)-first_size);
	  if (nbytes != sizeof(double)-first_size){ 
	    taucs_printf("taucs_close: Error writing data .\n");
	    return -1;
	  }
	  curr_file_offset = (double)(sizeof(double)-first_size);
      }
    }
    for(i=0;i<=h->last_created_file;i++){
      file_id=close(h->f[i]);
      if (file_id == -1) {
	sprintf(filename,"%s.%d",h->basename,i);
	taucs_printf("taucs_close: Could not close data file %s\n",filename);
	return -1;
      }
    }
    taucs_free(matrices);
  }
  
  taucs_free(f->type_specific);
  taucs_free(f);
  
  return 0;
}

taucs_io_handle* taucs_io_open_singlefile(char* filename)
{
  int f;
  ssize_t nbytes;
  /*int     nmatrices; omer*/
  taucs_io_handle* h;
  mode_t mode;
  taucs_io_handle_singlefile* hs;
  int i;

#ifdef OSTYPE_win32
  mode = _O_RDWR | _O_BINARY;
#else
  mode = O_RDWR;
#endif
  f  = open(filename,mode);
  if (f == -1) {
    taucs_printf("taucs_open: Could not open existed data file %s\n",filename);
    return NULL;
  }

  h = (taucs_io_handle*) taucs_malloc(sizeof(taucs_io_handle));
  if (!h) {
    taucs_printf("taucs_open: out of memory (4)\n");
    return NULL;
  }
  h->type      = IO_TYPE_SINGLEFILE;
  h->type_specific = (taucs_io_handle_singlefile*) taucs_malloc(sizeof(taucs_io_handle_singlefile));
  if (!h->type_specific) {
    taucs_printf("taucs_open: out of memory \n");
    taucs_free(h);
    return NULL;
  }
  hs = h->type_specific;
  hs->f = f;
 
  if (lseek(hs->f, strlen(TAUCS_FILE_SIGNATURE), SEEK_SET) == -1) {
    taucs_printf("taucs_open: lseek failed\n");
    return NULL;
  }
  nbytes = read(hs->f, &h->nmatrices, sizeof(int));
  if (nbytes != sizeof(int)) { 
    taucs_printf("taucs_open: Error read data .\n");
    return NULL;
  }
  nbytes = read(hs->f, &hs->last_offset, sizeof(int));
  if (nbytes != sizeof(int)) { 
    taucs_printf("taucs_open: Error read data .\n");
    return NULL;
  }

  hs->matrices = 
      (taucs_io_matrix_singlefile*) taucs_malloc((h->nmatrices)* sizeof(taucs_io_matrix_singlefile));

  /* seek of start offset of data */
  if (lseek(hs->f, hs->last_offset, SEEK_SET) == -1) {
    taucs_printf("taucs_open: lseek failed\n");
    return NULL;
  }
  /* reading metadata  for every matrix */
  for(i=0; i<h->nmatrices; i++){
    nbytes = read(hs->f,&hs->matrices[i].m, sizeof(int));
    if (nbytes != sizeof(int)) { 
      taucs_printf("taucs_open: Error writing data .\n");
      return NULL;
    }
    nbytes = read(hs->f,&hs->matrices[i].n, sizeof(int));
    if (nbytes != sizeof(int)) { 
      taucs_printf("taucs_open: Error writing data .\n");
      return NULL;
    }
    nbytes = read(hs->f,&hs->matrices[i].flags, sizeof(int));
    if (nbytes != sizeof(int)) { 
      taucs_printf("taucs_open: Error writing data .\n");
      return NULL;
    }
    nbytes = read(hs->f,&hs->matrices[i].offset, sizeof(off_t));
    if (nbytes != sizeof(off_t)) { 
      taucs_printf("taucs_open: Error writing data .\n");
      return NULL;
    }
  }
  return h;
}

taucs_io_handle* taucs_io_open_multifile(char* basename)
{
  int file_id;
  ssize_t nbytes;
  /*int     nmatrices; omer*/
  taucs_io_handle* h;
  mode_t mode;
  taucs_io_handle_multifile* hs;
  int i;
  char filename[256];
  int start_file_index;
  double curr_file_offset;
  int first_size;

  sprintf(filename,"%s.%d",basename,0);
#ifdef OSTYPE_win32
  mode = _O_RDWR | _O_BINARY;
#else
  mode = O_RDWR;
#endif
  file_id = open(filename,mode);

  if (file_id == -1) {
    taucs_printf("taucs_open: Could not open file %s\n",filename);
    return NULL;
  }

  h = (taucs_io_handle*) taucs_malloc(sizeof(taucs_io_handle));
  if (!h) {
    taucs_printf("taucs_open: out of memory (4)\n");
    return NULL;
  }
  h->type      = IO_TYPE_MULTIFILE;
  h->type_specific = (taucs_io_handle_multifile*) taucs_malloc(sizeof(taucs_io_handle_multifile));
  if (!h->type_specific) {
    taucs_printf("taucs_open: out of memory \n");
    taucs_free(h);
    return NULL;
  }
  hs = h->type_specific;
  hs->f[0] = file_id;
  strcpy(hs->basename,basename);
 
  if (lseek(hs->f[0], strlen(TAUCS_FILE_SIGNATURE), SEEK_SET) == -1) {
    taucs_printf("taucs_open: lseek failed\n");
    return NULL;
  }
  nbytes = read(hs->f[0], &h->nmatrices, sizeof(int));
  if (nbytes != sizeof(int)) { 
    taucs_printf("taucs_open: Error read data .\n");
    return NULL;
  }
  nbytes = read(hs->f[0], &hs->last_offset, sizeof(double));
  if (nbytes != sizeof(double)) { 
    taucs_printf("taucs_open: Error read data .\n");
    return NULL;
  }

  hs->matrices = 
      (taucs_io_matrix_multifile*) taucs_malloc((h->nmatrices)* sizeof(taucs_io_matrix_multifile));

  /* open all files before including start */
  start_file_index = (int)floor(((hs->last_offset)/(IO_FILE_RESTRICTION*1024*1024)));
  hs->last_created_file = start_file_index;
  for(i=0;i<=start_file_index;i++){
    sprintf(filename,"%s.%d",hs->basename,i);
    file_id = open(filename,mode);
    if (file_id == -1) {
      taucs_printf("taucs_open: Could not open data file %s\n",filename);
      return NULL;
    }
    hs->f[i] = file_id;
  }
  
  curr_file_offset = hs->last_offset - start_file_index*(IO_FILE_RESTRICTION*1024.0*1024.0);
  /* seek of start offset of data */
  if (lseek(hs->f[start_file_index],(off_t) curr_file_offset, SEEK_SET) == -1) {
    taucs_printf("taucs_open: lseek failed\n");
    return NULL;
  }
  /* reading metadata  for every matrix */
  for(i=0; i<h->nmatrices; i++){
    /* read m */
    if((int)(curr_file_offset+sizeof(int))<IO_FILE_RESTRICTION*1024*1024){
      nbytes = read(hs->f[start_file_index],&hs->matrices[i].m,sizeof(int));
      if (nbytes != sizeof(int)){ 
	taucs_printf("taucs_open: Error in open data .\n");
	return NULL;
      }
      curr_file_offset += (double)sizeof(int);
    }
    else{
      first_size = IO_FILE_RESTRICTION*1024*1024 - (int)curr_file_offset;
      nbytes = read(hs->f[start_file_index],&hs->matrices[i].m,first_size);
      /*if (nbytes != first_size) { omer*/
			if ((int)nbytes != first_size) { 
				taucs_printf("taucs_open: Error in open data .\n");
				return NULL;
      }
      start_file_index++;
      sprintf(filename,"%s.%d",hs->basename,start_file_index);
      file_id = open(filename,mode);
      if (file_id == -1) {
	taucs_printf("taucs_open: Could not open data file %s\n",filename);
	return NULL;
      }
      hs->f[start_file_index] = file_id;
      nbytes = read( hs->f[start_file_index],(char*)&hs->matrices[i].m+first_size,sizeof(int)-first_size);
      if (nbytes != sizeof(int)-first_size){ 
	taucs_printf("taucs_open: Error in open data .\n");
	    return NULL;
      }
      curr_file_offset = (double)(sizeof(int)-first_size);
    }

    /* read n */
    if((int)(curr_file_offset+sizeof(int))<IO_FILE_RESTRICTION*1024*1024){
      nbytes = read(hs->f[start_file_index],&hs->matrices[i].n,sizeof(int));
      if (nbytes != sizeof(int)){ 
	taucs_printf("taucs_open: Error in open data .\n");
	return NULL;
      }
      curr_file_offset += (double)sizeof(int);
    }
    else{
      first_size = IO_FILE_RESTRICTION*1024*1024 - (int)curr_file_offset;
      nbytes = read(hs->f[start_file_index],&hs->matrices[i].n,first_size);
      /*if (nbytes != first_size) { omer*/
			if ((int)nbytes != first_size) { 
				taucs_printf("taucs_open: Error in open data .\n");
				return NULL;
      }
      start_file_index++;
      sprintf(filename,"%s.%d",hs->basename,start_file_index);
      file_id = open(filename,mode);
      if (file_id == -1) {
	taucs_printf("taucs_open: Could not open data file %s\n",filename);
	return NULL;
      }
      hs->f[start_file_index] = file_id;
      nbytes = read( hs->f[start_file_index],(char*)&hs->matrices[i].n+first_size,sizeof(int)-first_size);
      if (nbytes != sizeof(int)-first_size){ 
	taucs_printf("taucs_open: Error in open data .\n");
	    return NULL;
      }
      curr_file_offset = (double)(sizeof(int)-first_size);
    }

    /* read flags */
    if((int)(curr_file_offset+sizeof(int))<IO_FILE_RESTRICTION*1024*1024){
      nbytes = read(hs->f[start_file_index],&hs->matrices[i].flags,sizeof(int));
      if (nbytes != sizeof(int)){ 
	taucs_printf("taucs_open: Error in open data .\n");
	return NULL;
      }
      curr_file_offset += (double)sizeof(int);
    }
    else{
      first_size = IO_FILE_RESTRICTION*1024*1024 - (int)curr_file_offset;
      nbytes = read(hs->f[start_file_index],&hs->matrices[i].flags,first_size);
      /*if (nbytes != first_size) { omer*/
			if ((int)nbytes != first_size) { 
	taucs_printf("taucs_open: Error in open data .\n");
	return NULL;
      }
      start_file_index++;
      sprintf(filename,"%s.%d",hs->basename,start_file_index);
      file_id = open(filename,mode);
      if (file_id == -1) {
	taucs_printf("taucs_open: Could not open data file %s\n",filename);
	return NULL;
      }
      hs->f[start_file_index] = file_id;
      nbytes = read( hs->f[start_file_index],(char*)&hs->matrices[i].flags+first_size,sizeof(int)-first_size);
      if (nbytes != sizeof(int)-first_size){ 
	taucs_printf("taucs_open: Error in open data .\n");
	    return NULL;
      }
      curr_file_offset = (double)(sizeof(int)-first_size);
    }

    /* read offset */
    if((int)(curr_file_offset+sizeof(double))<IO_FILE_RESTRICTION*1024*1024){
      nbytes = read(hs->f[start_file_index],&hs->matrices[i].offset,sizeof(double));
      if (nbytes != sizeof(double)){ 
	taucs_printf("taucs_open: Error in open data .\n");
	return NULL;
      }
      curr_file_offset += (double)sizeof(double);
    }
    else{
      first_size = IO_FILE_RESTRICTION*1024*1024 - (int)curr_file_offset;
      nbytes = read(hs->f[start_file_index],&hs->matrices[i].offset,first_size);
      /*if (nbytes != first_size) { omer*/
			if ((int)nbytes != first_size) { 
				taucs_printf("taucs_open: Error in open data .\n");
				return NULL;
      }
      start_file_index++;
      sprintf(filename,"%s.%d",hs->basename,start_file_index);
      file_id = open(filename,mode);
      if (file_id == -1) {
	taucs_printf("taucs_open: Could not open data file %s\n",filename);
	return NULL;
      }
      hs->f[start_file_index] = file_id;
      nbytes = read( hs->f[start_file_index],(char*)&hs->matrices[i].offset+first_size,sizeof(double)-first_size);
      if (nbytes != sizeof(double)-first_size){ 
	taucs_printf("taucs_open: Error in open data .\n");
	return NULL;
      }
      curr_file_offset = (double)(sizeof(double)-first_size);
    }
  }

  return h;
}

int taucs_io_delete(taucs_io_handle* f)
{
  int i;
  char filename[256];
  int  return_code = 0;

  taucs_printf("taucs_io_delete: starting\n");

  if (f->type == IO_TYPE_SINGLEFILE) {
    taucs_printf("taucs_io_delete: delete only works on multifile; delete singlefile directly\n");
    return -1;
  }
  if (f->type == IO_TYPE_MULTIFILE) {
    taucs_io_handle_multifile* h = ((taucs_io_handle_multifile*) f->type_specific);
    /*taucs_io_matrix_multifile* matrices; omer*/

    for (i=0; i <= h->last_created_file; i++) {
      close((h->f)[i]);
      sprintf(filename,"%s.%d",h->basename,i);
      if (unlink(filename) == -1) {
	taucs_printf("taucs_io_delete: could not delete <%s>\n",filename);
	return_code = -1;
      }
    }

    taucs_free(h->matrices);
  }
  
  taucs_free(f->type_specific);
  taucs_free(f);
  
  taucs_printf("taucs_io_delete: done\n");

  return return_code;
}

/*********************************************************/
/* GET_BASENAME                                          */
/* This routine is used in the ooc_lu code to generate   */
/* additional temporary files                            */
/*********************************************************/

char* taucs_io_get_basename(taucs_io_handle* f)
{
  if (f->type == IO_TYPE_SINGLEFILE) {
    taucs_printf("taucs_io_get_basename: WARNING: only works on multifile\n");
    return NULL;
  }
  if (f->type == IO_TYPE_MULTIFILE) {
    taucs_io_handle_multifile* h = ((taucs_io_handle_multifile*) f->type_specific);
    return h->basename;
  }
  return NULL;
}


/*************************************************************/
/*                                                           */
/*************************************************************/
