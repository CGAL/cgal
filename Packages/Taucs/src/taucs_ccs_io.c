/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

/*#pragma lang +C*/

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <math.h>

#include "taucs.h"

#ifdef OSTYPE_win32
#include <io.h> /*_telli64, _lseeki64*/
#else
#include <unistd.h>
#endif

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>


/*********************************************************/
/* read binary                                           */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL

taucs_ccs_matrix* 
taucs_ccs_read_binary(char* filename)
{
  taucs_ccs_matrix* A = NULL; /* warning*/
  int  nrows,ncols,flags,j;/*nnz, omer*/
  int     f;
  ssize_t bytes_read;
  int*    colptr;

  taucs_printf("taucs_ccs_binary: reading binary matrix %s\n",filename);
  
#ifdef OSTYPE_win32
  f = open(filename,_O_RDONLY |_O_BINARY);
#else
  f = open(filename,O_RDONLY);
#endif
  /*f = open(filename,O_RDONLY);*/

  bytes_read = read(f,&nrows,sizeof(int));
  bytes_read = read(f,&ncols,sizeof(int));
  bytes_read = read(f,&flags,sizeof(int));

  taucs_printf("\t%d-by-%d, flags = %08x\n",nrows,ncols,flags);
  taucs_printf("\t%d-by-%d, flags = %d  \n",nrows,ncols,flags);

  colptr = (int*) taucs_malloc((ncols+1) * sizeof(int));
  assert(colptr);
  
  bytes_read = read(f,colptr,(ncols+1)*sizeof(int));

  taucs_printf("colptr = [");
  for(j=0; j<min(ncols-1,10); j++)
    taucs_printf("%d,",colptr[j]);
  taucs_printf("...,%d]\n",colptr[ncols]);

	if ( 0 ) /* we need this so that we have 'else if' in each type */
	{}
#ifdef TAUCS_DOUBLE_IN_BUILD
  else if (flags & TAUCS_DOUBLE) {
    A = taucs_dccs_create(nrows,ncols,colptr[ncols]);
    if (!A) return NULL;
    bytes_read = read(f,A->rowind,colptr[ncols]*sizeof(int));
    bytes_read = read(f,A->values.d,colptr[ncols]*sizeof(taucs_double));
  }
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  else if (flags & TAUCS_SINGLE) {
    A = taucs_sccs_create(nrows,ncols,colptr[ncols]);
    if (!A) return NULL;
    bytes_read = read(f,A->rowind,colptr[ncols]*sizeof(int));
    bytes_read = read(f,A->values.s,colptr[ncols]*sizeof(taucs_single));
  }
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  else if (flags & TAUCS_DCOMPLEX) {
    A = taucs_zccs_create(nrows,ncols,colptr[ncols]);
    if (!A) return NULL;
    bytes_read = read(f,A->rowind,colptr[ncols]*sizeof(int));
    bytes_read = read(f,A->values.z,colptr[ncols]*sizeof(taucs_dcomplex));
  }
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  else if (flags & TAUCS_SCOMPLEX) {
    A = taucs_cccs_create(nrows,ncols,colptr[ncols]);
    if (!A) return NULL;
    bytes_read = read(f,A->rowind,colptr[ncols]*sizeof(int));
    bytes_read = read(f,A->values.c,colptr[ncols]*sizeof(taucs_scomplex));
  }
#endif
	else {
    assert(0);
  }

  A->flags = flags;
  
  for (j=0; j<=ncols; j++) (A->colptr)[j] = colptr[j];

  taucs_free(colptr);

  close(f);

  taucs_printf("taucs_ccs_read_binary: done reading\n");

  return A;
}

/*********************************************************/
/* read hb                                               */
/*********************************************************/


taucs_ccs_matrix* 
taucs_ccs_read_hb(char* filename,int flags)
{
  taucs_ccs_matrix* A = NULL;
  int  nrows,ncols,nnz,j;
  char fname[256];
  char type[3];
  
  for (j=0; j<256; j++) fname[j] = ' ';
  strcpy(fname,filename);

  taucs_printf("taucs_ccs_read_hb: reading HB matrix %s\n",filename);

  ireadhb_(fname,type,&nrows,&ncols,&nnz);

  if (type[0] == 'p' || type[0] == 'P') {

		if ( 0 ); /* we need this so that we have 'else if' in each type */
#ifdef TAUCS_DOUBLE_IN_BUILD
		else if (flags & TAUCS_DOUBLE) {
      A = taucs_dccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      dreadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.d/*taucs_values*/);
    }
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
		else if (flags & TAUCS_SINGLE) {
      A = taucs_sccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      sreadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.s/*taucs_values*/);
    }
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
		else if (flags & TAUCS_DCOMPLEX) {
      A = taucs_zccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      zreadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.z/*taucs_values*/);
    }
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
		else if (flags & TAUCS_SCOMPLEX) {
      A = taucs_cccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      creadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.c/*taucs_values*/);
    }
#endif
    else {
      assert(0);
    }
  }

  if (type[0] == 'r' || type[0] == 'R') {
		if ( 0 ); /* we need this so that we have 'else if' in each type */
#ifdef TAUCS_DOUBLE_IN_BUILD
		else if (flags & TAUCS_DOUBLE) {
      A = taucs_dccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      dreadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.d/*taucs_values*/);
    }
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
		else if (flags & TAUCS_SINGLE) {
      A = taucs_sccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      sreadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.s/*taucs_values*/);
    }
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
		else if (flags & TAUCS_DCOMPLEX) {
      taucs_printf("taucs_ccs_read_hb: warning: requested a complex type, matrix is real\n");
      A = taucs_dccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      dreadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.d/*taucs_values*/);
    }
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
		else if (flags & TAUCS_SCOMPLEX) {
      taucs_printf("taucs_ccs_read_hb: warning: requested a complex type, matrix is real\n");
      A = taucs_sccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      sreadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.s/*taucs_values*/);
    }
#endif
    else {
      assert(0);
    }
  }

  if (type[0] == 'c' || type[0] == 'C') {
		if ( 0 ); /* we need this so that we have 'else if' in each type */
#ifdef TAUCS_DCOMPLEX_IN_BUILD
		else if (flags & TAUCS_DCOMPLEX) {
      A = taucs_zccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      zreadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.z/*taucs_values*/);
    }
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
		else if (flags & TAUCS_SCOMPLEX) {
      taucs_printf("taucs_ccs_read_hb: warning: requested a complex type, matrix is real\n");
      A = taucs_cccs_create(nrows,ncols,nnz);
      if (!A) return NULL;
      creadhb_(fname,&nrows,&ncols,&nnz,
	       /*A->colptr,A->rowind,A->values); omer*/
				 A->colptr,A->rowind,A->values.c/*taucs_values*/);
    }
#endif
    else {
      assert(0);
    }
  }

  if (type[1] == 's' || type[1] == 'S')
    A->flags |= TAUCS_SYMMETRIC | TAUCS_LOWER;
  if (type[1] == 'h' || type[1] == 'H')
    A->flags |= TAUCS_HERMITIAN | TAUCS_LOWER;

  /* make indices 0-based */
  for (j=0; j<=ncols; j++) ((A->colptr)[j])--;
  for (j=0; j<nnz;    j++) ((A->rowind)[j])--;

  taucs_printf("taucs_ccs_read_hb: done reading\n");

  return A;
}


/*********************************************************/
/* write ijv                                             */
/*********************************************************/


int
taucs_ccs_write_ijv(taucs_ccs_matrix* m, char* ijvfilename)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (m->flags & TAUCS_DOUBLE)
    return taucs_dccs_write_ijv(m,ijvfilename);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (m->flags & TAUCS_SINGLE)
    return taucs_sccs_write_ijv(m,ijvfilename);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (m->flags & TAUCS_DCOMPLEX)
    return taucs_zccs_write_ijv(m,ijvfilename);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (m->flags & TAUCS_SCOMPLEX)
    return taucs_cccs_write_ijv(m,ijvfilename);
#endif
  
  assert(0);
  /*added omer*/
  return -1;
}

#endif /* TAUCS_CORE_GENERAL */


#ifndef TAUCS_CORE_GENERAL
int
taucs_dtl(ccs_write_ijv)(taucs_ccs_matrix* m, 
			 char* ijvfilename)
{
  int i,ip,j,n;
  taucs_datatype Aij;
  FILE* f;

  f = fopen(ijvfilename , "w");

  if (f == NULL) {
    taucs_printf("taucs_ccs_write_ijv: could not open ijv file %s\n",ijvfilename);
    return -1;
  }

  n = m->n;
  
  for (j=0; j<n; j++) {
    for (ip = (m->colptr)[j]; ip < (m->colptr[j+1]); ip++) {
      i   = (m->rowind)[ip];
      Aij = (m->taucs_values)[ip];

#ifdef TAUCS_CORE_DOUBLE
      fprintf(f,"%d %d %0.17e\n",i+1,j+1,Aij);
      if (i != j && ((m->flags) & TAUCS_SYMMETRIC))
	fprintf(f,"%d %d %0.17e\n",j+1,i+1,Aij);
#endif

#ifdef TAUCS_CORE_SINGLE
      fprintf(f,"%d %d %0.9e\n",i+1,j+1,Aij);
      if (i != j && ((m->flags) & TAUCS_SYMMETRIC))
	fprintf(f,"%d %d %0.9e\n",j+1,i+1,Aij);
#endif
      
#ifdef TAUCS_CORE_DCOMPLEX
      fprintf(f,"%d %d %0.17e+%0.17ei\n",i+1,j+1,taucs_re(Aij),taucs_im(Aij));
      if (i != j && ((m->flags) & TAUCS_SYMMETRIC))
	fprintf(f,"%d %d %0.17e+%0.17ei\n",j+1,i+1,taucs_re(Aij),taucs_re(Aij));
#endif
      
#ifdef TAUCS_CORE_SCOMPLEX
      fprintf(f,"%d %d %0.9e+%0.9ei\n",i+1,j+1,taucs_re(Aij),taucs_im(Aij));
      if (i != j && ((m->flags) & TAUCS_SYMMETRIC))
	fprintf(f,"%d %d %0.9e+%0.9ei\n",j+1,i+1,taucs_re(Aij),taucs_im(Aij));
#endif      

    }
  }

  fclose(f);

  return 0;
} 

#endif /*#ifndef TAUCS_CORE_GENERAL*/

/*********************************************************/
/* read ijv                                              */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL

taucs_ccs_matrix* 
taucs_ccs_read_ijv(char* ijvfilename,int flags)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (flags & TAUCS_DOUBLE)
    return taucs_dccs_read_ijv(ijvfilename,flags);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (flags & TAUCS_SINGLE)
    return taucs_sccs_read_ijv(ijvfilename,flags);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (flags & TAUCS_DCOMPLEX)
    return taucs_zccs_read_ijv(ijvfilename,flags);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (flags & TAUCS_SCOMPLEX)
    return taucs_cccs_read_ijv(ijvfilename,flags);
#endif
  
  assert(0);
  /*added omer*/
  return NULL;
}

#endif /* TAUCS_CORE_GENERAL */

#ifndef TAUCS_CORE_GENERAL

taucs_ccs_matrix* 
taucs_dtl(ccs_read_ijv)(char* ijvfilename,int flags)
{
  FILE* f;
  taucs_ccs_matrix*  m;
  int*    clen; 
  int*    is; 
  int*    js;
  taucs_datatype* vs;
  int ncols, nrows, nnz;
  int i,j,k,n;
  double         di,dj;
  taucs_datatype dv;

  f = fopen (ijvfilename , "r");

  if (f == NULL) {
    taucs_printf("taucs_ccs_read_ijv: could not open ijv file %s\n",ijvfilename);
    return NULL;
  }

  n = 10000;
  is = (int*)    taucs_malloc(n*sizeof(int));
  js = (int*)    taucs_malloc(n*sizeof(int));
  vs = (taucs_datatype*) taucs_malloc(n*sizeof(taucs_datatype));
  if (!is || !js || !vs) {
    taucs_printf("symccs_read_ijv: out of memory\n");
    taucs_free(is); taucs_free(js); taucs_free(vs); 
    return NULL; 
  }

  nnz = 0;
  nrows = ncols = 0;
  while (!feof(f)) {
    if (nnz == n) {
      n = (int) ( 1.25 * (double) n);
      taucs_printf("taucs_ccs_read_ijv: allocating %d ijv's\n",n);
      is = (int*)    taucs_realloc(is,n*sizeof(int));
      js = (int*)    taucs_realloc(js,n*sizeof(int));
      vs = (taucs_datatype*) taucs_realloc(vs,n*sizeof(taucs_datatype));
      if (!is || !js || !vs) { 
	taucs_printf("taucs_ccs_read_ijv: out of memory\n");
	taucs_free(is); taucs_free(js); taucs_free(vs); 
	return NULL; 
      }
    }

#ifdef TAUCS_CORE_DOUBLE
    if (fscanf(f, "%lg %lg %lg", &di, &dj, &dv) != 3) break;
#endif

#ifdef TAUCS_CORE_SINGLE
    if (fscanf(f, "%lg %lg %g", &di, &dj, &dv) != 3) break;
#endif

#ifdef TAUCS_CORE_COMPLEX
    {
      taucs_real_datatype dv_i;
      taucs_real_datatype dv_r;

#ifdef TAUCS_CORE_DCOMPLEX
      if (fscanf(f, "%lg %lg %lg+%lgi", &di, &dj, &dv_r,&dv_i) != 4) break;
#endif
#ifdef TAUCS_CORE_SCOMPLEX
      if (fscanf(f, "%lg %lg %g+%gi", &di, &dj, &dv_r, &dv_i) != 4) break;
#endif
      dv = taucs_complex_create(dv_r,dv_i);
    }
#endif

    is[nnz] = (int)di; js[nnz] = (int)dj; vs[nnz] = dv;/*omer*/
    /* we read the lower part */
    if ((flags & TAUCS_SYMMETRIC) && is[nnz] < js[nnz]) continue; 
    if ((flags & TAUCS_HERMITIAN) && is[nnz] < js[nnz]) continue; 
    nrows = max(is[nnz],nrows);
    ncols = max(js[nnz],ncols);
    nnz++;
   }

  fclose ( f );

  m = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!m) { 
    taucs_printf("taucs_ccs_read_ijv: out of memory\n");
    taucs_free(is); taucs_free(js); taucs_free(vs); 
    return NULL; 
  }
  m->n      = nrows;
  m->m      = ncols;
  m->flags  = 0;
  if (flags & TAUCS_SYMMETRIC) 
    m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER;
  if (flags & TAUCS_HERMITIAN) 
    m->flags  = TAUCS_HERMITIAN | TAUCS_LOWER;

#ifdef TAUCS_CORE_DOUBLE
  m->flags |= TAUCS_DOUBLE;
#endif

#ifdef TAUCS_CORE_SINGLE
  m->flags |= TAUCS_SINGLE;
#endif

#ifdef TAUCS_CORE_DCOMPLEX
  m->flags |= TAUCS_DCOMPLEX;
#endif

#ifdef TAUCS_CORE_SCOMPLEX
  m->flags |= TAUCS_SCOMPLEX;
#endif

  clen      = (int*)    taucs_malloc((ncols+1) * sizeof(int));
  m->colptr = (int*)    taucs_malloc((ncols+1) * sizeof(int));
  m->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  m->taucs_values = (taucs_datatype*) taucs_malloc(nnz * sizeof(taucs_datatype));
  if (!clen || !(m->colptr) || !(m->rowind) || !(m->rowind)) {
    taucs_printf("taucs_ccs_read_ijv: out of memory: ncols=%d nnz=%d\n",ncols,nnz);
    taucs_free(clen); taucs_free(m->colptr); taucs_free(m->rowind); 
    taucs_free(m->taucs_values);
    taucs_free (m); taucs_free(is); taucs_free(js); taucs_free(vs); 
    return NULL; 
  }

  for (j=0; j<ncols; j++) clen[j] = 0;
  for (k=0; k<nnz; k++) {
    i = is[k] - 1; /* make it 1-based */
    j = js[k] - 1; /* make it 1-based */
    ( clen[j] )++;
  }
  /* just check */
  k = 0;
  for (j=0; j<ncols; j++) 
    k += clen[j];
  assert(k == nnz);

  /* now compute column pointers */
  
  k = 0;
  for (j=0; j<ncols; j++) {
    int tmp;
    tmp =  clen[j];
    clen[j] = (m->colptr[j]) = k;
    k += tmp;
  }
  clen[ncols] = (m->colptr[ncols]) = k;
  assert(clen[ncols] == nnz);
  
  /* now read matrix into data structure */

  for (k=0; k<nnz; k++) {
    i = is[k] - 1; /* make it 1-based */
    j = js[k] - 1; /* make it 1-based */
    assert(i < nrows);
    assert(j < ncols);
    (m->taucs_values)[ clen[j] ] = vs[k];
    (m->rowind)[ clen[j] ] = i;
    clen[j] ++;
  }
  
  taucs_free(clen);
  taucs_free(vs);
  taucs_free(js);
  taucs_free(is);
  
  taucs_printf("taucs_ccs_read_ijv: read %s, n=%d\n",ijvfilename,m->n);

  return m;
} 

#endif /*#ifndef TAUCS_CORE_GENERAL*/

/*********************************************************/
/* read mtx                                              */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL

taucs_ccs_matrix* 
taucs_ccs_read_mtx(char* mtxfilename,int flags)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (flags & TAUCS_DOUBLE)
    return taucs_dccs_read_mtx(mtxfilename,flags);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (flags & TAUCS_SINGLE)
    return taucs_sccs_read_mtx(mtxfilename,flags);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (flags & TAUCS_DCOMPLEX)
    return taucs_zccs_read_mtx(mtxfilename,flags);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (flags & TAUCS_SCOMPLEX)
    return taucs_cccs_read_mtx(mtxfilename,flags);
#endif  
	
  assert(0);
  /*added omer*/
  return NULL;
}

#endif /* TAUCS_CORE_GENERAL */

#ifndef TAUCS_CORE_GENERAL

taucs_ccs_matrix* 
taucs_dtl(ccs_read_mtx)(char* filename,int flags)
{
  FILE* f;
  taucs_ccs_matrix*  m;
  int*    clen; 
  int*    is; 
  int*    js;
  taucs_datatype* vs;
  int ncols, nrows, nnz;
  int i,j,k,n;
  double di,dj;
  taucs_datatype dv;

  f = fopen (filename , "r");

  if (f == NULL) {
    taucs_printf("taucs_ccs_read_mtx: could not open mtx file %s\n",filename);
    return NULL;
  }

  if (fscanf(f, "%d %d %d", &nrows, &ncols, &nnz) != 3) {
    taucs_printf("taucs_ccs_read_mtx: wrong header\n");
    return NULL;
  }

  n = 10000;
  is = (int*)    taucs_malloc(n*sizeof(int));
  js = (int*)    taucs_malloc(n*sizeof(int));
  vs = (taucs_datatype*) taucs_malloc(n*sizeof(taucs_datatype));
  if (!is || !js || !vs) {
    taucs_printf("taucs_ccs_read_mtx: out of memory\n");
    taucs_free(is); taucs_free(js); taucs_free(vs); 
    return NULL; 
  }

  nnz = 0;
  nrows = ncols = 0;
  while (!feof(f)) {
    if (nnz == n) {
      n = (int) ( 1.25 * (double) n);
      taucs_printf("taucs_ccs_read_mtx: allocating %d ijv's\n",n);
      is = (int*)    taucs_realloc(is,n*sizeof(int));
      js = (int*)    taucs_realloc(js,n*sizeof(int));
      vs = (taucs_datatype*) taucs_realloc(vs,n*sizeof(taucs_datatype));
      if (!is || !js || !vs) { 
	taucs_printf("taucs_ccs_read_mtx: out of memory\n");
	taucs_free(is); taucs_free(js); taucs_free(vs); 
	return NULL; 
      }
    }

#ifdef TAUCS_CORE_DOUBLE
    if (fscanf(f, "%lg %lg %lg", &di, &dj, &dv) != 3) break;
#endif

#ifdef TAUCS_CORE_SINGLE
    if (fscanf(f, "%lg %lg %g", &di, &dj, &dv) != 3) break;
#endif

#ifdef TAUCS_CORE_COMPLEX
    {
      taucs_real_datatype dv_i;
      taucs_real_datatype dv_r;
#ifdef TAUCS_CORE_DCOMPLEX
      if (fscanf(f, "%lg %lg %lg+%lgi", &di, &dj, &dv_r,&dv_i) != 4) break;
#endif
#ifdef TAUCS_CORE_SCOMPLEX
      if (fscanf(f, "%lg %lg %g+%gi", &di, &dj, &dv_r,&dv_i) != 4) break;
#endif
      dv = taucs_complex_create(dv_r,dv_i);
    }
#endif

    is[nnz] = (int)di; js[nnz] = (int)dj; vs[nnz] = dv;/*omer*/
    /* upper or lower might be stored, we use lower */
    if ((flags & TAUCS_SYMMETRIC) && is[nnz] < js[nnz]) {
      int t = is[nnz];
      is[nnz] = js[nnz];
      js[nnz] = t;
    }

    if (flags & TAUCS_PATTERN) {
#ifdef TAUCS_CORE_DOUBLE
      if (is[nnz] == js[nnz]) vs[nnz] = (double) (nrows+1);
      else                    vs[nnz] = -1.0;
#endif

#ifdef TAUCS_CORE_SINGLE
      if (is[nnz] == js[nnz]) vs[nnz] = (float) (nrows+1);
      else                    vs[nnz] = -1.0;
#endif

#ifdef TAUCS_CORE_DCOMPEX
      assert(0);
#endif

#ifdef TAUCS_CORE_SCOMPLEX
      assert(0);
#endif
    }
    nrows = max(is[nnz],nrows);
    ncols = max(js[nnz],ncols);
    nnz++;
   }

  fclose ( f );

  m = (taucs_ccs_matrix*) taucs_malloc(sizeof(taucs_ccs_matrix));
  if (!m) { 
    taucs_printf("taucs_ccs_read_mtx: out of memory\n");
    taucs_free(is); taucs_free(js); taucs_free(vs); 
    return NULL; 
  }
  m->n      = nrows;
  m->m      = ncols;
  if (flags & TAUCS_SYMMETRIC) 
    m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER;
  else
    m->flags  = 0;

#ifdef TAUCS_CORE_DOUBLE
  m->flags |= TAUCS_DOUBLE;
#endif

#ifdef TAUCS_CORE_SINGLE
  m->flags |= TAUCS_SINGLE;
#endif

#ifdef TAUCS_CORE_DCOMPLEX
  m->flags |= TAUCS_DCOMPLEX;
#endif

#ifdef TAUCS_CORE_SCOMPLEX
  m->flags |= TAUCS_SCOMPLEX;
#endif

  clen      = (int*)    taucs_malloc((ncols+1) * sizeof(int));
  m->colptr = (int*)    taucs_malloc((ncols+1) * sizeof(int));
  m->rowind = (int*)    taucs_malloc(nnz       * sizeof(int));
  m->taucs_values = (taucs_datatype*) taucs_malloc(nnz * sizeof(taucs_datatype));
  if (!clen || !(m->colptr) || !(m->rowind) || !(m->rowind)) {
    taucs_printf("taucs_ccs_read_mtx: out of memory: ncols=%d nnz=%d\n",ncols,nnz);
    taucs_free(clen); taucs_free(m->colptr); taucs_free(m->rowind); 
    taucs_free(m->taucs_values);
    taucs_free (m); taucs_free(is); taucs_free(js); taucs_free(vs); 
    return NULL; 
  }

  for (j=0; j<ncols; j++) clen[j] = 0;
  for (k=0; k<nnz; k++) {
    i = is[k] - 1; /* make it 1-based */
    j = js[k] - 1; /* make it 1-based */
    ( clen[j] )++;
  }
  /* just check */
  k = 0;
  for (j=0; j<ncols; j++) 
    k += clen[j];
  assert(k == nnz);

  /* now compute column pointers */
  
  k = 0;
  for (j=0; j<ncols; j++) {
    int tmp;
    tmp =  clen[j];
    clen[j] = (m->colptr[j]) = k;
    k += tmp;
  }
  clen[ncols] = (m->colptr[ncols]) = k;
  assert(clen[ncols] == nnz);
  
  /* now read matrix into data structure */

  for (k=0; k<nnz; k++) {
    i = is[k] - 1; /* make it 1-based */
    j = js[k] - 1; /* make it 1-based */
    assert(i < nrows);
    assert(j < ncols);
    (m->taucs_values)[ clen[j] ] = vs[k];
    (m->rowind)[ clen[j] ] = i;
    clen[j] ++;
  }
  
  taucs_free(clen);
  taucs_free(vs);
  taucs_free(js);
  taucs_free(is);
  
  taucs_printf("taucs_ccs_read_mtx: read %s, n=%d\n",filename,m->n);

  return m;
} 

#endif /* #ifndef TAUCS_CORE_GENERAL */

/*********************************************************/
/* read ccs                                              */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL

taucs_ccs_matrix* 
taucs_ccs_read_ccs(char* ccsfilename,int flags)
{

#ifdef TAUCS_DOUBLE_IN_BUILD
  if (flags & TAUCS_DOUBLE)
    return taucs_dccs_read_ccs(ccsfilename,flags);
#endif

#ifdef TAUCS_SINGLE_IN_BUILD
  if (flags & TAUCS_SINGLE)
    return taucs_sccs_read_ccs(ccsfilename,flags);
#endif

#ifdef TAUCS_DCOMPLEX_IN_BUILD
  if (flags & TAUCS_DCOMPLEX)
    return taucs_zccs_read_ccs(ccsfilename,flags);
#endif

#ifdef TAUCS_SCOMPLEX_IN_BUILD
  if (flags & TAUCS_SCOMPLEX)
    return taucs_cccs_read_ccs(ccsfilename,flags);
#endif  
	
  assert(0);
  /*added omer*/
  return NULL;
}

#endif /* TAUCS_CORE_GENERAL */

#ifndef TAUCS_CORE_GENERAL

taucs_ccs_matrix* 
taucs_dtl(ccs_read_ccs)(char* filename,int flags)
{
  FILE* f;
  taucs_ccs_matrix*  m;

  /*
  int*    clen; 
  int*    is; 
  int*    js;
  taucs_datatype* vs;
  int ncols, nrows, nnz;
  int i,ip,j,k,n;
  */
  /* taucs_datatype dv;*/
  /* double         di,dj;*/

  int i,ip,j,N,*pointers;

  f = fopen(filename ,"r");

  if (f == NULL) {
    taucs_printf("taucs_ccs_read_ccs: could not open ccs file %s\n",filename);
    return NULL;
  }

  fscanf(f,"%d",&N);

  pointers = (int*) taucs_malloc((N+1)*sizeof(int));
  for(i=0; i<N+1; ++i) {
    fscanf(f,"%d",&pointers[i]);
  }

  m = taucs_dtl(ccs_create)(N, N, pointers[N]);
  for (i=0; i<=N; i++) (m->colptr)[i] = pointers[i];

  for(i=0; i<pointers[N]; ++i)
    fscanf(f,"%d",(m->rowind)+i);

#ifdef TAUCS_CORE_DOUBLE  
  for(i=0; i<pointers[N]; ++i)
    fscanf(f,"%lg",(m->taucs_values)+i);
#endif
  
#ifdef TAUCS_CORE_SINGLE  
  for(i=0; i<pointers[N]; ++i)
    fscanf(f,"%g",(m->taucs_values)+i);
#endif
  
#ifdef TAUCS_CORE_DCOMPLEX  
  for(i=0; i<pointers[N]; ++i) {
    taucs_real_datatype dv_r;
    taucs_real_datatype dv_i;
    fscanf(f,"%lg+%lgi",&dv_r,&dv_i);
    (m->taucs_values)[i] = taucs_complex_create(dv_r,dv_i);
  }
#endif
  
#ifdef TAUCS_CORE_SCOMPLEX
  for(i=0; i<pointers[N]; ++i) {
    taucs_real_datatype dv_r;
    taucs_real_datatype dv_i;
    fscanf(f,"%g+%gi",&dv_r,&dv_i);
    (m->taucs_values)[i] = taucs_complex_create(dv_r,dv_i);
  }
#endif
  
  if (flags & TAUCS_SYMMETRIC) {
    m->flags  = TAUCS_SYMMETRIC | TAUCS_LOWER;
    for (j=0; j<N; j++) {
      for (ip=(m->colptr)[j]; ip<(m->colptr)[j+1]; ip++) {
	i = (m->rowind)[ip];
	assert(i >= j);
      }
    }
  } else
    m->flags  = 0;

#ifdef TAUCS_CORE_DOUBLE
  m->flags |= TAUCS_DOUBLE;
#endif

#ifdef TAUCS_CORE_SINGLE
  m->flags |= TAUCS_SINGLE;
#endif

#ifdef TAUCS_CORE_DCOMPLEX
  m->flags |= TAUCS_DCOMPLEX;
#endif

#ifdef TAUCS_CORE_SCOMPLEX
  m->flags |= TAUCS_SCOMPLEX;
#endif

  taucs_free(pointers);
  
  taucs_printf("taucs_ccs_read_ccs: read %s, n=%d\n",filename,m->n);

  return m;
} 

#endif /*#ifndef TAUCS_CORE_GENERAL*/

/*********************************************************/
/* vector io                                             */
/*********************************************************/

#ifdef TAUCS_CORE_GENERAL

void* 
taucs_vec_read_binary(int n, int flags, char* filename)
{
  void* v = NULL; /* warning */
  /*int   nrows; omer*/
  ssize_t bytes_read;
  int f;

  taucs_printf("taucs_vec_read_binary: reading binary vector %s\n",filename);
  
#ifdef OSTYPE_win32
  f = open(filename,_O_RDONLY |_O_BINARY);
#else
  f = open(filename,O_RDONLY);
#endif
  /*f = open(filename,O_RDONLY);*/

  if (flags & TAUCS_DOUBLE) {
    v = taucs_malloc(n * sizeof(taucs_double));
    if (!v) return NULL;
    bytes_read = read(f,v,n*sizeof(taucs_double));
  } else if (flags & TAUCS_SINGLE) {
    v = taucs_malloc(n * sizeof(taucs_single));
    if (!v) return NULL;
    bytes_read = read(f,v,n*sizeof(taucs_single));
  } else if (flags & TAUCS_DCOMPLEX) {
    v = taucs_malloc(n * sizeof(taucs_dcomplex));
    if (!v) return NULL;
    bytes_read = read(f,v,n*sizeof(taucs_dcomplex));
  } else if (flags & TAUCS_SCOMPLEX) {
    v = taucs_malloc(n * sizeof(taucs_scomplex));
    if (!v) return NULL;
    bytes_read = read(f,v,n*sizeof(taucs_scomplex));
  } else {
    assert(0);
  }

  close(f);

  taucs_printf("taucs_vec_read_binary: done reading\n");

  return v;
}

int
taucs_vec_write_binary(int n, int flags, void* v, char* filename)
{
  /*int   nrows; omer*/
  ssize_t bytes_read;
  int f;

  taucs_printf("taucs_vec_write_binary: writing binary vector %s\n",filename);
  
#ifdef OSTYPE_win32
  f = open(filename,
	   _O_WRONLY | _O_CREAT | _O_BINARY, 
	   _S_IREAD | _S_IWRITE | _S_IEXEC);
#else
  f = open(filename,O_WRONLY | O_CREAT | O_TRUNC, S_IRWXO | S_IRWXG | S_IRWXU);
#endif

  if (flags & TAUCS_DOUBLE) {
    bytes_read = write(f,v,n*sizeof(taucs_double));
  } else if (flags & TAUCS_SINGLE) {
    bytes_read = write(f,v,n*sizeof(taucs_single));
  } else if (flags & TAUCS_DCOMPLEX) {
    bytes_read = write(f,v,n*sizeof(taucs_dcomplex));
  } else if (flags & TAUCS_SCOMPLEX) {
    bytes_read = write(f,v,n*sizeof(taucs_scomplex));
  } else {
    assert(0);
  }

  close(f);

  taucs_printf("taucs_vec_read_binary: done reading\n");

  return 0;
}

#endif /* TAUCS_CORE_GENERAL */

/*********************************************************/
/*                                                       */
/*********************************************************/

