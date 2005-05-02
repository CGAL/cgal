/*********************************************************/
/* TAUCS Test Suite - test                               */
/* Author: Omer Meshar                                   */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifndef WIN32
#include <unistd.h>
#include <pthread.h>
#endif
#define TAUCS_CORE_DOUBLE
#include <taucs.h>

int actual_main(int argc, char* argv[]);

int main(int argc, char* argv[])
{
  return actual_main(argc,argv);
}

int actual_main(int argc, char* argv[])
{
	int m = 4,n = 4,nnz = 4, i;
	taucs_ccs_matrix * pMatrix = taucs_ccs_create( m, n, nnz, TAUCS_DOUBLE );
	printf("2");

  pMatrix->colptr[0] = 0;
  pMatrix->colptr[1] = 4;  

	printf("3");
  for ( i = 0; i < 4; i++ )
	{
    pMatrix->rowind[i] = i;
		pMatrix->taucs_values[i] = i;
	}
	printf("4");
	i = taucs_ccs_write_ijv(pMatrix, "test.txt" );
	printf("5");
	taucs_dccs_free(pMatrix);
	return i;
}



