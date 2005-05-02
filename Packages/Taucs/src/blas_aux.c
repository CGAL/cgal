/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

int lsame(char* ca, char* cb)
{
  return (tolower(*ca) == tolower(*cb));
}

void xerbla(char* srname, int* info)
{
  fprintf(stderr,"** On entry to %.6s parameter number %d had an illegal value\n",
	  srname,*info);
  fprintf(stdout,"** On entry to %.6s parameter number %d had an illegal value\n",
	  srname,*info);
  exit(1);
}

int lsame_(char* ca, char* cb)
{
  return (tolower(*ca) == tolower(*cb));
}

void xerbla_(char* srname, int* info)
{
  fprintf(stderr,"** On entry to %.6s parameter number %d had an illegal value\n",
	  srname,*info);
  fprintf(stdout,"** On entry to %.6s parameter number %d had an illegal value\n",
	  srname,*info);
  exit(1);
}
