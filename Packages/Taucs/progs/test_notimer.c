/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/* 

TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG BASE
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END

*/

#include <stdio.h>
#include "taucs.h"

int main()
{
  double x,y;
  int    i,j;
  
  if (taucs_wtime() != 0.0) {
    printf("test failed, taucs_wtime was supposed to return 0.0\n");
    return 1;
  }
  if (taucs_ctime() != 0.0) {
    printf("test failed, taucs_ctime was supposed to return 0.0\n");
    return 1;
  }

  printf("this may take a while\n");
  x = 1.0 + 1e-10;
  y = 1.0 + 1e-10;
  for (j=0; j<32; j++) {
    for (i=0; i<1000000; i++) x = x*y;
    for (i=0; i<1000000; i++) x = x/y;
  }
  printf("ignore %.4e %.4e\n",x,y);
  
  if (taucs_wtime() != 0.0) {
    printf("test failed, taucs_wtime was supposed to return 0.0\n");
    return 1;
  }
  if (taucs_ctime() != 0.0) {
    printf("test failed, taucs_ctime was supposed to return 0.0\n");
    return 1;
  }

  return 0;
}
