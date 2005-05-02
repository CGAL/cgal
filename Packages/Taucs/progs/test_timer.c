/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/* 

TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG BASE
TAUCS_CONFIG TIMING
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END

*/

#include <stdio.h>
#include "taucs.h"

int main()
{
  double x,y;
  int    i,j;
  double w1,w2,c1,c2;
  
  w1 = taucs_wtime();
  c1 = taucs_ctime();

  printf("this may take a while\n");
  x = 1.0 + 1e-10;
  y = 1.0 + 1e-10;
  for (j=0; j<32; j++) {
    for (i=0; i<1000000; i++) x = x*y;
    for (i=0; i<1000000; i++) x = x/y;
  }
  printf("ignore %.4e %.4e\n",x,y);
  
  w2 = taucs_wtime();
  c2 = taucs_ctime();

  if (w1==w2) {
    printf("test failed, taucs_wtime did not advance at all\n");
    return 1;
  }
  if (c1==c2) {
    printf("test failed, taucs_ctime did not advance at all\n");
    return 1;
  }

  printf("test succeeded, wtime=%.2e seconds, ctime=%.2e seconds\n",
	  w2-w1,c2-c1);
  return 0;
}
