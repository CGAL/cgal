/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/*
TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG MALLOC_STUBS
TAUCS_CONFIG ADVANCED_MEMORY_OPS
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END
*/

#include <stdio.h>
#include <setjmp.h>
#include <signal.h>

#include <taucs.h>

double recursive(int n, int d)
{
  if (d==n) return (double) d;
  else      return (double) d + recursive(n,d+1);
}

int main()
{
  int n;

  taucs_logfile("stdout");

  taucs_maximize_stacksize();

  for (n=1; n<2000000; n *= 2) {
    printf("test_stack: depth %d, starting\n",n);
    recursive(n,0);
    printf("test_stack: depth %d, done\n",n);
  }

  return 0;
}
