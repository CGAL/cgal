/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/
/*
TAUCS_CONFIG_BEGIN
TAUCS_CONFIG_DEFAULT OFF
TAUCS_CONFIG BASE
TAUCS_CONFIG ADVANCED_MEMORY_OPS
TAUCS_CONFIG AD_HOC_TEST
TAUCS_CONFIG_END
*/

#include <stdio.h>
#include <taucs.h>

int main()
{
  int failed = 0;
  double system;
  double available;


  system    = taucs_system_memory_size();
  available = taucs_available_memory_size();

  printf("system    memory = %.0f MBtytes\n",system   /1048576.0);
  printf("available memory = %.0f MBtytes\n",available/1048576.0);
  
  if (system==0.0 || available==0.0) failed=1;

  if (failed) {
    printf("test failed\n");
    return -1;
  } else {
    printf("test succeeded\n");
    return 1;
  }
}
