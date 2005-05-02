/*********************************************************/
/* TAUCS                                                 */
/* Author: Sivan Toledo                                  */
/*********************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>
#include <math.h>
#include <stdarg.h>

#include "taucs.h"

#ifndef OSTYPE_win32
#include <unistd.h>
#endif

/*********************************************************/
/* logging                                               */
/*********************************************************/

#define LOG_NONE 0
#define LOG_STDERR 1
#define LOG_STDOUT 2
#define LOG_FILE   3

static char log_file_prefix[256];
static int  log_file_type = LOG_NONE;
static int  first_time = 0;

void
taucs_logfile(char* file_prefix)
{
  if (!strcmp(file_prefix,"stderr")) {
    log_file_type = LOG_STDERR;
  } else if (!strcmp(file_prefix,"stdout")) {
    log_file_type = LOG_STDOUT;
  } else if (!strcmp(file_prefix,"none")) {
    log_file_type = LOG_NONE;
  } else {
    strcpy(log_file_prefix,file_prefix);
    log_file_type = LOG_FILE;
    first_time = 1;
  }
}

int
taucs_printf(char *fmt, ...)
{
  static FILE* logf;
  va_list      ap;

  if (log_file_type == LOG_NONE) return 0;

  if (first_time && log_file_type == LOG_FILE) {
    char filename[256];

    sprintf(filename,"%s",log_file_prefix);

    if ((logf = fopen(filename,"w")) == NULL) {
      fprintf(stderr,"could not open log file %s, exiting\n",filename);
      exit(1);
    }
    first_time = 0;
  }

  if (log_file_type == LOG_STDERR) logf = stderr;
  if (log_file_type == LOG_STDOUT) logf = stdout;

  va_start(ap, fmt);

  vfprintf(logf, fmt, ap);

  fflush(logf);

  va_end(ap);

  return 0;
}

/*********************************************************/
/*                                                       */
/*********************************************************/
