#include <string.h>

#include "gis.h" 
#include "inr.h"

/* get a string from a file and discard the ending newline character
   if any */
char *fgetns(char *str, int n,  _image *im );
