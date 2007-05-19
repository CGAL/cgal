#include <string.h>

#include "gis.h" 
#include "inr.h"

/* get a string from a file and discard the ending newline character
   if any */
char *fgetns(char *str, int n,  _image *im ) {
  char *ret;
  int l;

  memset( str, 0, n );
  ret = ImageIO_gets( im, str, n );

  if(!ret) return NULL;

  l = strlen(str);
  if(l > 0 && str[l-1] == '\n') str[l-1] = '\0';
  return ret;
}
