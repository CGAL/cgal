/* *******************************************************

   cc_unique.c

   Reads a text from stdin line by line, eliminates 
   duplicated lines and write the other lines to stdout.

   V1.1   27.07.1998   Lutz Kettner

******************************************************** */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#define MaxLineLength 32000

char buf1Storage[MaxLineLength] = "",
     buf2Storage[MaxLineLength];

int main( void){
    char *buf1 = buf1Storage,
         *buf2 = buf2Storage,
         *tmp;

    while( gets(buf2) != NULL) {
        if ( strcmp( buf1, buf2)) {   /* strings are different - print line */
            puts( buf2);
            /* Swap buffers */
            tmp = buf1; buf1 = buf2; buf2 = tmp;
        }                             /* else ignore string */
    }
    return(0);
}
