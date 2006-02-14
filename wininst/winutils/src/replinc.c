/* the idea is shamelessly stolen from LEDA distribution */
/* the implementation is different; no restriction on the file length */
/*
** author        : Dima Pasechnik <dima@cs.uu.nl> December 1999
*/

  /* replace "include" by "!include" */
  /* but leave "!include" intact */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#define BUFLEN 256 /* maximal length of a line 
		      (longer lines will be truncated */

main(int argc, char** argv)
{ 
  char *buf;
  FILE *file, *tmpf;
  int n;
  int i;

  buf = (char*)malloc(BUFLEN);

  if (argc < 2)
  { fprintf(stderr,"usage: replinc 'file'");
    exit(1);
   }

  file = fopen(argv[1],"r");

  if (!file)
  { fprintf(stderr,"replinc: cannot open file '%s'",argv[1]);
    exit(1);
   }

  tmpf = tmpfile();
  for(n = 0; fgets(buf, BUFLEN, file); n++) {
    if (strncmp(buf,"include",7) == 0) putc('!', tmpf);
    fputs(buf, tmpf);
  }
  fclose(file);
  rewind(tmpf);
  file = fopen(argv[1],"w");
  for(i = 0; i < n; i++) {
    fgets(buf, BUFLEN, tmpf);
    fputs(buf, file);
  }
  fclose(file);
  _rmtmp();

  return 0;
}
