#ifndef PNM_H
#define PNM_H

#include <stdio.h>

#include "imageio/ImageIO.h"


int readPpmImage(const char *name,_image *im);
int writePpmImage(char *name, _image *im);
int readPgmAsciiImage(const char *name,_image *im);
int readPgmImage(const char *name,_image *im);
int writePgmImage(char *name,  _image *im);
int testPgmAsciiHeader(char *magic,const char *name);
int testPgmHeader(char *magic,const char *name);
int testPpmHeader(char *magic,const char *name);
PTRIMAGE_FORMAT createPgmFormat();
PTRIMAGE_FORMAT createPgmAscIIFormat();
PTRIMAGE_FORMAT createPpmFormat();


#endif
