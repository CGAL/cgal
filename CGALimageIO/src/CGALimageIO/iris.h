#ifndef IRIS_H
#define IRIS_H

#include "imageio/ImageIO.h"



PTRIMAGE_FORMAT createIrisFormat();
int readIrisImage( const char *name, _image *im );
int testIrisHeader(char *magic,const char *name);


#endif
