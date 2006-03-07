#ifndef IRIS_H
#define IRIS_H

#include <imageio/ImageIO.h>


#ifdef __cplusplus
extern "C" {
#endif

PTRIMAGE_FORMAT createIrisFormat();
int readIrisImage( const char *name, _image *im );
int testIrisHeader(char *magic,const char *name);

#ifdef __cplusplus
}
#endif

#endif
