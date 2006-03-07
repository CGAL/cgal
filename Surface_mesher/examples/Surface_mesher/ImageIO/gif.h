#ifndef GIF_H
#define GIF_H


#ifdef __cplusplus
extern "C" {
#endif
#include <stdio.h>
#include <imageio/ImageIO.h>
int readGifImage(const char *name,_image *im);

int testGifHeader(char *magic,const char *name);
/** creates an return the file format structure associated with the Gif file format */
PTRIMAGE_FORMAT createGifFormat();

#ifdef __cplusplus
}
#endif

#endif
