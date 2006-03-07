
#ifndef INR_H
#define INR_H


#include <imageio/ImageIO.h>


#ifdef __cplusplus
extern "C" {
#endif

/** Writes the given inrimage header in an already opened file.
    @param im image descriptor
    @param end image data endianness (END_UNKNOWN to use architecture endianness) */
int _writeInrimageHeader(const _image *im,
			 ENDIANNESS end);

/** Writes the given image body in an already opened file.
    @param im image descriptor */
int _writeInrimageData(const _image *im);

/** read header from an opened inrimage file
    @param im opened inrmage descriptor */
int readInrimageHeader(const char *,_image *im);

/** test if an image file has the Inrimage format, given the magic string
    and filename
    @param magic this string is of size 5 and corresponds to the first 4 characters of the file (terminated with a null character)
    @param filename the file name
    @return -1 if it is not an Inrimage file and 0 otherwise */
int testInrimageHeader(char *magic,const char *filename); 

/** calls _writeInrimageHeader and _writeInrimageData
    @param basename the file name without any extension 
    @param im structure
    @return same as  _writeInrimageHeader*/
int writeInrimage(char *basename,_image *im);
/** creates an return the file format structure associated with the Inrimage file format */
PTRIMAGE_FORMAT createInrimageFormat();

#ifdef __cplusplus
}
#endif

#endif
