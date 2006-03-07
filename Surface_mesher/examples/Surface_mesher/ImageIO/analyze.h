#ifndef ANALYZE_H
#define ANALYZE_H

#ifdef _MSC_VER
#pragma warning ( disable : 4068 4786 4081 )
#endif



#ifdef __cplusplus
extern "C" {
#endif
#include <imageio/ImageIO.h>

/* read analyse format header

   return:
   -1: error
   0: success
 */
int readAnalyzeHeader(const char* name,_image *im);

int testAnalyzeHeader(char *magic,const char *name);

/** creates an return the file format structure associated with the Analyze file format */
PTRIMAGE_FORMAT createAnalyzeFormat();

/* 
   return:
   -1: error
    1: success
 */
int writeAnalyze( char *basename, _image* im ) ;


/* 
   return:
   -1: error
    1: success
 */
int writeAnalyzeHeader( const _image* im ) ;



/* 
   return:
   -1: error
    1: success
 */
int writeAnalyzeData( const _image* im ) ;

#ifdef __cplusplus
}
#endif

#endif
