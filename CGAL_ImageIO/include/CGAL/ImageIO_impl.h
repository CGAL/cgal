// Copyright (c) 2005, 2006 ASCLEPIOS Project, INRIA Sophia-Antipolis (France)
// Copyright (c) 2007 Geometrica Project, INRIA Sophia-Antipolis (France)
// Copyright (c) 2008 GeometryFactory, Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of the ImageIO Library, and as been adapted for CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#ifdef _MSC_VER
// Suppress deprecated warning for fileno and strdup
#  pragma warning(disable: 4127 4706 4996)

#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <io.h>
#include <stdio.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <stdlib.h>
#include <string.h>

/* formats actuellement lus

   format     |   extension(s)   |  lecture  | ecriture
   INRIMAGE   |  .inr[.gz]       |     X     |     X     -> + .gradient[.gz] + .gradient_direction[.gz]
   GIS        |  .dim, .ima[.gz] |     X     |     X
   ANALYZE    |  .hdr, .img[.gz] |     X     |     X
   PNM        |  .ppm, .pgm      |     X     |     X
   GIF        |  .gif            |     X     |
   BMP        |  .gif            |     X     |
*/
#include <CGAL/ImageIO/inr.h>
#include <CGAL/ImageIO/gif.h>
#include <CGAL/ImageIO/gis.h>
#include <CGAL/ImageIO/pnm.h>
#include <CGAL/ImageIO/bmp.h>
#include <CGAL/ImageIO/iris.h>
#include <CGAL/ImageIO/analyze.h>
#ifdef MINC_FILES
#include <CGAL/ImageIO/mincio.h>
#endif


#ifdef CGAL_USE_ZLIB
#define ZLIB_MAJOR_VERSION ZLIB_VERNUM>>8
#define ZLIB_MINOR_VERSION (ZLIB_VERNUM-(ZLIB_MAJOR_VERSION <<8))
#define CGAL_USE_GZFWRITE  ZLIB_MAJOR_VERSION > 0x12 || ((ZLIB_MAJOR_VERSION == 0x12) && (ZLIB_MINOR_VERSION >= 0x90))
#endif
struct Remove_supported_file_format {
  ~Remove_supported_file_format()
  {
    removeSupportedFileFormat();
  }
};

#ifdef CGAL_HEADER_ONLY

inline PTRIMAGE_FORMAT & get_static_firstFormat()
{
  static PTRIMAGE_FORMAT firstFormat = nullptr;
  return firstFormat;
}

inline PTRIMAGE_FORMAT & get_static_inrimageFormat()
{
  static PTRIMAGE_FORMAT inrimageFormat = nullptr;
  return inrimageFormat;
}

inline Remove_supported_file_format & get_static_rsff()
{
  static Remove_supported_file_format rsff;
  return rsff;
}
// Dummy call to get_static_rsff(), otherwise it would not get instanced
CGAL_UNUSED static Remove_supported_file_format &rsff_dummy_ref = get_static_rsff();


#else // not header-only

/** the first file format is initialized to null */
static PTRIMAGE_FORMAT firstFormat = nullptr;
inline PTRIMAGE_FORMAT & get_static_firstFormat()
{
  return firstFormat;
}

/** the Inrimage file format (default format) is initialized to null */
static PTRIMAGE_FORMAT InrimageFormat = nullptr;
inline PTRIMAGE_FORMAT & get_static_inrimageFormat()
{
  return InrimageFormat;
}

static Remove_supported_file_format rsff;
inline Remove_supported_file_format & get_static_rsff()
{
  return rsff;
}

#endif

/*--------------------------------------------------
 *
 * mimics standard routines
 *
 --------------------------------------------------*/



extern "C" {
  /* default allocation routine */
  static void *(*allocRoutine)(size_t) = 0;
  /* default deallocation routine */
  static void (*deleteRoutine)(void *) = 0;
}

CGAL_INLINE_FUNCTION
void *ImageIO_alloc(size_t s) {
  if(!allocRoutine) allocRoutine = malloc;
  return ( (*allocRoutine)(s) );
}
/* call allocation routine */
CGAL_INLINE_FUNCTION
void ImageIO_free(void *m) {
  if(!deleteRoutine) deleteRoutine = free;
  (*deleteRoutine)(m);
}
/* call deallocation routine */


CGAL_INLINE_FUNCTION
unsigned int ImageIO_limit_len(size_t to_be_read)
{
  return (unsigned int)(std::min)(to_be_read, size_t(1u<<30));
}

/* mimics fwrite() function.
   According to _openWriteImage(), openMode will has one
   of the following value:
   - OM_STD (for stdout)
   - OM_GZ
   - OM_FILE
*/
CGAL_INLINE_FUNCTION
size_t ImageIO_write(const _image *im, const void *buf, size_t len) {
  size_t to_be_written = len;
  std::ptrdiff_t l = -1;
  char *b = (char*)buf;

  switch(im->openMode) {
  default :
  case OM_CLOSE :
    return 0;
  case OM_STD :
#ifdef CGAL_USE_ZLIB
    while ( (to_be_written > 0) && ((l = gzwrite(im->fd, (void *) b, ImageIO_limit_len(to_be_written))) > 0) ) {
      to_be_written -= l;
      b += l;
    }
#else
    while ( (to_be_written > 0) && ((l = fwrite( b, 1, ImageIO_limit_len(to_be_written), im->fd )) > 0) ) {
      to_be_written -= l;
      b += l;
    }
#endif
    return ( len - to_be_written );
#ifdef CGAL_USE_ZLIB
  case OM_GZ :
    while ( (to_be_written > 0) && ((l = gzwrite(im->fd, (void *) b, ImageIO_limit_len(to_be_written))) > 0) ) {
      to_be_written -= l;
      b += l;
    }
    if(l<0)
    {
      int errnum;
      fprintf(stderr, "zlib error: %s\n", gzerror(im->fd, &errnum));
    }
    return ( len - to_be_written );
#if CGAL_USE_GZFWRITE
  case OM_FILE :
    while ( (to_be_written > 0) && ((l = gzfwrite( b, sizeof(char), ImageIO_limit_len(to_be_written), im->fd )) > 0) ) {
      to_be_written -= l;
      b += l;
    }
    return ( len - to_be_written );
#endif // CGAL_USE_GZFWRITE
#else
  case OM_FILE :
    while ( (to_be_written > 0) && ((l = fwrite( b, 1, ImageIO_limit_len(to_be_written), im->fd )) > 0) ) {
      to_be_written -= l;
      b += l;
    }
    return ( len - to_be_written );
#endif
  }
}

CGAL_INLINE_FUNCTION
size_t ImageIO_limit_read(size_t to_be_read)
{
  return (std::min)(to_be_read, size_t(1u<<30));
}

/* mimics fread() function.
   According to _openReadImage(), openMode will has one
   of the following value:
   - OM_STD (for stdin)
   - OM_GZ *or* OM_FILE
*/
CGAL_INLINE_FUNCTION
size_t ImageIO_read(const _image *im, void *buf, size_t len)
{
  size_t to_be_read = len;
  int l = -1;
  char *b = (char*)buf;

  switch(im->openMode) {
  default :
  case OM_CLOSE :
    return 0;
  case OM_STD :
#ifdef CGAL_USE_ZLIB
    while ( (to_be_read > 0) && ((l = gzread(im->fd, (void *) b, ImageIO_limit_len(to_be_read))) > 0) ) {
      to_be_read -= l;
      b += l;
    }
#else
    while ( (to_be_read > 0) && ((l = fread( b, 1, ImageIO_limit_len(to_be_read), im->fd )) > 0) ) {
      to_be_read -= l;
      b += l;
    }
#endif
    return ( len - to_be_read );
#ifdef CGAL_USE_ZLIB
  case OM_GZ :
#if CGAL_USE_GZFWRITE
  case OM_FILE :
#endif// CGAL_USE_GZFWRITE
    while ( (to_be_read > 0) && ((l = gzread(im->fd, (void *) b, ImageIO_limit_len(to_be_read))) > 0) ) {
      to_be_read -= l;
      b += l;
    }
    if(l<0)
    {
      int errnum;
      fprintf(stderr, "zlib error: %s\n", gzerror(im->fd, &errnum));
    }
    return ( len - to_be_read );
#else
  case OM_FILE :
    while ( (to_be_read > 0) && ((l = fread( b, 1, ImageIO_limit_len(to_be_read), im->fd )) > 0) ) {
      to_be_read -= l;
      b += l;
    }
    return ( len - to_be_read );
#endif
  }

  //return 0;
}



/* mimics fgets() function.
   According to _openReadImage(), openMode will has one
   of the following value:
   - OM_STD (for stdout)
   - OM_GZ *or* OM_FILE
*/
CGAL_INLINE_FUNCTION
char *ImageIO_gets( const _image *im, char *str, int size )
{
  char *ret = nullptr;
  switch(im->openMode) {
  default :
  case OM_CLOSE :
    return nullptr;
  case OM_STD :
#ifdef CGAL_USE_ZLIB
    ret = (char *) gzgets(im->fd, str, size );
#else
    ret = fgets(str, size, im->fd);
#endif
    break;
#ifdef CGAL_USE_ZLIB
  case OM_GZ :
#if CGAL_USE_GZFWRITE
  case OM_FILE :
#endif // CGAL_USE_GZFWRITE
    ret = (char *) gzgets(im->fd, str, size);
    break;

#else
  case OM_FILE :
    ret = fgets(str, size, im->fd);
    break;
#endif
  }
  return ret;
}



CGAL_INLINE_FUNCTION
long ImageIO_seek( const _image *im, long offset, int whence ) {
  switch(im->openMode) {
  case OM_CLOSE :
  default :
    return -1;
#ifdef CGAL_USE_ZLIB
  case OM_GZ:
#if CGAL_USE_GZFWRITE
  case OM_FILE:
#endif //CGAL_USE_GZFWRITE
    return gzseek(im->fd, offset, whence );
#else
  case OM_FILE:
    return fseek( (FILE*)im->fd, offset, whence );
#endif
  }
}

/* return non 0 in case of error
 */
CGAL_INLINE_FUNCTION
int ImageIO_error( const _image *im )
{
  switch(im->openMode) {
  case OM_CLOSE :
  default :
    return 0;
#ifdef CGAL_USE_ZLIB
  case OM_GZ :
#if CGAL_USE_GZFWRITE
  case OM_FILE :
#endif //CGAL_USE_GZFWRITE
    static int errnum;
    (void)gzerror(im->fd, &errnum);
    return( (errnum != Z_OK) || gzeof(im->fd) );
#else
  case OM_FILE :
    return( ferror( (FILE*)im->fd ) || feof( (FILE*)im->fd ) );
#endif
  }
  //return 0;
}



/* Upon successful completion 0 is returned.

   Closing the standard output with gzclose()
   is necessary as it will flush the pending output.
 */
CGAL_INLINE_FUNCTION
int ImageIO_close( _image* im )
{
  int ret=0;

  switch ( im->openMode ) {
  default :
  case OM_CLOSE :
    break;
#ifdef CGAL_USE_ZLIB
  case OM_GZ :
  case OM_STD :
#if CGAL_USE_GZFWRITE
  case OM_FILE :
#endif//CGAL_USE_GZFWRITE
    ret = gzclose( im->fd );
    break;
#else
  case OM_STD :
    break;
  case OM_FILE :
    ret = fclose( (FILE*)im->fd );
#endif
  }

  im->fd = nullptr;
  im->openMode = OM_CLOSE;

  return ret;
}







/* given an initialized file descriptor and a file name,
   open file from stdin (if name == nullptr, or name == "-", or name == "<"),
   or a standard/gzipped file otherwise (gzipped files are handled assuming
   that it is compiled and linked with zlib).
   openMode will have one of the following value:
   - OM_STD (for stdin)
   - OM_GZ *or* OM_FILE
*/
CGAL_INLINE_FUNCTION
void _openReadImage(_image* im, const char *name) {
  if(im->openMode == OM_CLOSE) {

    /* open from stdin */
    if( name == nullptr || name[0] == '\0'
        || (name[0] == '-' && name[1] == '\0')
        || (name[0] == '<' && name[1] == '\0') ) {
#ifdef CGAL_USE_ZLIB
      im->fd = gzdopen(fileno(stdin), "rb");
#else
      im->fd = fdopen(fileno(stdin), "rb");
#endif
      im->openMode = OM_STD;
    }

    else {
#ifdef CGAL_USE_ZLIB
      im->fd = gzopen(name, "rb");
      if(im->fd) im->openMode = OM_GZ;
#else
      im->fd = fopen(name, "rb");
      if(im->fd) im->openMode = OM_FILE;
#endif

    }

  }
}





/* given an initialized file descriptor and a file name,
   open file from stdout (if name == nullptr, or name == "-", or name == ">"),
   a gzipped pipe (if name got the extension ".gz")
   or a standard file otherwise.
   openMode will have one of the following value:
   - OM_STD (for stdout)
   - OM_GZ
   - OM_FILE
*/
CGAL_INLINE_FUNCTION
void _openWriteImage(_image* im, const char *name)
{
  im->openMode = OM_CLOSE;

  if( name == nullptr || name[0] == '\0'
      || (name[0] == '-' && name[1] == '\0')
      || (name[0] == '>' && name[1] == '\0') ) {

#ifdef CGAL_USE_ZLIB
#if (defined _LINUX_) || (defined _SOLARIS_) || (defined _SGI_)
    im->fd = gzdopen(1, "wb");
#else
    im->fd = gzdopen(fileno(stdout), "wb");
#endif
#else
    im->fd = (_ImageIO_file) stdout;
#endif
    im->openMode = OM_STD;
  }

  else{
#ifdef CGAL_USE_ZLIB

    /* from gzopen() doc:
       ... The mode parameter is as in fopen ("rb" or "wb") but can
       also include a compression level ("wb9") or a strategy: 'f' for
       filtered data as in "wb6f", 'h' for Huffman only compression as
       in "wb1h" ...
       However, a small .gz header will be written ... thus gz(d)open can not
       be used for rgular files.
    */

    if( !strncmp(name+strlen(name)-3, ".gz", 3) )
      {
#ifdef _MSC_VER
        int ffd=_open(name,_O_RDWR | _O_CREAT| _O_TRUNC | _O_BINARY, _S_IREAD|_S_IWRITE);
        im->fd = gzdopen( ffd, "wb" );
#else
        im->fd = gzopen( name, "wb" );
#endif
        im->openMode = OM_GZ;
      }
    else
    {
#if CGAL_USE_GZFWRITE
      im->fd = (_ImageIO_file) gzopen(name, "wb");
      im->openMode = OM_FILE;
#else
      fprintf(stderr, "_openWriteImage: error: zlib version 1.2.9 or later\n"
                      "is required to save in non-compressed files\n");
      return;
#endif// CGAL_USE_GZFWRITE
    }

#else //CGAL_USE_ZLIB
    {
      im->fd = (_ImageIO_file) fopen(name, "wb");
      im->openMode = OM_FILE;
    }
#endif
  }
}



/* set allocation and deallocation routines */
CGAL_INLINE_FUNCTION
void setImageIOAllocationRoutines(ALLOCATION_FUNCTION alloc,
                                  DEALLOCATION_FUNCTION del) {
  if(alloc != nullptr) allocRoutine = alloc;
  if(del != nullptr) deleteRoutine = del;
}





/*--------------------------------------------------
 *
 *
 *
 --------------------------------------------------*/





CGAL_INLINE_FUNCTION
ENDIANNESS _getEndianness()
{
  union {
    unsigned char uc[2];
    unsigned short us;
  } twobytes;
  twobytes.us = 255;
  /* on linux or dec
   */
  if ( twobytes.uc[1] == 0 ) return( END_LITTLE );
  /* on solaris or sgi
   */
  return( END_BIG );
}







/* Allocates and initializes an image descriptor */
CGAL_INLINE_FUNCTION
_image *_initImage() {
  _image *im;

  im = (_image *) ImageIO_alloc(sizeof(_image));
  if ( im == nullptr ) return( im );

  /* default image size is 1*1*1 */
  im->xdim = im->ydim = im->zdim = im->vdim = 1;
  /* default image voxel size is 1.0*1.0*1.0 */
  im->vx = im->vy = im->vz = 1.0;

  /* default image center  is 0 0 0 */
  im->cx = im->cy = im->cz = 0;

  /* default image offset  is 0 0 0 */
  im->tx = im->ty = im->tz = 0.0;

  /* default image rotation  is 0 0 0 */
  im->rx = im->ry = im->rz = 0.0;

  /* no data yet */
  im->data = nullptr;

  /* no file associated to image */
  im->fd = nullptr;
  im->openMode = OM_CLOSE;
  im->endianness = END_UNKNOWN;

  /* unknown data kind
     default is binary
   */
  im->dataMode = DM_BINARY;

  /* no user string */
  im->user = nullptr;
  im->nuser = 0;

  /* unknown word kind */
  im->wdim = 0;
  im->wordKind = WK_UNKNOWN;
  im->vectMode = VM_SCALAR;
  im->sign = SGN_UNKNOWN;
  im->imageFormat = nullptr;

  /** eventually initializes the supported file formats */
  if (get_static_firstFormat()==nullptr)
    initSupportedFileFormat();
  /* return image descriptor */
  return im;
}

CGAL_INLINE_FUNCTION
_image *_createImage(std::size_t x, std::size_t y, std::size_t z, std::size_t v,
                     double vx, double vy, double vz, std::size_t w,
                     WORD_KIND wk, SIGN sgn)
{
  _image *im;

  im = (_image *) ImageIO_alloc(sizeof(_image));
  if ( im == nullptr ) return( im );

  im->xdim = x;
  im->ydim = y;
  im->zdim = z;
  im->vdim = v;
  im->vx = vx;
  im->vy = vy;
  im->vz = vz;

  /* default image center  is 0 0 0 */
  im->cx = im->cy = im->cz = 0;

  /* default image offset  is 0 0 0 */
  im->tx = im->ty = im->tz = 0.0;

  /* default image rotation  is 0 0 0 */
  im->rx = im->ry = im->rz = 0.0;

  /* no data yet */
  im->data = ImageIO_alloc(std::size_t(x)*std::size_t(y)*std::size_t(z)*std::size_t(v)*std::size_t(w));

  /* no file associated to image */
  im->fd = nullptr;
  im->openMode = OM_CLOSE;
  im->endianness = END_UNKNOWN;

  /* unknown data kind
     default is binary
   */
  im->dataMode = DM_BINARY;

  /* no user string */
  im->user = nullptr;
  im->nuser = 0;

  /* unknown word kind */
  im->wdim = w;
  im->wordKind = wk;
  im->vectMode = VM_SCALAR;
  im->sign = sgn;
  im->imageFormat = nullptr;

  /** eventually initializes the supported file formats */
  if (get_static_firstFormat()==nullptr)
    initSupportedFileFormat();
  /* return image descriptor */
  return im;
}

/* return the bounding box of the image */
CGAL_INLINE_FUNCTION
void _get_image_bounding_box(_image* im,
                             double* x_min, double* y_min, double* z_min,
                             double* x_max, double* y_max, double* z_max) {
  *x_min = im->tx;
  *y_min = im->ty;
  *z_min = im->tz;
  *x_max = (double(im->xdim) - 1.0)*(im->vx) + *x_min ;
  *y_max = (double(im->ydim) - 1.0)*(im->vy) + *y_min ;
  *z_max = (double(im->zdim) - 1.0)*(im->vz) + *z_min ;
}

/* Free an image descriptor */
CGAL_INLINE_FUNCTION
void _freeImage(_image *im) {
  unsigned int i;

  if ( im == nullptr ) return;

  /* close image if opened */
  if(im->openMode != OM_CLOSE) ImageIO_close(im);

  /* free data if any */
  if(im->data != nullptr) ImageIO_free(im->data);
  im->data = nullptr;

  /* free user string array if any */
  if( (im->nuser > 0) && (im->user != nullptr) ) {
    for(i = 0; i < im->nuser; i++)
      if ( im->user[i] != nullptr ) ImageIO_free(im->user[i]);
    ImageIO_free(im->user);
  }
  im->nuser = 0;
  im->user = nullptr;

  /* free given descriptor */
  ImageIO_free(im);
}





/* Reads an image from a file and returns an image descriptor or nullptr if
   reading failed.
   Reads from stdin if image name is nullptr. */
CGAL_INLINE_FUNCTION
_image* _readImage(const char *name) {
  _image *im;


  /* read header */
  im = _readImageHeader( name );

  if(im != nullptr && im->openMode != OM_CLOSE) {
    /* read body */
    if(_readImageData(im) < 0) {
      fprintf(stderr, "_readImage: error: invalid data encountered in \'%s\'\n",
              name);
      _freeImage(im);
      return nullptr;
    }
    ImageIO_close(im);
  }

  return im;
}

// raw
CGAL_INLINE_FUNCTION
_image* _readImage_raw(const char *name,
                       const unsigned int rx,
                       const unsigned int ry,
                       const unsigned int rz,
                       const double vx,
                       const double vy,
                       const double vz,
                       const unsigned int offset,
                       const std::size_t wdim,
                       WORD_KIND wk,
                       SIGN sgned
                       )
{
  _image *im = nullptr;
  im = (_image *) ImageIO_alloc(sizeof(_image));
  if ( im == nullptr )
    return nullptr;

  im->xdim = rx;
  im->ydim = ry;
  im->zdim = rz;
  im->vdim = 1;
  im->vx = vx;
  im->vy = vy;
  im->vz = vz;

  // image center
  im->cx = im->cy = im->cz = 0;

  // image offset
  im->tx = im->ty = im->tz = 0.0;

  // image rotation
  im->rx = im->ry = im->rz = 0.0;


  im->fd = nullptr;
  im->openMode = OM_CLOSE;
  im->endianness = END_UNKNOWN;

  im->dataMode = DM_BINARY;

  // no user string
  im->user = nullptr;
  im->nuser = 0;

  // word type (unsigned byte)
  im->wdim = wdim;
  im->wordKind = wk;
  im->vectMode = VM_SCALAR;
  im->sign = sgned;
  im->imageFormat = nullptr;

  // read file
  ::_openReadImage(im, name);
  if(!im->fd) {
    fprintf(stderr, "_readImage_raw: error: unable to open file \'%s\'\n", name);
    _freeImage(im);
    return nullptr;
  }

  // read offset
  if(offset > 0) {
    im->data = ImageIO_alloc(offset+1);
    ImageIO_read(im, im->data, offset);
    ImageIO_free(im->data);
  }
  // allocate memory
  im->data = ImageIO_alloc(rx*ry*rz*wdim);
  if(im->data == nullptr)
    return nullptr;

  // read
  ImageIO_read(im, im->data, rx*ry*rz*wdim);

  ImageIO_close(im);
  /*
    unsigned int i,j,k;
    unsigned char *data = (unsigned char *)im->data;
    int dimxy = rx * ry;
    for(i=0;i<rx;i++)
    for(j=0;j<ry;j++)
    for(k=0;k<rz;k++)
    {
    unsigned char voxel;
    fread(&voxel,1,1,pFile);
    data[k * dimxy + j * rx + i] = voxel;
    }
  */

  return im;
}

/* Reads an image from a file and returns an image descriptor or nullptr if<br>
   reading failed.<br>
   Reads from stdin if image name is nullptr.
   If the image is vectorial, it is uninterlaced. */
CGAL_INLINE_FUNCTION
_image* _readNonInterlacedImage(const char *name) {
  _image *im;

  /* read header */
  im = _readImageHeader(name);

  if(im && im->openMode != OM_CLOSE) {
    /* read scalar image body */
    if(im->vdim == 1) {
      if(_readImageData(im) < 0) {
        fprintf(stderr, "_readImage: error: invalid data encountered in \'%s\'\n",
                name);
        _freeImage(im);
        return nullptr;
      }
    }
    /* read vectorial image body */
    else {
      im->vectMode = VM_NON_INTERLACED;
      if(_readNonInterlacedImageData(im) < 0) {
        fprintf(stderr, "_readImage: error: invalid data encountered in \'%s\'\n",
                name);
        _freeImage(im);
        return nullptr;
      }
    }
    ImageIO_close(im);
   }

  return im;
}











/* Write inrimage given in inr in file name. If file name's suffix is
   .gz, the image is gziped. If file name's suffix is .hdr, the image
   is written in ANALYZE format. If file name is nullptr, image is written
   on stdout */
CGAL_INLINE_FUNCTION
int _writeImage(_image *im, const char *name_to_be_written ) {
  int r = ImageIO_NO_ERROR;
  std::size_t length = 0;
  char *name = nullptr;
  char *baseName = nullptr;

  if ( im == nullptr ) return -1;

  /* different conventions for the standard input
   */
  if ( name_to_be_written == nullptr || name_to_be_written[0] == '\0'
       || (name_to_be_written[0] == '-' && name_to_be_written[1] == '\0')
       || (name_to_be_written[0] == '>' && name_to_be_written[1] == '\0') ) {
    name = nullptr;
  }
  else {
    name = strdup( name_to_be_written );
  }

  initSupportedFileFormat();

  /* what is the wanted format
   */
  if ( name == nullptr ) {
    im->imageFormat = get_static_inrimageFormat();
  } else {
    std::size_t i,extLength;
    PTRIMAGE_FORMAT f;
    char ext[IMAGE_FORMAT_NAME_LENGTH];
    char *ptr;


    /* scan all formats; */
    im->imageFormat=nullptr;
    length=strlen(name);

    for(f=get_static_firstFormat();(f!=nullptr)&& (im->imageFormat==nullptr);f=f->next) {
      /* scan all extensions for that format */
      ptr=&f->fileExtension[0];

      do {
        /* get next file extension */
        i=0;
        for(i=0;((*ptr)!=',' && (*ptr)!='\0');i++,ptr++) {
          ext[i]=(*ptr);
        }
        if ((*ptr)==',') {
          ext[i]='\0';
          ptr++;
        }
        else {
          ext[i]='\0';
        }
        extLength=strlen(ext);

        /* test if the tail of name matches the extension */
        if ( (length > extLength) && (!strcmp( name + length - extLength, ext)) ) {
          im->imageFormat=f;
          /* copy original name and removes extension */
          baseName=strdup(name);
          for(i= length - extLength;i<length;++i)
            baseName[i]='\0';
        }

      } while (((*ptr)!='\0') && (im->imageFormat==nullptr));
    }

    if (!im->imageFormat) {
      fprintf(stderr, "_writeImage: warning : unknown extension in %s = assuming Inrimage\n",name);
      im->imageFormat=get_static_inrimageFormat();
      baseName=strdup(name);
    }
  }


  /* open file descriptor */
  /* _openWriteImage( im, name ) ;



  if(!im->fd) {
     fprintf(stderr, "_writeImage: error: open failed\n");
     if ( name != nullptr ) free( name );
     if ( baseName != nullptr ) free( baseName );
     return ImageIO_OPENING;
  }
  */

  if (im->imageFormat) {

    if (im->imageFormat->writeImage==nullptr) {
      im->imageFormat=get_static_inrimageFormat();
    }

    if ( 0 ) {
      fprintf(stderr, "_writeImage: will write '%s' with '%s' format\n",
              name, im->imageFormat->realName );
    }

    if ((*im->imageFormat->writeImage)(name, im)<0) {
      fprintf(stderr, "_writeImage: error: unable to write \'%s\'\n",
              name);
      r = ImageIO_WRITING_HEADER;
    }

  }



  /* close file descriptor */
  ImageIO_close( im );

  im->fd = nullptr;
  im->openMode = OM_CLOSE;

  if ( baseName != nullptr ) free( baseName );
  if ( name != nullptr ) free( name );

  return r;
}












/* read header from an image file

   if standard input, it's an inrimage
   if not, get a magic string
   and try to find the good format

   if data are in a separate file,
   the header reading procedure will open
   the data file.

   error:
   0  success
   -1 unknown image type
   -2 error while opening
   -3 error while reading header
   -4 error while reading header or data
 */
CGAL_INLINE_FUNCTION
_image *_readImageHeader( const char *name ) {
  int error = 0;
  return( _readImageHeaderAndGetError( name, &error ) );
}



CGAL_INLINE_FUNCTION
_image *_readImageHeaderAndGetError( const char *name_to_be_read, int *error )
{
  _image *im;
  char magic[5];
  char *name = nullptr;
  PTRIMAGE_FORMAT f;
  int res;

  *error = ImageIO_NO_ERROR;

  /* open image file */
  im = _initImage();
  if ( name_to_be_read == nullptr || name_to_be_read[0] == '\0'
       || (name_to_be_read[0] == '-' && name_to_be_read[1] == '\0')
       || (name_to_be_read[0] == '<' && name_to_be_read[1] == '\0') ) {
    name = nullptr;
  }
  else {
    name = strdup( name_to_be_read );
  }


  _openReadImage(im, name);

  if(!im->fd) {
    if(name == nullptr) {
      fprintf(stderr, "_readImageHeaderAndGetError: error: nullptr file name\n");
    }  else {
      fprintf(stderr, "_readImageHeaderAndGetError: error: unable to open file \'%s\'\n", name);
    }
    _freeImage(im);
    *error = ImageIO_OPENING;
    if ( name != nullptr ) free( name );
    return nullptr;
  }

  initSupportedFileFormat();

  /* what is the wanted format ?
     assume that stdin is inrimage
   */
  if(im->openMode == OM_STD) {
    im->imageFormat=get_static_inrimageFormat();
  }
  else {
   /* get magic string for disk files
    */
    ImageIO_read(im, magic, 4);
    magic[4] = '\0';
    ImageIO_seek(im, 0L, SEEK_SET);
    /** test each format */
    for(f=get_static_firstFormat();(f!=nullptr)&& (im->imageFormat==nullptr);f=f->next) {
      /* test if it is the correct format based on magic and file extension */
      if (((*f->testImageFormat)(magic, name)) >=0) {
        im->imageFormat=f;
      }
    }
  }

  if ( im->imageFormat == nullptr ) {
    fprintf(stderr, "_readImageHeaderAndGetError: does not find image format for \'%s\'\n", name);
    ImageIO_close( im );
    _freeImage(im);
    *error = ImageIO_UNKNOWN_TYPE;
    if ( name != nullptr ) free( name );
    return nullptr;
  }

  /* now tests if the header can be read correctly */

  res=(*(im->imageFormat)->readImageHeader)(name,im);
  /* could read header only */
  if (res == 0) {
    if ( name != nullptr ) free( name );
    return( im );
  }
  /* could read header and data */
  else if ( res > 0 ) {
    ImageIO_close(im);
    if ( name != nullptr ) free( name );
    return im;
  }

  /* could not read error : throw error */
  fprintf(stderr, "_readImageHeaderAndGetError: an error occurs when reading image\n" );
  if ( name == nullptr || im->openMode == OM_STD) {
    fprintf(stderr, "\t from \'standard input\'" );
  }
  else {
    fprintf(stderr, "\t from file \'%s\'", name );
  }
  fprintf(stderr, " using format \'%s\'\n", (im->imageFormat)->realName );
  ImageIO_close( im );
  _freeImage(im);
  *error = ImageIO_READING_HEADER;
  if ( name != nullptr ) free( name );
  return nullptr;

}







CGAL_INLINE_FUNCTION
static void _swapImageData( _image *im )
{
  unsigned char *ptr1, *ptr2, b[8];
  unsigned short int si, *ptr3, *ptr4;
  unsigned int        i, *ptr5, *ptr6;
  std::size_t size, length;

  if( _getEndianness() != im->endianness) {

    size = std::size_t(im->xdim) * im->ydim * im->zdim * im->vdim * im->wdim;
    if ( size <= 0 ) return;

    length = size / im->wdim;
    ptr1 = ptr2 = (unsigned char *) im->data;

    /* 2 bytes swap */
    if(im->wdim == 2) {
      /*
        while(length--) {
        b[0] = *ptr1++;
        b[1] = *ptr1++;
        *ptr2++ = b[1];
        *ptr2++ = b[0];
        }
      */
      ptr3 = ptr4 = (unsigned short int *) im->data;
      while( length-- ) {
        si = *ptr3++;
        *ptr4++ = (unsigned short int)(((si >> 8) & 0xff) | (si << 8));
      }
    }

    /* 4 bytes swap */
    else if(im->wdim == 4) {
      /*
        while(length--) {
        b[0] = *ptr1++;
        b[1] = *ptr1++;
        b[2] = *ptr1++;
        b[3] = *ptr1++;
        *ptr2++ = b[3];
        *ptr2++ = b[2];
        *ptr2++ = b[1];
        *ptr2++ = b[0];
        }
      */
      ptr5 = ptr6 = (unsigned int *) im->data;
      while( length-- ) {
        i = *ptr5++;
        *ptr6++ =  (i << 24) | ((i & 0xff00) << 8) | ((i >> 8) & 0xff00) | ((i >> 24) & 0xff);
      }
    }
    /* 8 bytes swap */
    else if(im->wdim == 8) {
      while(length--) {
        b[0] = *ptr1++;
        b[1] = *ptr1++;
        b[2] = *ptr1++;
        b[3] = *ptr1++;
        b[4] = *ptr1++;
        b[5] = *ptr1++;
        b[6] = *ptr1++;
        b[7] = *ptr1++;
        *ptr2++ = b[7];
        *ptr2++ = b[6];
        *ptr2++ = b[5];
        *ptr2++ = b[4];
        *ptr2++ = b[3];
        *ptr2++ = b[2];
        *ptr2++ = b[1];
        *ptr2++ = b[0];
      }
    }
  }
}






/* Read data of an inrimage.
   If im->data is not nullptr, assume that the buffer was previously allocated
   Swap bytes depending on the endianness and the current architecture  */
CGAL_INLINE_FUNCTION
int _readImageData(_image *im) {
  unsigned long size, nread;

  if(im->openMode != OM_CLOSE) {
    size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;

    if ( size <= 0 ) return -3;

    if(!im->data) {
      im->data = (unsigned char *) ImageIO_alloc(size);
      if(!im->data) return -2;
    }

    nread = ImageIO_read(im, im->data, size);
    if(nread != size) return -1;


    /* architecture is big endian and data little endian
       length = nb of points
     */
    _swapImageData( im );


  }

  return 1;
}





/* Read data of a vectorial inrimage, making the resulting buffer non-
   inerlaced.
   If im->data is not nullptr, assume that the buffer was previously allocated
   Swap bytes depending on the endianness and the current architecture. */
CGAL_INLINE_FUNCTION
int _readNonInterlacedImageData(_image *im) {
  unsigned long size, nread;
  unsigned char **vp, *buf;
  unsigned int i, j, k, v, w;

  if(im->vdim == 1) return _readImageData(im);

  if(im->openMode != OM_CLOSE) {
    size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;

    if ( size <= 0 ) return -3;

    if(!im->data) {
      im->data = (unsigned char *) ImageIO_alloc(size);
      if(!im->data) return -2;
    }

    vp = (unsigned char **) ImageIO_alloc(im->vdim * sizeof(unsigned char *));
    buf = (unsigned char *) ImageIO_alloc(im->vdim * im->wdim);
    size = im->xdim * im->ydim * im->zdim * im->wdim;
    for(v = 0; v < im->vdim; v++)
      vp[v] = (unsigned char *) im->data + v * size;

    for(k = 0; k < im->zdim; k++) {
      for(j = 0; j < im->ydim; j++) {
        for(i = 0; i < im->xdim; i++) {
          nread = ImageIO_read(im, buf, im->vdim * im->wdim);
          if(nread != im->vdim * im->wdim) return -1;
          for(v = 0; v < im->vdim; v++)
            for(w = 0; w < im->wdim; w++)
              *vp[v]++ = *buf++;
          buf -= im->vdim * im->wdim;
        }
      }
    }

    ImageIO_free(buf);
    ImageIO_free(vp);

    /* architecture is big endian and data little endian */
    _swapImageData( im );


    /* reorder lines */
    /* no non-interlaced data for ANALYZE. But if ever... */
/*     if( im->imageFormat == IF_ANALYZE ) { */
/*        int v ; */
/*        int vdim = im->vdim ; */
/*        int lineSize = im->wdim * im->xdim ; */
/*        int vsize = lineSize * im->ydim * im->zdim ; */
/*        char* swapped = ImageIO_alloc(lineSize) ; */
/*        for( v = 0 ; v < vdim ; ++v ) */
/*        { */
/*           char* buf1 = (char*)im->data + v*vsize ; */
/*           char* buf2 = buf1 + vsize - lineSize ; */

/*           while( buf1 < buf2 ) */
/*           { */
/*              memcpy( swapped, buf1, lineSize ) ; */
/*              memcpy( buf1, buf2, lineSize ) ; */
/*              memcpy( buf2, swapped, lineSize ) ; */
/*              buf1 += lineSize ; */
/*              buf2 -= lineSize ; */
/*           } */

/*           ImageIO_free( swapped ) ; */
/*        } */
/*     } */
  }

  return 1;
}


/* Reads body from a non-interlaced vectorial inrimage whose header has
   been read by _readImageHeader. The image buffer is interlaced. */
CGAL_INLINE_FUNCTION
int _readNonInterlacedFileData(_image *im) {
  unsigned long size, nread;
  unsigned char *ptr1, *vp, *buf;
  unsigned int i, j, k, v, w;

  if(im->vdim == 1) return _readImageData(im);

  if(im->openMode != OM_CLOSE) {
    size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;

    if ( size <= 0 ) return -3;

    if(!im->data) {
      im->data = (unsigned char *) ImageIO_alloc(size);
      if(!im->data) return -2;
    }

    size = im->xdim * im->ydim * im->zdim * im->wdim;
    buf = ptr1 = (unsigned char *) ImageIO_alloc(size);

    for(v = 0; v < im->vdim; v++) {
      buf = ptr1;
      nread = ImageIO_read(im, buf, size);
      if(nread != size) return -1;
      vp = (unsigned char *) im->data + (v * im->wdim);
      for(k = 0; k < im->zdim; k++) {
        for(j = 0; j < im->ydim; j++) {
          for(i = 0; i < im->xdim; i++) {
            for(w = 0; w < im->wdim; w++) *vp++ = *buf++;
            vp += (im->vdim - 1) * im->wdim;
          }
        }
      }
    }

    ImageIO_free(buf);

    /* architecture is big endian and data little endian */
    _swapImageData( im );


    /* reorder lines */
    /* no non-interlaced data for ANALYZE. But if ever... */
/*     if( im->imageFormat == IF_ANALYZE ) { */
/*        int v ; */
/*        int vdim = im->vdim ; */
/*        int lineSize = im->wdim * im->xdim ; */
/*        int vsize = lineSize * im->ydim * im->zdim ; */
/*        char* swapped = ImageIO_alloc(lineSize) ; */
/*        for( v = 0 ; v < vdim ; ++v ) */
/*        { */
/*           char* buf1 = (char*)im->data + v*vsize ; */
/*           char* buf2 = buf1 + vsize - lineSize ; */

/*           while( buf1 < buf2 ) */
/*           { */
/*              memcpy( swapped, buf1, lineSize ) ; */
/*              memcpy( buf1, buf2, lineSize ) ; */
/*              memcpy( buf2, swapped, lineSize ) ; */
/*              buf1 += lineSize ; */
/*              buf2 -= lineSize ; */
/*           } */

/*           ImageIO_free( swapped ) ; */
/*        } */
/*     } */
  }

  return 1;
}










/*--------------------------------------------------
 *
 * ?????
 *
 --------------------------------------------------*/





/* check the type of image in fileName */
CGAL_INLINE_FUNCTION
PTRIMAGE_FORMAT imageType(const char *fileName) {
  _ImageIO_file f;
  char magic[5];
  PTRIMAGE_FORMAT format;

  if(!fileName) {
#ifdef CGAL_USE_ZLIB
    f = gzdopen(fileno(stdin), "rb");
#else
    f = fdopen(fileno(stdin), "rb");
#endif
  }
  else {
#ifdef CGAL_USE_ZLIB
    f = gzopen(fileName, "rb");
#else
    f = fopen(fileName, "rb");
#endif
  }

  if(!f) return nullptr;

#ifdef CGAL_USE_ZLIB
  gzread( f, (void *) magic, 4);
#else
  fread( (void *) magic, 1, 4, f );
#endif


  magic[4] = '\0';

#ifdef CGAL_USE_ZLIB
  gzclose( f );
#else
  if(fileName) fclose( f );
#endif

  if (get_static_firstFormat()==nullptr)
    initSupportedFileFormat();

  for(format=get_static_firstFormat();(format!=nullptr);format=format->next) {
    /* test if it is the correct header based on magic and file extension */
    if (((*format->testImageFormat)(magic,fileName)) >=0) {
      return format;
    }
  }
  return 0;

}





/*--------------------------------------------------
 *
 * Image Format Management
 *
 --------------------------------------------------*/





/** adds a format at the beginning of the list of image formats.
    Test if all mandatory fields have been filled */
CGAL_INLINE_FUNCTION
int addImageFormat( PTRIMAGE_FORMAT format)
{
  if ( (format->testImageFormat) &&
       (format->readImageHeader) &&
       (strlen(format->fileExtension)>0) &&
       (strlen(format->realName)>0) ) {

    format->next=get_static_firstFormat();
    get_static_firstFormat()=format;

    return 0;

  }
  else {
    fprintf(stderr,"addImageFormat: information missing in file format %s\n",
            format->realName);
    return -1;
  }
}

/** adds a format at the end of the list of image formats.
    Test if all mandatory fields have been filled */
CGAL_INLINE_FUNCTION
int addImageFormatAtEnd( PTRIMAGE_FORMAT format)
{
  PTRIMAGE_FORMAT f;
  if ( (format->testImageFormat) &&
       (format->readImageHeader) &&
       (strlen(format->fileExtension)>0) &&
       (strlen(format->realName)>0) ) {

    format->next = nullptr;

    if (get_static_firstFormat() == nullptr) {
      get_static_firstFormat()=format;
    }
    else {
      for(f=get_static_firstFormat();(f->next!=nullptr);f=f->next)
        ;
      f->next=format;
    }

    return 0;

  }
  else {
    fprintf(stderr,"addImageFormatAtEnd: information missing in file format %s\n",
            format->realName);
    return -1;
  }
}


/** creates supported image formats */
CGAL_INLINE_FUNCTION
void initSupportedFileFormat()
{
  PTRIMAGE_FORMAT f;
  if ( get_static_inrimageFormat() == nullptr ) {
    f = createAnalyzeFormat();
    addImageFormatAtEnd( f );
    f = createBMPFormat();
    addImageFormatAtEnd( f );
    f = createGifFormat();
    addImageFormatAtEnd( f );
    f = createGisFormat();
    addImageFormatAtEnd( f );
    f = createIrisFormat();
    addImageFormatAtEnd( f );
    f = createPgmFormat();
    addImageFormatAtEnd( f );
    f = createPgmAscIIFormat();
    addImageFormatAtEnd( f );
    f = createPpmFormat();
    addImageFormatAtEnd( f );
    get_static_inrimageFormat() = createInrimageFormat();
    addImageFormat( get_static_inrimageFormat() );
  }
}



CGAL_INLINE_FUNCTION
PTRIMAGE_FORMAT firstImageFormat() {
  return get_static_firstFormat();
}



/** prints supported image formats */
CGAL_INLINE_FUNCTION
void printSupportedFileFormat() {
  PTRIMAGE_FORMAT f;
  int i;

  initSupportedFileFormat();

  for(i=0, f=get_static_firstFormat();(f!=nullptr);i++, f=f->next) {
    if ( (f->testImageFormat) &&
         (f->readImageHeader) &&
         (strlen(f->fileExtension)>0) &&
         (strlen(f->realName)>0)) {
      fprintf( stderr, "#%2d: format name ='%s', extensions='%s'",
              i, f->realName, f->fileExtension );
      if (f->readImageHeader)
        fprintf( stderr, ", read" );
      if (f->writeImage)
        fprintf( stderr, ", write" );
      fprintf( stderr, "\n" );
    }
  }
}


/** remove supported image formats */
CGAL_INLINE_FUNCTION
void removeSupportedFileFormat() {
  PTRIMAGE_FORMAT f=get_static_firstFormat();

  while( f != nullptr) {
    PTRIMAGE_FORMAT f_old = f;
    f = f->next;
    ImageIO_free( f_old);
  }
  get_static_inrimageFormat()=nullptr;

}


/** trilinear interpolation in an _image float type
 */
CGAL_INLINE_FUNCTION
float triLinInterp(const _image* image,
                   float posx,
                   float posy,
                   float posz,
                   float value_outside /*= 0.f */)
{
  const std::size_t dimx = image->xdim;
  const std::size_t dimy = image->ydim;
  const std::size_t dimz = image->zdim;
  const std::size_t dimxy = dimx*dimy;

  if(posx < 0.f || posy < 0.f || posz < 0.f )
    return value_outside;

  posx = static_cast<float>(posx /(image->vx));
  posy = static_cast<float>(posy /(image->vy));
  posz = static_cast<float>(posz /(image->vz));

  //patch suggested by J.Cugnoni to prevent integer overflow
  if(posz >= float(dimz-1) || posy >= float(dimy-1) || posx >= float(dimx-1))
    return value_outside;

  const int i1 = (int)(posz);
  const int j1 = (int)(posy);
  const int k1 = (int)(posx);

  const int i2 = i1 + 1;
  const int j2 = j1 + 1;
  const int k2 = k1 + 1;

  const float KI2 = float(i2)-posz;
  const float KI1 = posz-float(i1);
  const float KJ2 = float(j2)-posy;
  const float KJ1 = posy-float(j1);

  CGAL_IMAGE_IO_CASE
    (image,
     Word *array = (Word *) image->data;
     return (((float)array[i1 * dimxy + j1 * dimx + k1] * KI2 +
              (float)array[i2 * dimxy + j1 * dimx + k1] * KI1) * KJ2 +
             ((float)array[i1 * dimxy + j2 * dimx + k1] * KI2 +
              (float)array[i2 * dimxy + j2 * dimx + k1] * KI1) * KJ1) * (float(k2)-posx)+
     (((float)array[i1 * dimxy + j1 * dimx + k2] * KI2 +
       (float)array[i2 * dimxy + j1 * dimx + k2] * KI1) * KJ2 +
      ((float)array[i1 * dimxy + j2 * dimx + k2] * KI2 +
       (float)array[i2 * dimxy + j2 * dimx + k2] * KI1) * KJ1) * (posx-float(k1));
     );
  return 0.f;
}

// Gives the value of the image at pixel (i,j,k), converted in float.
CGAL_INLINE_FUNCTION
float evaluate(const _image* image,
               const std::size_t i,
               const std::size_t j,
               const std::size_t k)
{
  using CGAL::IMAGEIO::static_evaluate;

  CGAL_IMAGE_IO_CASE(image, return (float)static_evaluate<Word>(image, i, j, k); );

  return 0.f;
}

/** convert the data of the image to float
*/
CGAL_INLINE_FUNCTION
void convertImageTypeToFloat(_image* image){
  if(image->wordKind == WK_FLOAT && image->wdim == 4)
    return;

  const std::size_t dimx = image->xdim;
  const std::size_t dimy = image->ydim;
  const std::size_t dimz = image->zdim;

  float * array = (float*)ImageIO_alloc (dimx * dimy * dimz *sizeof(float));
  if (array == nullptr ) {
    fprintf ( stderr, "allocation error\n" );
    return;
  }

  CGAL_IMAGE_IO_CASE
    (image,
     Word * typedArray  = (Word *)(image->data);
     for(unsigned int i = 0; i<dimx * dimy * dimz;++i)
       array[i] = (float)(typedArray[i]);
     );

  ImageIO_free ( image->data );
  image->data = array;

  image->wordKind = WK_FLOAT;
  image->wdim = 4;
}

