// Copyright (c) 2005-2008 ASCLEPIOS Project, INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of the ImageIO Library, and as been adapted for CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
//
//
// Author(s)     :  ASCLEPIOS Project (INRIA Sophia-Antipolis), Laurent Rineau

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <stdio.h>
#include <stdlib.h>
#include <CGAL/ImageIO/bmptypes.h>
#include <CGAL/ImageIO/bmpendian.h>
#include <CGAL/ImageIO/bmpread.h>






/** Magic header for bmp */
#define BMP_MAGIC "BM"


#ifdef CGAL_HEADER_ONLY

inline int& get_static_verbose_bmp()
{
  static int _VERBOSE_BMP_ = 1;
  return _VERBOSE_BMP_;
}

#else // CGAL_HEADER_ONLY

static int _VERBOSE_BMP_ = 1;
inline int& get_static_verbose_bmp()
{ return _VERBOSE_BMP_; }

#endif // CGAL_HEADER_ONLY


CGAL_INLINE_FUNCTION
int testBmpHeader(char *magic,const char *) {
  if( !strncmp(magic, BMP_MAGIC, strlen(BMP_MAGIC)))
    return 0;
  else
    return -1;
}

CGAL_INLINE_FUNCTION
PTRIMAGE_FORMAT createBMPFormat() {
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));

  f->testImageFormat=&testBmpHeader;
  f->readImageHeader=&readBmpImage;
  f->writeImage=0;
  strcpy(f->fileExtension,".bmp");
  strcpy(f->realName,"BMP");
  return f;
}

CGAL_INLINE_FUNCTION
int readBmpImage( const char *name,_image *im)
{
  int dimx, dimy, dimv;

  im->data = _readBmpImage( name, &dimx,  &dimy,  &dimv );
  if ( im->data == nullptr ) {
    fprintf( stderr, "readBmpImage: unable to read \'%s\'\n", name );
    return( -1 );
  }

  im->xdim = dimx;
  im->ydim = dimy;
  im->vdim = dimv;

  im->zdim = 1;

  im->wdim = 1;
  im->wordKind = WK_FIXED;
  im->sign = SGN_UNSIGNED;

  return 1;
}



CGAL_INLINE_FUNCTION
void *_readBmpImage( const char *name,
                     int *dimx, int *dimy, int *dimz )
{
  const char *proc="_readBmpImage";
  void *buf = (void*)nullptr;
  unsigned char *myBuf = nullptr;

  FILE *fp = nullptr;
  RGB **argbs = nullptr;
  char **xorMasks = nullptr, **andMasks = nullptr;
  CGAL_UINT32 *heights = nullptr, *widths = nullptr, row = 0, col = 0;
  CGAL_UINT16 fileType = 0;
  long filePos = 0;
  int numImages = 0, i = 0;
  int rc = 0;

    fp = fopen(name, "rb");
    if (fp == nullptr) {
      if ( get_static_verbose_bmp() )
        fprintf( stderr, "%s: error in opening %s\n", proc, name );
      return( (void*)nullptr );
    }


    /*
     * Read the first two bytes as little-endian to determine the file type.
     * Preserve the file position.
     */
    filePos = ftell(fp);
    rc = readUINT16little(fp, &fileType);
    if (rc != 0) {
      fclose(fp);
      if ( get_static_verbose_bmp() )
        fprintf( stderr, "%s: error in getting file type %s\n", proc, name );
      return( (void*)nullptr );
    }
    fseek(fp, filePos, SEEK_SET);

    /*
     * Read the images.
     */
    switch (fileType) {
    case TYPE_ARRAY:
        /*
         * If this is an array of images, read them.  All the arrays we need
         * will be allocated by the reader function.
         */
        rc = readMultipleImage(fp, &argbs, &xorMasks, &andMasks, &heights,
                               &widths, &numImages);
        break;
    case TYPE_BMP:
    case TYPE_ICO:
    case TYPE_ICO_COLOR:
    case TYPE_PTR:
    case TYPE_PTR_COLOR:
        /*
         * If this is a single-image file, we've a little more work.  In order
         * to make the output part of this test program easy to write, we're
         * going to allocate dummy arrays that represent what
         * readMultipleImage would have allocated.  We'll read the data into
         * those arrays.
         */
        argbs = (RGB **)calloc(1, sizeof(RGB *));
        if (argbs == nullptr)
        {
            rc = 1005;
            break;
        }
        xorMasks = (char **)calloc(1, sizeof(char *));
        if (xorMasks == nullptr)
        {
            free(argbs);
            rc = 1005;
            break;
        }
        andMasks = (char **)calloc(1, sizeof(char *));
        if (andMasks == nullptr)
        {
            free(argbs);
            free(xorMasks);
            rc = 1005;
            break;
        }
        heights = (CGAL_UINT32 *)calloc(1, sizeof(CGAL_UINT32));
        if (heights == nullptr)
        {
            free(argbs);
            free(xorMasks);
            free(andMasks);
            rc = 1005;
            break;
        }
        widths = (CGAL_UINT32 *)calloc(1, sizeof(CGAL_UINT32));
        if (widths == nullptr)
        {
            free(argbs);
            free(xorMasks);
            free(andMasks);
            free(heights);
            rc = 1005;
            break;
        }
        numImages = 1;

        /*
         * Now that we have our arrays allocted, read the image into them.
         */
        switch (fileType) {
        case TYPE_BMP:
            rc = readSingleImageBMP(fp, argbs, widths, heights);
            break;
        case TYPE_ICO:
        case TYPE_PTR:
            rc = readSingleImageICOPTR(fp, xorMasks, andMasks, widths,
                                       heights);
            break;
        case TYPE_ICO_COLOR:
        case TYPE_PTR_COLOR:
            rc = readSingleImageColorICOPTR(fp, argbs, xorMasks, andMasks,
                                            widths, heights);
            break;
        }
        break;
    default:
        rc = 1000;
    }
    fclose(fp);



    /*
     * At this point, everything's been read.  Display status messages based
     * on the return values.
     */
    switch (rc) {
    case 1000:
    case 1006:
      if ( get_static_verbose_bmp() )
        fprintf( stderr, "%s: File is not a valid bitmap file\n", proc );
      break;
    case 1001:
      if ( get_static_verbose_bmp() )
        fprintf( stderr, "%s: Illegal information in an image\n", proc );
      break;
    case 1002:
      if ( get_static_verbose_bmp() )
        fprintf( stderr, "%s: Legal information that I can't handle yet in an image\n", proc );
      break;
    case 1003:
    case 1004:
    case 1005:
      if ( get_static_verbose_bmp() )
        fprintf( stderr, "%s: Ran out of memory\n", proc );
      break;
    case 0:
      if ( get_static_verbose_bmp() > 1 )
        fprintf( stderr, "%s: Got good data from file, writing results\n", proc );
      break;
    default:
      if ( get_static_verbose_bmp() )
        fprintf( stderr, "%s: Error reading file rc=%d\n", proc,rc );
      break;
    }

    /*
     * If the return value wasn't 0, something went wrong.
     */
    if (rc != 0)
    {
        if (rc != 1000 && rc != 1005)
        {
            for (i=0; i<numImages; i++)
            {
                if (argbs[i] != nullptr)
                    free(argbs[i]);
                if (andMasks[i] != nullptr)
                    free(andMasks[i]);
                if (xorMasks[i] != nullptr)
                    free(xorMasks[i]);
            }
            free(argbs);
            free(andMasks);
            free(xorMasks);
            free(widths);
            free(heights);
        }
        return( (void*)nullptr );
    }


    /*
     * Dump the images.
     */
    if ( get_static_verbose_bmp() > 1 )
      fprintf (stderr, "%s: There are %d images in the file\n", proc, numImages);

    if ( numImages >= 2 )
      fprintf (stderr, "%s: read only first image among %d\n", proc, numImages );

    /*
     * my stuff:
     * just reading one bmp image
     */
    if ( (numImages > 0) &&
         (argbs[0] != nullptr) ) {

      buf = (void*)malloc( widths[0]*heights[0]*3 * sizeof( unsigned char ) );
      if ( buf == (void*)nullptr ) {
        if ( get_static_verbose_bmp() )
          fprintf( stderr, "%s: error in allocating data buffer for %s\n", proc, name );

        for (i=0; i<numImages; i++) {
          if (argbs[i] != nullptr)
            free(argbs[i]);
          if (andMasks[i] != nullptr)
            free(andMasks[i]);
          if (xorMasks[i] != nullptr)
            free(xorMasks[i]);
        }
        free(argbs);
        free(andMasks);
        free(xorMasks);
        free(widths);
        free(heights);

        return( (void*)nullptr );
      }

      myBuf = (unsigned char*)buf;
      i = 0;
      for (row = 0; row < heights[0]; row++)
      for (col = 0; col < widths[0]; col++,i+=3) {
        myBuf[i  ] = argbs[0][row * widths[0] + col].red;
        myBuf[i+1] = argbs[0][row * widths[0] + col].green;
        myBuf[i+2] = argbs[0][row * widths[0] + col].blue;
      }

      *dimx = int(widths[0]);
      *dimy = int(heights[0]);
      *dimz = 3;

    } else {
      if ( get_static_verbose_bmp() )
        fprintf( stderr, "%s: no image or null image\n", proc );

      for (i=0; i<numImages; i++) {
        if (argbs[i] != nullptr)
          free(argbs[i]);
        if (andMasks[i] != nullptr)
          free(andMasks[i]);
        if (xorMasks[i] != nullptr)
          free(xorMasks[i]);
      }
      free(argbs);
      free(andMasks);
      free(xorMasks);
      free(widths);
      free(heights);

      return( (void*)nullptr );
    }



    for (i=0; i<numImages; i++) {
      if (argbs[i] != nullptr)
        free(argbs[i]);
      if (andMasks[i] != nullptr)
        free(andMasks[i]);
      if (xorMasks[i] != nullptr)
        free(xorMasks[i]);
    }
    free(argbs);
    free(andMasks);
    free(xorMasks);
    free(widths);
    free(heights);

    return( buf );





#if 0
    /*
     * old stuff from test.c
     */

    for (i=0; i<numImages; i++)
    {
        /*
         * Loop through all the images that were returned.
         */
      if ( get_static_verbose_bmp() ) {
        fprintf (stderr, "%s: Doing image number %d\n\n", proc, i+1);
        fprintf (stderr, "%s: Image dimensions: (%ld,%ld)\n", proc, widths[i], heights[i]);
      }

      if (argbs[i] != nullptr) {
        /*
         * If the image has colors, dump them (BMP, color ICO and color
         * PTR files
         */
        fprintf(stderr, "Colors");
        for (row = 0; row < heights[i]; row++)
          {
            fprintf (stderr, "\n\nRow %ld pixels (R,G,B), hex values:\n",
                     row);
            for (col = 0; col < widths[i]; col++)
              {
                fprintf (stderr, "(%2.2x,%2.2x,%2.2x)",
                         argbs[i][row * widths[i] + col].red,
                         argbs[i][row * widths[i] + col].green,
                         argbs[i][row * widths[i] + col].blue);
              }
          }
      } else {
        /*
         * If there image has no colors, say so.  (monochrome ICO and PTR
         * files)
         */
        fprintf (stderr, "No color image\n");
      }



      if (xorMasks[i] != nullptr) {
        /*
         * If the image has an xor mask, dump it.  (ICO and PTR files)
         */
        fprintf (stderr, "\nXOR mask\n");
        for (row = 0; row < heights[i]; row++)
          {
            for (col = 0; col < widths[i]; col++)
              {
                fprintf (stderr, "%c",
                         xorMasks[i][row * widths[i] + col] ? '@' : '.');
              }
            fprintf (stderr, "\n");
          }
      } else {
        /*
         * If the image has no xor mask, say so.  (BMP files).
         */
        fprintf (stderr, "No xor mask\n");
      }



      if (andMasks[i] != nullptr) {
        /*
         * If the image has an and mask, dump it.  (ICO and PTR files)
         */
        fprintf (stderr, "\nAND mask\n");
        for (row = 0; row < heights[i]; row++)
          {
            for (col = 0; col < widths[i]; col++)
              {
                fprintf (stderr, "%c",
                         andMasks[i][row * widths[i] + col] ? '@' : '.');
              }
            fprintf (stderr, "\n");
          }
      } else {
        /*
         * If the image has noand mask, say so.  (BMP files)
         */
        fprintf (stderr, "No and mask\n");
      }


      if (i != numImages-1)
        fprintf (stderr, "\n------------------------------------------\n\n");

    }

    /*
     * Dumping is complete.  Free all the arrays and quit
     */
    for (i=0; i<numImages; i++)
    {
        if (argbs[i] != nullptr)
            free(argbs[i]);
        if (andMasks[i] != nullptr)
            free(andMasks[i]);
        if (xorMasks[i] != nullptr)
            free(xorMasks[i]);
    }
    free(argbs);
    free(andMasks);
    free(xorMasks);
    free(widths);
    free(heights);

    return( buf );
#endif
}

CGAL_INLINE_FUNCTION
void IoBmp_verbose ( )
{
  if ( get_static_verbose_bmp() <= 0 )
    get_static_verbose_bmp() = 1;
  else
    get_static_verbose_bmp() += 1;
}

CGAL_INLINE_FUNCTION
void IoBmp_noverbose ( )
{
  get_static_verbose_bmp() = 0;
}
