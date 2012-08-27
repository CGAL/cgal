// Copyright (c) 2005-2008 ASCLEPIOS Project, INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of the ImageIO Library, and as been adapted for
// CGAL (www.cgal.org).
// You can redistribute it and/or  modify it under the terms of the
// GNU Lesser General Public License as published by the Free Software Foundation;
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// These files are provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
//
//
// Author(s)     :  ASCLEPIOS Project (INRIA Sophia-Antipolis), Laurent Rineau

/*
 * xviris.c - load routine for IRIS 'rgb' format pictures
 *
 * LoadIRIS()
 * WriteIRIS()
 */
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "iris.h"



/** Magic header for RGB files */
#define IRIS_MAGIC 0732



typedef unsigned char byte;

#define BPPMASK			0x00ff
#define ITYPE_VERBATIM		0x0000
#define ITYPE_RLE		0x0100
#define ISRLE(type)		(((type) & 0xff00) == ITYPE_RLE)
#define ISVERBATIM(type)	(((type) & 0xff00) == ITYPE_VERBATIM)
#define BPP(type)		((type) & BPPMASK)
#define RLE(bpp)		(ITYPE_RLE | (bpp))
#define VERBATIM(bpp)		(ITYPE_VERBATIM | (bpp))






#define TAGLEN	(5)

#define RINTLUM (79)
#define GINTLUM (156)
#define BINTLUM (21)

#define OFFSET_R	3	/* this is byte order dependent */
#define OFFSET_G	2
#define OFFSET_B	1
#define OFFSET_A	0

#define ILUM(r,g,b)     ((int)(RINTLUM*(r)+GINTLUM*(g)+BINTLUM*(b))>>8)
#define CHANOFFSET(z)	(3-(z))	/* this is byte order dependent */

static byte    *getimagedata  (const _image *im, unsigned short, int, int, int);
static void     interleaverow (byte *, byte *, int, int);
static void     expandrow     (byte *, byte *, int);
static void     readtab       (const _image *im, unsigned long *, int);
static void     addimgtag     (byte *, int, int);


static unsigned short  getshort      (const _image *im);
static unsigned long   getlong       (const _image *im);



int testIrisHeader(char *magic,const char *) {
  if((((unsigned char *)magic)[0]<<8) + ((unsigned char *)magic)[1]
     == IRIS_MAGIC)
    return 0;
  else 
    return -1;
}

PTRIMAGE_FORMAT createIrisFormat() {
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));

  f->testImageFormat=&testIrisHeader;
  f->readImageHeader=&readIrisImage;
  f->writeImage=0;
  strcpy(f->fileExtension,".rgb");
  strcpy(f->realName,"IRIS");
  return f;
}



/******************************************/
static void interleaverow(byte *lptr, byte *cptr, int z, int n)
{
  lptr += z;
  while(n--) {
    *lptr = *cptr++;
    lptr += 4;
  }
}


/******************************************/
static void expandrow(byte *optr, byte *iptr, int z)
{
  byte pixel, count;

  optr += z;
  while (1) {
    pixel = *iptr++;
    if ( !(count = (pixel & 0x7f)) ) return;
    if (pixel & 0x80) {
      while (count>=8) {
	optr[0*4] = iptr[0];
	optr[1*4] = iptr[1];
	optr[2*4] = iptr[2];
	optr[3*4] = iptr[3];
	optr[4*4] = iptr[4];
	optr[5*4] = iptr[5];
	optr[6*4] = iptr[6];
	optr[7*4] = iptr[7];
	optr += 8*4;
	iptr += 8;
	count -= 8;
      }
      while(count--) {
	*optr = *iptr++;
	optr+=4;
      }
    }
    else {
      pixel = *iptr++;
      while(count>=8) {
	optr[0*4] = pixel;
	optr[1*4] = pixel;
	optr[2*4] = pixel;
	optr[3*4] = pixel;
	optr[4*4] = pixel;
	optr[5*4] = pixel;
	optr[6*4] = pixel;
	optr[7*4] = pixel;
	optr += 8*4;
	count -= 8;
      }
      while(count--) {
	*optr = pixel;
	optr+=4;
      }
    }
  }
}

/****************************************************/
static void readtab(const _image *im, unsigned long *tab, int n)
{
  while (n) {
    *tab++ = getlong(im);
    n--;
  }
}


/*****************************************************/
static void addimgtag(byte *dptr, int xsize, int ysize)
{
  /* this is used to extract image data from core dumps. 
     I doubt this is necessary...  --jhb */

  dptr    = dptr + (xsize * ysize * 4);
  dptr[0] = 0x12;  dptr[1] = 0x34;  dptr[2] = 0x56;  dptr[3] = 0x78;
  dptr += 4;

  dptr[0] = 0x59;  dptr[1] = 0x49;  dptr[2] = 0x33;  dptr[3] = 0x33;
  dptr += 4;

  dptr[0] = 0x69;  dptr[1] = 0x43;  dptr[2] = 0x42;  dptr[3] = 0x22;
  dptr += 4;

  dptr[0] = (xsize>>24)&0xff;
  dptr[1] = (xsize>>16)&0xff;
  dptr[2] = (xsize>> 8)&0xff;
  dptr[3] = (xsize    )&0xff;
  dptr += 4;

  dptr[0] = (ysize>>24)&0xff;
  dptr[1] = (ysize>>16)&0xff;
  dptr[2] = (ysize>> 8)&0xff;
  dptr[3] = (ysize    )&0xff;
}

/* byte order independent read/write of shorts and longs. */
/*****************************************************/
static unsigned short getshort( const _image *im)
{
  byte buf[2];
  ImageIO_read( im, buf, (size_t) 2);
  return (buf[0]<<8)+(buf[1]<<0);
}

/*****************************************************/
static unsigned long getlong( const _image *im )
{
  byte buf[4];
  ImageIO_read( im, buf, (size_t) 4);
  return (((unsigned long) buf[0])<<24) + (((unsigned long) buf[1])<<16)
       + (((unsigned long) buf[2])<<8) + buf[3];
}

/*****************************************************/
int readIrisImage( const char *, _image *im ) {
  byte   *rawdata, *rptr;
  byte   *pic824,  *bptr, *iptr;
  int     i, j, size;
  unsigned short imagic, type;
  int xsize, ysize, zsize;


  /* read header information from file */
  imagic = getshort( im );
  type   = getshort( im );
  getshort( im );
  xsize  = getshort( im );
  ysize  = getshort( im );
  zsize  = getshort( im );

  if ( ImageIO_error(im) ) return( 0 );

  if (imagic != IRIS_MAGIC) return( 0 );

  rawdata = getimagedata(im, type, xsize, ysize, zsize);
  if (!rawdata) return( 0 );

  if ( ImageIO_error(im) ) return( 0 );   /* probably truncated file */


  /*  t1=texture_alloc();
  (void) strcpy(t1->name,fname);
  t1->type=IRIS;
  
  t1->nrows = ysize;
  t1->ncols = xsize;
    
  t1->image = create_int_array(t1->nrows,t1->ncols);*/

  /* got the raw image data.  Convert to an XV image (1,3 bytes / pix) */
  if (zsize < 3) 
    {  /* grayscale */
      im->xdim = xsize;
      im->ydim = ysize;
      im->zdim = 1;
      im->vdim = 1;
      im->wdim = 1;
      im->wordKind = WK_FIXED;
      im->sign = SGN_UNSIGNED;
      im->data = ImageIO_alloc(xsize * ysize);

      pic824 = (byte *) ImageIO_alloc((size_t) xsize * ysize);
      if (!pic824) exit(-1);
      
      /* copy plane 3 from rawdata into pic824, inverting pic vertically */
      for (i = 0, bptr = pic824; i < (int) ysize; i++) 
	{
	  rptr = rawdata + 3 + ((ysize - 1) - i) * (xsize * 4);
	  for (j = 0; j < (int) xsize; j++, bptr++, rptr += 4)
	    *bptr = *rptr;
	}
      
      
      size = im->xdim * im->ydim;
      for (bptr = pic824, iptr = (unsigned char *) im->data,
	     i = 0; i < size; ++i, ++iptr, ++bptr)
	{
	  *iptr = *bptr;
	}
      
      ImageIO_free(pic824);
    
  }

  else {  /* truecolor */
    im->xdim = xsize;
    im->ydim = ysize;
    im->zdim = zsize / 3;
    im->vdim = 4;
    im->wdim = 1;
    im->wordKind = WK_FIXED;
    im->sign = SGN_UNSIGNED;
    im->data = ImageIO_alloc(xsize * ysize * im->zdim * 4);

    pic824 = (byte *) ImageIO_alloc((size_t) xsize * ysize * 3);
    if (!pic824) 
      exit(1);
    
    /* copy plane 3 from rawdata into pic824, inverting pic vertically */
    for (i = 0, bptr = pic824; i<(int) ysize; i++) 
      {
	rptr = rawdata + ((ysize - 1) - i) * (xsize * 4);

	for (j=0; j<(int) xsize; j++, rptr += 4) {
	  *bptr++ = rptr[3];
	  *bptr++ = rptr[2];
	  *bptr++ = rptr[1];
	}
      }
    
    
    size = im->xdim * im->ydim;
    for (bptr = pic824, iptr = (unsigned char *) im->data, i = 0;
	 i < size; ++i, iptr += 4, bptr += 3)
      {
	if ( _getEndianness() == END_LITTLE ) {
	  iptr[0] = 0xFF;
	  iptr[1] = bptr[2];
	  iptr[2] = bptr[1];
	  iptr[3] = bptr[0];
	}
	else {
	  iptr[0] = bptr[0];
	  iptr[1] = bptr[1];
	  iptr[2] = bptr[2];
	  iptr[3] = 0xFF;
	}
      }
      ImageIO_free(pic824);
  }

  ImageIO_free(rawdata);

  return 1;
}     




/****************************************************/
static byte *getimagedata(const _image *im, unsigned short type, int xsize, int ysize,
			  int zsize)
{
  /* read in a B/W RGB or RGBA iris image file and return a 
     pointer to an array of 4-byte pixels, arranged ABGR, NULL on error */

  byte   *base, *lptr;
  byte   *verdat;
  int     y, z, tablen;
  int     bpp, rle, badorder;
  int     rlebuflen;
  byte *rledat;
  unsigned long *starttab, *lengthtab, cur;


  rle     = ISRLE(type);
  bpp     = BPP(type);

  if (bpp != 1) {
    return (byte *) NULL;
  }

  if (rle) {

    rlebuflen = 2 * xsize + 10;
    tablen    = ysize * zsize;
    starttab  = (unsigned long *) ImageIO_alloc((size_t) tablen * sizeof(long));
    lengthtab = (unsigned long *) ImageIO_alloc((size_t) tablen * sizeof(long));
    rledat    = (byte *)   ImageIO_alloc((size_t) rlebuflen);

    if (!starttab || !lengthtab || !rledat) 
      exit(1);

    ImageIO_seek( im, 512L, SEEK_SET );
    readtab(im, starttab,  tablen);
    readtab(im, lengthtab, tablen);

    if ( ImageIO_error(im) ) {
      ImageIO_free(starttab);  ImageIO_free(lengthtab);  ImageIO_free(rledat);
      return (byte *) NULL;
    }


    /* check data order */
    cur = 0;
    badorder = 0;
    for (y=0; y<ysize && !badorder; y++) {
      for (z=0; z<zsize && !badorder; z++) {
	if (starttab[y+z*ysize] < cur) badorder = 1;
	else cur = starttab[y+z*ysize];
      }
    }

    ImageIO_seek( im, (long) (512 + 2*tablen*4), 0);

    cur = 512 + 2*tablen*4;

    base = (byte *) ImageIO_alloc((size_t) (xsize*ysize+TAGLEN) * 4);
    if (!base) 
      exit(1);

    addimgtag(base,xsize,ysize);

    if (badorder) {
      for (z=0; z<zsize; z++) {
	lptr = base;
	for (y=0; y<ysize; y++) {
	  if (cur != starttab[y+z*ysize]) {
	    ImageIO_seek( im, (long) starttab[y+z*ysize], 0);
	    cur = starttab[y+z*ysize];
	  }

	  if (lengthtab[y+z*ysize]>(unsigned long)rlebuflen) {
	    ImageIO_free(starttab); ImageIO_free(lengthtab); ImageIO_free(rledat); ImageIO_free(base);
	    return (byte *) NULL;
	  }
	  ImageIO_read(im, rledat, (size_t) lengthtab[y+z*ysize]);
	  cur += lengthtab[y+z*ysize];
	  expandrow(lptr,rledat,3-z);
	  lptr += (xsize * 4);
	}
      }
    }
    else {
      lptr = base;
      for (y=0; y<ysize; y++) {
	for (z=0; z<zsize; z++) {
	  if (cur != starttab[y+z*ysize]) {
	    ImageIO_seek( im, (long) starttab[y+z*ysize], 0);
	    cur = starttab[y+z*ysize];
	  }
	  
	  ImageIO_read(im, rledat, (size_t) lengthtab[y+z*ysize]);
	  cur += lengthtab[y+z*ysize];
	  expandrow(lptr,rledat,3-z);
	}
	lptr += (xsize * 4);
      }
    }

    ImageIO_free(starttab);
    ImageIO_free(lengthtab);
    ImageIO_free(rledat);
    return base;
  }      /* end of RLE case */

  else {  /* not RLE */
    verdat = (byte *) ImageIO_alloc((size_t) xsize);
    base   = (byte *) ImageIO_alloc((size_t) (xsize*ysize+TAGLEN) * 4);
    if (!base || !verdat) 
      exit(1);

    addimgtag(base,xsize,ysize);

    ImageIO_seek( im, 512L,0);

    for (z=0; z<zsize; z++) {
      lptr = base;
      for (y=0; y<ysize; y++) {
	ImageIO_read(im, verdat, (size_t) xsize);
	interleaverow(lptr,verdat,3-z,xsize);
	lptr += (xsize * 4);
      }
    }

    ImageIO_free(verdat);
    return base;
  }
}
















