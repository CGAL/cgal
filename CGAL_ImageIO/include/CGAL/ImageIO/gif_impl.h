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
// SPDX-License-Identifier: LGPL-3.0+
//
//
// Author(s)     :  ASCLEPIOS Project (INRIA Sophia-Antipolis), Laurent Rineau

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <string.h>

#ifdef WIN32
#include <io.h>
#endif

#ifdef CGAL_HEADER_ONLY
  
  #define CGAL_GLOBAL_STATE_VAR(TYPE, NAME, VALUE)  \
    inline TYPE & get_static_##NAME()               \
    {                                               \
      static TYPE NAME = VALUE;                     \
      return NAME;                                  \
    }
  
#else // CGAL_HEADER_ONLY

  #define CGAL_GLOBAL_STATE_VAR(TYPE, NAME, VALUE)  \
    TYPE NAME;                                      \
    inline TYPE& get_static_##NAME()                \
    {                                               \
      return NAME;                                  \
    }

#endif // CGAL_HEADER_ONLY


/** Magic header for GIF files */
#define GIF_MAGIC "GIF8"



typedef unsigned char byte;
#define TRUE 1
#define FALSE 0

#define NEXTBYTE      (*ptr++)
#define EXTENSION     0x21
#define IMAGESEP      0x2c
#define TRAILER       0x3b
#define INTERLACEMASK 0x40
#define COLORMAPMASK  0x80


#define DEBUG 0

//FILE *fp;
//int   gif89 = 0;
static const char *id87 = "GIF87a";
static const char *id89 = "GIF89a";

static int EGApalette[16][3] = {
  {0,0,0},       {0,0,128},     {0,128,0},     {0,128,128}, 
  {128,0,0},     {128,0,128},   {128,128,0},   {200,200,200},
  {100,100,100}, {100,100,255}, {100,255,100}, {100,255,255},
  {255,100,100}, {255,100,255}, {255,255,100}, {255,255,255} };
  

static int  ReadCode();
static void DoInterlace(byte);
static int GifError(const char *);

CGAL_GLOBAL_STATE_VAR(byte *, Raster, NULL) /* The raster data stream, unblocked */
CGAL_GLOBAL_STATE_VAR(byte *, RawGIF, NULL)
CGAL_GLOBAL_STATE_VAR(byte *, r, NULL)
CGAL_GLOBAL_STATE_VAR(byte *, g, NULL)
CGAL_GLOBAL_STATE_VAR(byte *, b, NULL)     /* The colormap */
CGAL_GLOBAL_STATE_VAR(int, BitOffset, 0)    /* Bit Offset of next code */
CGAL_GLOBAL_STATE_VAR(int, XC, 0)
CGAL_GLOBAL_STATE_VAR(int, YC, 0)	    /* Output X and Y coords of current pixel */
CGAL_GLOBAL_STATE_VAR(int, CodeSize, 0)	    /* Code size, read from GIF header */
CGAL_GLOBAL_STATE_VAR(int, ReadMask, 0)     /* Code AND mask for current code size */
CGAL_GLOBAL_STATE_VAR(int, Pass, 0)	    /* Used by output routine if interlaced pic */
CGAL_GLOBAL_STATE_VAR(int, Width, 0)
CGAL_GLOBAL_STATE_VAR(int, Height, 0)       /* image dimensions */
CGAL_GLOBAL_STATE_VAR(unsigned char *, org, NULL)
CGAL_GLOBAL_STATE_VAR(unsigned char *, buf, NULL)

CGAL_INLINE_FUNCTION
int testGifHeader(char *magic,const char *) {
  if (!strcmp(magic, GIF_MAGIC))
    return 0;
  else 
    return -1;
}
CGAL_INLINE_FUNCTION
PTRIMAGE_FORMAT createGifFormat() {
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));

  f->testImageFormat=&testGifHeader;
  f->readImageHeader=&readGifImage;
  f->writeImage=0;
  strcpy(f->fileExtension,".gif");
  strcpy(f->realName,"Gif");
  return f;
}

CGAL_INLINE_FUNCTION
void use(int) {} // warning killer

/*****************************/
CGAL_INLINE_FUNCTION
int readGifImage(const char *name,_image *im) {
FILE *fp;
int   gif89 = 0;

  byte ch, ch1;
  byte *ptr, *ptr1;
  int i, block;
  int npixels, maxpixels, aspect;
  float normaspect;
  int OutCount = 0,		/* Decompressor output 'stack count' */
    RWidth, RHeight,		/* screen dimensions */
    /*LeftOfs, TopOfs,		 image offset */
    BitsPerPixel,		/* Bits per pixel, read from GIF header */
    ColorMapSize,		/* number of colors */
    Background,		         /* background color */
    InitCodeSize,		/* Starting code size, used during Clear */
    Code,			/* Value returned by ReadCode */
    MaxCode,			/* limiting value for current code size */
    ClearCode,			/* GIF clear code */
    EOFCode,			/* GIF end-of-information code */
    CurCode, OldCode=0, InCode,	/* Decompressor variables */
    FirstFree,			/* First free code, generated per GIF spec */
    FreeCode,			/* Decompressor,next free slot in hash table */
    FinChar=0,			/* Decompressor variable */
    BitMask,			/* AND mask for data size */
    Misc;                       /* miscellaneous bits (interlace, local cmap)*/
  int Interlace, HasColormap;
  /* not used
   */
  /* char header[10]; */

  /* The hash table used by the decompressor */
  int Prefix[4096];
  int Suffix[4096];
  /* An output array used by the decompressor */
  int OutCode[4097];

  /* initialize variables */
  get_static_BitOffset() =
    get_static_XC() =
    get_static_YC() =
    get_static_Pass() =
    OutCount =
    npixels =
    maxpixels = 0;
  get_static_RawGIF() = get_static_Raster() = NULL;
  gif89 = 0;

#ifdef WIN32
 fp = fopen(name,"rb");
#else
  fp = fopen(name,"r");
#endif
  fp = fopen(name,"rb");
  if (!fp) {
    return(GifError("could not open a GIF file"));
  }
  /* find the size of the file */
  fseek(fp, 0L, 2);
  long filesize = ftell(fp);
  fseek(fp, 0L, 0);

 /* the +256's are so we can read truncated GIF files without fear of 
     segmentation violation */
  if (!(ptr = get_static_RawGIF() = (byte *) ImageIO_alloc(filesize+256)))
    return( GifError("not enough memory to read gif file") );
  
  if (!(get_static_Raster() = (byte *) ImageIO_alloc(filesize+256)))    
    return( GifError("not enough memory to read gif file") );
  
  if (fread(ptr, filesize, 1, fp) != 1)
    return( GifError("GIF data read failed") );
   if      (strncmp((char *) ptr, id87, 6)==0) gif89 = 0;
  else if (strncmp((char *) ptr, id89, 6)==0) gif89 = 1;
  else    return( GifError("not a GIF file"));
  
  ptr += 6;
 /* Get variables from the GIF screen descriptor */
  
  ch = NEXTBYTE;
  RWidth = ch + 0x100 * NEXTBYTE;	/* screen dimensions... not used. */
  ch = NEXTBYTE;
  RHeight = ch + 0x100 * NEXTBYTE;
  use(RWidth);
  use(RHeight);
  
  ch = NEXTBYTE;
  HasColormap = ((ch & COLORMAPMASK) ? TRUE : FALSE);
  
  BitsPerPixel = (ch & 7) + 1;
  ColorMapSize = 1 << BitsPerPixel;
  BitMask = ColorMapSize - 1;
  
  Background = NEXTBYTE;		/* background color... not used. */
  use(Background);


  aspect = NEXTBYTE;
  if (aspect) {
    if (!gif89) return(GifError("corrupt GIF file (screen descriptor)"));
    else normaspect = (float) (aspect + 15) / 64.0f;   /* gif89 aspect ratio */
    if (DEBUG) fprintf(stderr,"GIF89 aspect = %f\n", normaspect);
  }


  /* Read in global colormap. */
  if (HasColormap)
    {
      get_static_r() = (byte *) ImageIO_alloc(ColorMapSize * sizeof(byte));
      get_static_g() = (byte *) ImageIO_alloc(ColorMapSize * sizeof(byte));
      get_static_b() = (byte *) ImageIO_alloc(ColorMapSize * sizeof(byte));

      for (i = 0; i < ColorMapSize; i++) {
	get_static_r()[i] = NEXTBYTE;
	get_static_g()[i] = NEXTBYTE;
        get_static_b()[i] = NEXTBYTE;
      }
    }
  else {
    /* no colormap in GIF file */
    /* put std EGA palette (repeated 16 times) into colormap, for lack of
       anything better to do */
    ColorMapSize = 256;
    get_static_r() = (byte *) ImageIO_alloc(256 * sizeof(byte));
    get_static_g() = (byte *) ImageIO_alloc(256 * sizeof(byte));
    get_static_b() = (byte *) ImageIO_alloc(256 * sizeof(byte));
   
    for (i = 0; i < 256; i++) {
      get_static_r()[i] = byte(EGApalette[i&15][0]);
      get_static_g()[i] = byte(EGApalette[i&15][1]);
      get_static_b()[i] = byte(EGApalette[i&15][2]);
    }
  }

  /* possible things at this point are:
   *   an application extension block
   *   a comment extension block
   *   an (optional) graphic control extension block
   *       followed by either an image
   *	   or a plaintext extension
   */
  while (1) {
    block = NEXTBYTE;
   
    if (block == EXTENSION) {  /* parse extension blocks */
      int i, fn, blocksize, aspnum, aspden;
	  
      /* read extension block */
      fn = NEXTBYTE;
	  
      if (DEBUG) fprintf(stderr,"GIF extension type 0x%02x\n", fn);
	  
      if (fn == 'R') {                  /* GIF87 aspect extension */
	blocksize = NEXTBYTE;
	if (blocksize == 2) {
	  aspnum = NEXTBYTE;
	  aspden = NEXTBYTE;
	  if (aspden>0 && aspnum>0) 
	    normaspect = (float) aspnum / (float) aspden;
	  else { normaspect = 1.0;  aspnum = aspden = 1; }
			
	  if (DEBUG) fprintf(stderr,"GIF87 aspect extension: %d:%d = %f\n\n", 
			     aspnum, aspden,normaspect);
	}
	else {
	  for (i=0; i<blocksize; i++) (void)NEXTBYTE;
	}
      }
	  
      else if (fn == 0xFE) {  /* Comment Extension.  just eat it */
	int ch, j, sbsize;
		 
	if (DEBUG) fprintf(stderr,"Comment extension:  ");
	/* read (and ignore) data sub-blocks */
	do {
	  j = 0;  sbsize = NEXTBYTE;
	  while (j<sbsize) {
	    ch = NEXTBYTE;  j++;
	    if (DEBUG) fprintf(stderr,"%c", ch);
	  }
	} while (sbsize);
	if (DEBUG) fprintf(stderr,"\n\n");
      }
	  
      else if (fn == 0x01) {  /* PlainText Extension */
	int j,sbsize,ch;
	int tgLeft, tgTop, tgWidth, tgHeight, cWidth, cHeight, fg, bg;
		 
	/*	SetISTR(ISTR_WARNING, 
		"PlainText extension found in GIF file.  Ignored.");*/
		 
	sbsize   = NEXTBYTE;
	tgLeft   = NEXTBYTE;  tgLeft   += (NEXTBYTE)<<8;
	tgTop    = NEXTBYTE;  tgTop    += (NEXTBYTE)<<8;
	tgWidth  = NEXTBYTE;  tgWidth  += (NEXTBYTE)<<8;
	tgHeight = NEXTBYTE;  tgHeight += (NEXTBYTE)<<8;
	cWidth   = NEXTBYTE;
	cHeight  = NEXTBYTE;
	fg       = NEXTBYTE;
	bg       = NEXTBYTE;
	i=12;
	for ( ; i<sbsize; i++) (void)NEXTBYTE;   /* read rest of first subblock */
		 
	if (DEBUG) fprintf(stderr,
			   "PlainText: tgrid=%d,%d %dx%d  cell=%dx%d  col=%d,%d\n",
			   tgLeft, tgTop, tgWidth, tgHeight, cWidth, cHeight,
			   fg, bg);
		 
	/* read (and ignore) data sub-blocks */
	do {
	  j = 0;
	  sbsize = NEXTBYTE;
	  while (j<sbsize) {
	    ch = NEXTBYTE;  j++;
	    if (DEBUG) fprintf(stderr,"%c", ch);
	  }
	} while (sbsize);
	if (DEBUG) fprintf(stderr,"\n\n");
      }
	  
	  
      else if (fn == 0xF9) {  /* Graphic Control Extension */
	int j, sbsize;
		 
	if (DEBUG) fprintf(stderr,"Graphic Control extension\n\n");
		 
	/*	SetISTR(ISTR_WARNING, 
		"Graphic Control Extension in GIF file.  Ignored.");*/
		 
	/* read (and ignore) data sub-blocks */
	do {
	  j = 0; sbsize = NEXTBYTE;
	  while (j<sbsize) { (void)NEXTBYTE;  j++; }
	} while (sbsize);
      }
      
	  
      else { /* unknown extension */
	int j, sbsize;
		 
	if (DEBUG) fprintf(stderr,"unknown GIF extension 0x%02x\n\n", fn);
		 
	/*	SetISTR(ISTR_WARNING, 
		"Unknown extension 0x%02x in GIF file.  Ignored.",fn);*/
		 
	/* read (and ignore) data sub-blocks */
	do {
	  j = 0; sbsize = NEXTBYTE;
	  while (j<sbsize) { (void)NEXTBYTE;  j++; }
	} while (sbsize);
      }
    }
   
    else if (block == IMAGESEP) break;   /* read an image */

    else if (block == TRAILER) {
      return( GifError("no image data found in GIF file") );
    }
   
    else return (GifError("Unknown block type found in file."));
  }


  /* read in values from the image descriptor */
  ch = NEXTBYTE;
  /* LeftOfs = ch + 0x100 * NEXTBYTE;*/
  ch = NEXTBYTE;
  /* TopOfs = ch + 0x100 * NEXTBYTE; */
  ch = NEXTBYTE;
  ch = NEXTBYTE;
  ch = NEXTBYTE;
  get_static_Width() = ch + 0x100 * NEXTBYTE;
  ch = NEXTBYTE;
  get_static_Height() = ch + 0x100 * NEXTBYTE;

  Misc = NEXTBYTE;
  Interlace = ((Misc & INTERLACEMASK) ? TRUE : FALSE);

  if (Misc & 0x80) {
    for (i=0; i< 1 << ((Misc&7)+1); i++) {
      get_static_r()[i] = NEXTBYTE;
      get_static_g()[i] = NEXTBYTE;
      get_static_b()[i] = NEXTBYTE;
    }
  }


  if (!HasColormap && !(Misc&0x80)) {
    /* no global or local colormap */
    /*    SetISTR(ISTR_WARNING,
	  "No colormap in this GIF file.  Assuming EGA colors.");*/
  }



  /* Start reading the raster data. First we get the intial code size
   * and compute decompressor constant values, based on this code size.
   */
  
  /*  SetISTR(ISTR_FORMAT, "GIF%s, %d bits per pixel, %sinterlaced.  (%d bytes)",
      (gif89) ? "89" : "87", BitsPerPixel, 
      Interlace ? "" : "non-", filesize);*/

  get_static_CodeSize() = NEXTBYTE;

  ClearCode = (1 << get_static_CodeSize());
  EOFCode = ClearCode + 1;
  FreeCode = FirstFree = ClearCode + 2;

  /* The GIF spec has it that the code size is the code size used to
   * compute the above values is the code size given in the file, but the
   * code size used in compression/decompression is the code size given in
   * the file plus one. (thus the ++).
   */
  get_static_CodeSize()++;
  InitCodeSize = get_static_CodeSize();
  MaxCode = (1 << get_static_CodeSize());
  get_static_ReadMask() = MaxCode - 1;

  /* UNBLOCK:
   * Read the raster data.  Here we just transpose it from the GIF array
   * to the Raster array, turning it from a series of blocks into one long
   * data stream, which makes life much easier for ReadCode().
   */
  
  ptr1 = get_static_Raster();
  do {
    ch = ch1 = NEXTBYTE;
    while (ch--) { 
		*ptr1 = NEXTBYTE; 
		ptr1++; }
    if ((ptr - get_static_RawGIF()) > filesize) {
      /*      SetISTR(ISTR_WARNING,
	      "This GIF file seems to be truncated.  Winging it.");*/
      break;
    }
  } while(ch1);
  ImageIO_free(get_static_RawGIF());  get_static_RawGIF() = NULL;


  if (DEBUG) {
    fprintf(stderr,"xv: LoadGIF() - picture is %dx%d, %d bits, %sinterlaced\n",
	    get_static_Width(), get_static_Height(), BitsPerPixel, Interlace ? "" : "non-");
  }
  

  /* Allocate the 'pic' */
  maxpixels = get_static_Width()*get_static_Height();
  im->xdim = get_static_Width();
  im->ydim = get_static_Height();
  im->zdim = 1;
  im->vdim = 3;
  im->wdim = 1;
  im->wordKind = WK_FIXED;
  im->sign = SGN_UNSIGNED;
  im->data = ImageIO_alloc(get_static_Width() * get_static_Height() * 3);
  get_static_org() = get_static_buf() = (unsigned char *) im->data;

  if (!get_static_org())
    return( GifError("not enough memory for image buffer") );


  /* Decompress the file, continuing until you see the GIF EOF code.
   * One obvious enhancement is to add checking for corrupt files here.
   */
  Code = ReadCode();
  while (Code != EOFCode) {
   /* Clear code sets everything back to its initial value, then reads the
     * immediately subsequent code as uncompressed data.
     */
   
	  if (Code == ClearCode) { 
      get_static_CodeSize() = InitCodeSize;
      MaxCode = (1 << get_static_CodeSize());
      get_static_ReadMask() = MaxCode - 1;
      FreeCode = FirstFree;
      Code = ReadCode();
      CurCode = OldCode = Code;
      FinChar = CurCode & BitMask;
      if (!Interlace) {
	*get_static_buf()++ = get_static_r()[FinChar];
	*get_static_buf()++ = get_static_g()[FinChar];
	*get_static_buf()++ = get_static_b()[FinChar];
      }
      else DoInterlace((byte)FinChar);
      npixels++;
    }

    else {
      /* If not a clear code, must be data: save same as CurCode and InCode */
	   
      /* if we're at maxcode and didn't get a clear, stop loading */
      if (FreeCode>=4096) { 
		  printf("freecode blew up\n"); 
	      break;
	  }
	   
      CurCode = InCode = Code;
	   
      /* If greater or equal to FreeCode, not in the hash table yet;
       * repeat the last character decoded
       */
	   
      if (CurCode >= FreeCode) {
	CurCode = OldCode;
	if (OutCount > 4096) {   
		printf("outcount1 blew up\n");  break; }
	OutCode[OutCount++] = FinChar;
      }
	   
      /* Unless this code is raw data, pursue the chain pointed to by CurCode
       * through the hash table to its end; each code in the chain puts its
       * associated output code on the output queue.
       */
	   
      while (CurCode > BitMask) {
		  if (OutCount > 4096) {
			  fprintf(stderr,"outcount2 blew up\n"); break;}   /* corrupt file */
	OutCode[OutCount++] = Suffix[CurCode];
	CurCode = Prefix[CurCode];
      }
	   
      if (OutCount > 4096) {  
		  printf("outcount blew up\n");  break; }
	   
      /* The last code in the chain is treated as raw data. */
	   
      FinChar = CurCode & BitMask;
      OutCode[OutCount++] = FinChar;
	   
      /* Now we put the data out to the Output routine.
       * It's been stacked LIFO, so deal with it that way...
       */

      /* safety thing:  prevent exceeding range of 'pic' */
      if (npixels + OutCount > maxpixels) OutCount = maxpixels-npixels;
	   
      npixels += OutCount;
      if (!Interlace) for (i=OutCount-1; i>=0; i--) {
          *get_static_buf()++ = get_static_r()[OutCode[i]];
          *get_static_buf()++ = get_static_g()[OutCode[i]];
          *get_static_buf()++ = get_static_b()[OutCode[i]];
      }
      else  for (i=OutCount-1; i>=0; i--) DoInterlace((byte)OutCode[i]);
      OutCount = 0;
	   
      /* Build the hash table on-the-fly. No table is stored in the file. */
      
      Prefix[FreeCode] = OldCode;
      Suffix[FreeCode] = FinChar;
      OldCode = InCode;
	   
      /* Point to the next slot in the table.  If we exceed the current
       * MaxCode value, increment the code size unless it's already 12.  If it
       * is, do nothing: the next code decompressed better be CLEAR
       */
	   
      FreeCode++;
      if (FreeCode >= MaxCode) {
	if (get_static_CodeSize() < 12) {
	  get_static_CodeSize()++;
	  MaxCode *= 2;
	  get_static_ReadMask() = (1 << get_static_CodeSize()) - 1;
	}
      }
    }
    Code = ReadCode();
    if (npixels >= maxpixels) break;
  }
  ImageIO_free(get_static_Raster());  get_static_Raster() = NULL;

  if (npixels != maxpixels) {
    /*    SetISTR(ISTR_WARNING,"This GIF file seems to be truncated.  Winging it.");*/
    if (!Interlace)
      memset(get_static_buf(), 0, 3*(maxpixels-npixels)); /* clear to EOBuffer */
  }
  /*  SetDirRButt(F_FORMAT, F_GIF);
      SetDirRButt(F_COLORS, F_FULLCOLOR);*/
  return 1;
}


/* Fetch the next code from the raster data stream.  The codes can be
 * any length from 3 to 12 bits, packed into 8-bit bytes, so we have to
 * maintain our location in the Raster array as a BIT Offset.  We compute
 * the byte Offset into the raster array by dividing this by 8, pick up
 * three bytes, compute the bit Offset into our 24-bit chunk, shift to
 * bring the desired code to the bottom, then mask it off and return it. 
 */

CGAL_INLINE_FUNCTION
static int ReadCode()
{
  int RawCode, ByteOffset;
  
  ByteOffset = get_static_BitOffset() / 8;
  RawCode = get_static_Raster()[ByteOffset] + (get_static_Raster()[ByteOffset + 1] << 8);
  if (get_static_CodeSize() >= 8)
    RawCode += ( ((int) get_static_Raster()[ByteOffset + 2]) << 16);
  RawCode >>= (get_static_BitOffset() % 8);
  get_static_BitOffset() += get_static_CodeSize();

  return(RawCode & get_static_ReadMask());
}


/***************************/
CGAL_INLINE_FUNCTION
static void DoInterlace(byte Index) {
  static byte *ptr = NULL;
  static int   oldYC = -1;
  
  if (oldYC != get_static_YC()) {
    ptr = get_static_org() + 3 * get_static_YC() * get_static_Width();
    oldYC = get_static_YC();
  }
  
  if (get_static_YC() < get_static_Height()) {
    *ptr++ = get_static_r()[Index];
    *ptr++ = get_static_g()[Index];
    *ptr++ = get_static_b()[Index];
  }
  
  /* Update the X-coordinate, and if it overflows, update the Y-coordinate */
  if (++get_static_XC() == get_static_Width()) {
    
    /* deal with the interlace as described in the GIF
     * spec.  Put the decoded scan line out to the screen if we haven't gone
     * past the bottom of it
     */
    
    get_static_XC() = 0;
    
    switch (get_static_Pass()) {
    case 0:
      get_static_YC() += 8;
      if (get_static_YC() >= get_static_Height()) { get_static_Pass()++; get_static_YC() = 4; }
      break;
      
    case 1:
      get_static_YC() += 8;
      if (get_static_YC() >= get_static_Height()) { get_static_Pass()++; get_static_YC() = 2; }
      break;
      
    case 2:
      get_static_YC() += 4;
      if (get_static_YC() >= get_static_Height()) { get_static_Pass()++; get_static_YC() = 1; }
      break;
      
    case 3:
      get_static_YC() += 2;  break;
      
    default:
      break;
    }
  }
}


      
/*****************************/
CGAL_INLINE_FUNCTION
static int GifError(const char *st) {
  fprintf(stderr,"readGifImage: error: %s\n",st);

  if (get_static_RawGIF() != NULL) ImageIO_free(get_static_RawGIF());
  if (get_static_Raster() != NULL) ImageIO_free(get_static_Raster());
  
  return -1;
}

#undef CGAL_GLOBAL_STATE_VAR

