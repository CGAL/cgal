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

#include <string.h>

#include "analyze.h"

/* compile time endianness */
/* replaced by _getEndianness()
   see below
#if (defined (_ALPHA_) || defined (_LINUX_))
#define ARCHITECTURE_ENDIANNESS END_LITTLE
#else
#define ARCHITECTURE_ENDIANNESS END_BIG
#endif
*/

/** Magic header for ANALYZE files written in little endian format */
#define ANALYZE_LE_MAGIC "\000\000\001\134"
/** Magic header for ANALYZE files written in big endian format */
#define ANALYZE_BE_MAGIC "\134\001\000\000"

#define DT_NONE			0 
#define DT_UNKNOWN              0 /*Unknown data type*/ 
#define DT_BINARY               1 /*Binary (1 bit per voxel)*/ 
#define DT_UNSIGNED_CHAR        2 /*Unsigned character (8 bits per voxel)*/ 
#define DT_SIGNED_SHORT         4 /*Signed short (16 bits per voxel)*/ 
#define DT_SIGNED_INT           8 /*Signed integer (32 bits per voxel)*/ 
#define DT_FLOAT                16 /*Floating point (32 bits per voxel)*/ 
#define DT_COMPLEX              32 /*Complex (64 bits per voxel; 2 floating point numbers) */
#define DT_DOUBLE               64 /*Double precision (64 bits per voxel)*/ 
#define DT_RGB                  128 /* */
#define DT_ALL                  255 /* */

struct header_key                       /*      header_key       */
    {                                           /* off + size*/
        int sizeof_hdr;                         /* 0 + 4     */
        char data_type[10];                     /* 4 + 10    */
        char db_name[18];                       /* 14 + 18   */
        int extents;                            /* 32 + 4    */
        short int session_error;                /* 36 + 2    */
        char regular;                           /* 38 + 1    */
        char hkey_un0;                          /* 39 + 1    */
    };                                          /* total=40  */

struct image_dimension                  /*      image_dimension  */
    {                                           /* off + size*/
        short int dim[8];                       /* 0 + 16    */
        char vox_units[4];                      /* 16 + 4    */
        char cal_units[8];                      /* 20 + 4    */
        short int unused1;                      /* 24 + 2    */
        short int datatype;                     /* 30 + 2    */
        short int bitpix;                       /* 32 + 2    */
        short int dim_un0;                      /* 34 + 2    */
        float pixdim[8];                        /* 36 + 32   */
                        /* 
                                pixdim[] specifies the voxel dimensions:
                                pixdim[1] - voxel width
                                pixdim[2] - voxel height
                                pixdim[3] - interslice distance
                                        ..etc
                        */
        float vox_offset;                       /* 68 + 4    */
        float funused1;                         /* 72 + 4    */
        float funused2;                         /* 76 + 4    */
        float funused3;                         /* 80 + 4    */
        float cal_max;                          /* 84 + 4    */
        float cal_min;                          /* 88 + 4    */
        int compressed;                         /* 92 + 4    */
        int verified;                           /* 96 + 4    */
        int glmax, glmin;                       /* 100 + 8   */
    };                                          /* total=108 */
         
struct data_history                     /*      data_history     */
    {                                           /* off + size*/
        char descrip[80];                       /* 0 + 80    */
        char aux_file[24];                      /* 80 + 24   */
        char orient;                            /* 104 + 1   */
        char originator[10];                    /* 105 + 10  */
        char generated[10];                     /* 115 + 10  */
        char scannum[10];                       /* 125 + 10  */
        char patient_id[10];                    /* 135 + 10  */
        char exp_date[10];                      /* 145 + 10  */
        char exp_time[10];                      /* 155 + 10  */
        char hist_un0[3];                       /* 165 + 3   */
        int views;                              /* 168 + 4   */
        int vols_added;                         /* 172 + 4   */
        int start_field;                        /* 176 + 4   */
        int field_skip;                         /* 180 + 4   */
        int omax,omin;                          /* 184 + 8   */
        int smax,smin;                          /* 192 + 8   */
    };                                        


struct dsr                                   /*      dsr              */
   {                                         /* off + size*/
     struct header_key hk;                   /* 0 + 40    */
     struct image_dimension dime;            /* 40 + 108  */
     struct data_history hist;               /* 148 + 200 */
   };                                        /* total=348 */

 


/*---------------- _swapLong ------------------------------------------------*/

static void _swapLong( unsigned char *pntr)
{
  unsigned char b0, b1, b2, b3;
  
  b0 = *pntr;
  b1 = *(pntr+1);
  b2 = *(pntr+2);
  b3 = *(pntr+3);
  
  *pntr = b3;
  *(pntr+1) = b2;
  *(pntr+2) = b1;
  *(pntr+3) = b0;
}
        
/*---------------- _swapShort -----------------------------------------------*/

static void _swapShort( unsigned char *pntr)
{
  unsigned char b0, b1;
  
  b0 = *pntr;
  b1 = *(pntr+1);

  *pntr = b1;
  *(pntr+1) = b0;
}

/*---------------- _swapAnalyzeHdr ------------------------------------------*/

static void _swapAnalyzeHdr( struct dsr *pntr)
{
  
  _swapLong((unsigned char*) &pntr->hk.sizeof_hdr) ;
  _swapLong((unsigned char*) &pntr->hk.extents) ;
  _swapShort((unsigned char*) &pntr->hk.session_error) ;
  _swapShort((unsigned char*) &pntr->dime.dim[0]) ;
  _swapShort((unsigned char*) &pntr->dime.dim[1]) ;
  _swapShort((unsigned char*) &pntr->dime.dim[2]) ;
  _swapShort((unsigned char*) &pntr->dime.dim[3]) ;
  _swapShort((unsigned char*) &pntr->dime.dim[4]) ;
  _swapShort((unsigned char*) &pntr->dime.dim[5]) ;
  _swapShort((unsigned char*) &pntr->dime.dim[6]) ;
  _swapShort((unsigned char*) &pntr->dime.dim[7]) ;
  _swapShort((unsigned char*) &pntr->dime.unused1) ;
  _swapShort((unsigned char*) &pntr->dime.datatype) ;
  _swapShort((unsigned char*) &pntr->dime.bitpix) ;
  _swapLong((unsigned char*) &pntr->dime.pixdim[0]) ;
  _swapLong((unsigned char*) &pntr->dime.pixdim[1]) ;
  _swapLong((unsigned char*) &pntr->dime.pixdim[2]) ;
  _swapLong((unsigned char*) &pntr->dime.pixdim[3]) ;
  _swapLong((unsigned char*) &pntr->dime.pixdim[4]) ;
  _swapLong((unsigned char*) &pntr->dime.pixdim[5]) ;
  _swapLong((unsigned char*) &pntr->dime.pixdim[6]) ;
  _swapLong((unsigned char*) &pntr->dime.pixdim[7]) ;
  _swapLong((unsigned char*) &pntr->dime.vox_offset) ;
  _swapLong((unsigned char*) &pntr->dime.funused1) ;
  _swapLong((unsigned char*) &pntr->dime.funused2) ;
  _swapLong((unsigned char*) &pntr->dime.cal_max) ;
  _swapLong((unsigned char*) &pntr->dime.cal_min) ;
  _swapLong((unsigned char*) &pntr->dime.compressed) ;
  _swapLong((unsigned char*) &pntr->dime.verified) ;
  _swapShort((unsigned char*) &pntr->dime.dim_un0) ;
  _swapLong((unsigned char*) &pntr->dime.glmax) ;
  _swapLong((unsigned char*) &pntr->dime.glmin) ;
}
PTRIMAGE_FORMAT createAnalyzeFormat() {
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));

  f->testImageFormat=&testAnalyzeHeader;
  f->readImageHeader=&readAnalyzeHeader;
  f->writeImage=&writeAnalyze;
  strcpy(f->fileExtension,".hdr,.hdr.gz,.img,.img.gz");
  strcpy(f->realName,"Analyze");
  return f;
}

int testAnalyzeHeader(char *magic,const char *) {
   /* opened image is an ANALYZE */
  if( !memcmp(magic,ANALYZE_LE_MAGIC,4) ||
      !memcmp(magic,ANALYZE_BE_MAGIC,4)) 
    return 0;
  else 
    return -1;
}

int writeAnalyze( char *name, _image* im) {
  char *outputName;
  int length, extLength=0, res;


  length=strlen(name);
  outputName= (char *)ImageIO_alloc(length+8);
  
  if ( strncmp( name+length-4, ".hdr", 4 ) == 0 ) {
    extLength = 4;
  }
  else if ( strncmp( name+length-4, ".img", 4 ) == 0 ) {
    extLength = 4;
  }
  else if ( strncmp( name+length-7, ".img.gz", 7 ) == 0 ) {
    extLength = 7;
  }
  else if ( strncmp( name+length-7, ".hdr.gz", 7 ) == 0 ) {
    extLength = 7;
  }

  strncpy( outputName, name, length-extLength );
  if ( strncmp( name+length-7, ".hdr.gz", 7 ) == 0 )
    strcpy( outputName+length-extLength, ".hdr.gz" );
  else
    strcpy( outputName+length-extLength, ".hdr" );
  
  _openWriteImage(im, outputName);
  if( !im->fd ) {
    fprintf(stderr, "writeAnalyze: error: unable to open file \'%s\'\n", outputName);
    if ( outputName != NULL ) ImageIO_free( outputName );
    return ImageIO_OPENING;
  }

  res = writeAnalyzeHeader(im);
  if ( res < 0 ) {
    fprintf(stderr, "writeAnalyze: error: unable to write header of \'%s\'\n",
	    outputName);
    if ( outputName != NULL ) ImageIO_free( outputName );
    ImageIO_close( im );
    im->fd = NULL;
    im->openMode = OM_CLOSE;
    return( res );
  }
  
  ImageIO_close(im);
  
  strncpy( outputName, name, length-extLength );
  if ( strncmp( name+length-3, ".gz", 3 ) == 0 ) {
    strcpy( outputName+length-extLength, ".img.gz" );
  }
  else {
    strcpy( outputName+length-extLength, ".img" );
  }

  _openWriteImage(im, outputName);

  if( !im->fd ) {
    fprintf(stderr, "writeAnalyze: error: unable to open file \'%s\'\n", outputName);
    if ( outputName != NULL ) ImageIO_free( outputName );
    return ImageIO_OPENING;
  }

  res = writeAnalyzeData(im);
  if (res < 0) {
    fprintf(stderr, "writeAnalyze: error: unable to write data in \'%s\'\n",
	    outputName );
    ImageIO_close( im );
    im->fd = NULL;
    im->openMode = OM_CLOSE;
    return( res );
  }

  if ( outputName != NULL ) ImageIO_free( outputName );
  ImageIO_close( im );
  im->fd = NULL;
  im->openMode = OM_CLOSE;

  return ( res );
}

/* 
   return:
   -1: error
   0: success
 */

int _readAnalyzeHeader( _image* im, const char* name, 
			struct dsr *analyzeHeader )
{

   unsigned int i ;
   /* compile time endianness */
   ENDIANNESS ARCHITECTURE_ENDIANNESS = _getEndianness();

   if(im->openMode != OM_CLOSE) {
      
     ImageIO_read( im, analyzeHeader, sizeof(struct dsr) );
      
      if( analyzeHeader->hk.sizeof_hdr == sizeof(struct dsr) )
      {
	 im->endianness = ARCHITECTURE_ENDIANNESS ;
      }
      else
      {
	_swapAnalyzeHdr( analyzeHeader );
	 if( analyzeHeader->hk.sizeof_hdr != sizeof(struct dsr) )
	 {
	    fprintf (stderr,
		     "_readAnalyzeHeader: error: unknown magic (%d)...\n",
		     analyzeHeader->hk.sizeof_hdr );
	    return -1;
	 }
	 if( ARCHITECTURE_ENDIANNESS == END_LITTLE )
	 {
	    im->endianness = END_BIG;
	 }
	 else
	 {
	    im->endianness = END_LITTLE;
	 }
      }
      
      if ( analyzeHeader->dime.dim[0] > 4 ) {
	 fprintf (stderr,
		  "_readAnalyzeHeader: error: dimensionality not supported (%d)...\n",
		  analyzeHeader->dime.dim[0] );
	 return -1;
      }
      
      im->xdim = analyzeHeader->dime.dim[1];
      im->ydim = analyzeHeader->dime.dim[2];
      im->zdim = analyzeHeader->dime.dim[3];

      /* 0 time-points is a convention for one volume only at MNI */
      /* Corrected by X. Pennec following a bug report by Irina Kezele at MNI */
      if( analyzeHeader->dime.dim[4] == 0 ){
	fprintf (stderr,
		 "_readAnalyzeHeader: warning: time dimension / number if volume (dim[4]) is 0. Assuming this means 1 (otherwise there would be no image...). \n");
	analyzeHeader->dime.dim[4] = 1; 
      }

      /* Analyze doesn't support vector images.
	 The forth dimension relates to time. */
      if( analyzeHeader->dime.dim[4] != 1 )
      {
	 fprintf (stderr,
		  "_readAnalyzeHeader: error: time dimension not supported (%d)...\n",
		  analyzeHeader->dime.dim[4] );
	 return -1;
      }
      im->vectMode = VM_SCALAR;
      
      im->vx = analyzeHeader->dime.pixdim[1];
      im->vy = analyzeHeader->dime.pixdim[2];
      im->vz = analyzeHeader->dime.pixdim[3];

      if( im->vx == 0.0 ) im->vx = 1.0 ;
      if( im->vy == 0.0 ) im->vy = im->vx ;
      if( im->vz == 0.0 ) im->vz = im->vy ;
      
      switch(analyzeHeader->dime.datatype)
      {
	 case DT_BINARY:
	 case DT_UNSIGNED_CHAR:
	 case DT_SIGNED_SHORT:
	 case DT_SIGNED_INT:
	 case DT_FLOAT:
	 case DT_COMPLEX:
	 case DT_DOUBLE:
	    im->vdim = 1;
	    break ;
	    
	 case DT_RGB:
	    im->vdim = 3;
	    break ;
	    
	 default:
	    fprintf (stderr,
		     "_readAnalyzeHeader: error: data type not supported (%d)...\n",
		     analyzeHeader->dime.datatype );
	    return -1;
      }
      
      switch(analyzeHeader->dime.datatype)
      {
	 case DT_BINARY:
	 case DT_UNSIGNED_CHAR:
	 case DT_SIGNED_SHORT:
	 case DT_SIGNED_INT:
	 case DT_RGB:
	    im->wordKind = WK_FIXED;
	    break ;
	    
	 case DT_FLOAT:
	 case DT_COMPLEX:
	 case DT_DOUBLE:
	    im->wordKind = WK_FLOAT;
	    break ;
	    
	 default:
	    fprintf (stderr,
		     "_readAnalyzeHeader: error: data type not supported (%d)...\n",
		     analyzeHeader->dime.datatype );
	    return -1;
      }
      
      switch(analyzeHeader->dime.datatype)
      {
	 case DT_BINARY:
	 case DT_UNSIGNED_CHAR:
	 case DT_RGB:
	    im->sign = SGN_UNSIGNED;
	    break ;
	    
	 case DT_SIGNED_SHORT:
	 case DT_SIGNED_INT:
	 case DT_FLOAT:
	 case DT_COMPLEX:
	 case DT_DOUBLE:
	    im->sign = SGN_SIGNED;
	    break ;
	    
	 default:
	    fprintf (stderr,
		     "_readAnalyzeHeader: error: data type not supported (%d)...\n",
		     analyzeHeader->dime.datatype );
	    return -1;
      }
      
      im->wdim = analyzeHeader->dime.bitpix;
      if( analyzeHeader->dime.datatype == DT_RGB )
      {
	 im->wdim /= 3 ;
      }
      if(im->wdim != 8 && im->wdim != 16 && im->wdim != 32 && im->wdim != 64)
      {
	 fprintf (stderr,
		  "_readAnalyzeHeader: error: pixel size not supported (%d)...\n",
		  analyzeHeader->dime.bitpix );
	 return -1;
      }
      im->wdim >>= 3 ;
      
      /* There are 17 optional data fields 
	 be careful in the allocation
       */
      im->nuser = 1 + 17 ;
      im->user = (char **) ImageIO_alloc(im->nuser * sizeof(char *));
      for ( i=0; i<im->nuser; i++ ) im->user[i] = NULL;
      i = 0 ;
      
      im->user[i] = (char *) ImageIO_alloc((strlen("Data lost in the Analyze -> ImageIO conversion:") + 1));
      sprintf( im->user[i++], "Data lost in the Analyze -> ImageIO conversion:" );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  descrip: ") + 1 + strlen(analyzeHeader->hist.descrip) ));
      sprintf( im->user[i++], "  descrip: %s", analyzeHeader->hist.descrip );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  aux_file: ") + 1 + strlen(analyzeHeader->hist.descrip) ));
      sprintf( im->user[i++], "  aux_file: %s", analyzeHeader->hist.descrip );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  orient: ") + 1+ 2));
      sprintf( im->user[i++], "  orient: %d", analyzeHeader->hist.orient );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  originator: ") + 1 + strlen(analyzeHeader->hist.originator) ));
      sprintf( im->user[i++], "  originator: %s", analyzeHeader->hist.originator );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  generated: ") + 1 + strlen(analyzeHeader->hist.generated) ));
      sprintf( im->user[i++], "  generated: %s", analyzeHeader->hist.generated );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  scannum: ") + 1 + strlen(analyzeHeader->hist.scannum) ));
      sprintf( im->user[i++], "  scannum: %s", analyzeHeader->hist.scannum );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  patient_id: ") + 1 + strlen(analyzeHeader->hist.patient_id) ));
      sprintf( im->user[i++], "  patient_id: %s", analyzeHeader->hist.patient_id );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  exp_date: ") + 1 + strlen(analyzeHeader->hist.exp_date) ));
      sprintf( im->user[i++], "  exp_date: %s", analyzeHeader->hist.exp_date );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  exp_time: ") + 1 + strlen(analyzeHeader->hist.exp_time) ));
      sprintf( im->user[i++], "  exp_time: %s", analyzeHeader->hist.exp_time );

      /* A 32 bit int doesn't print on more than 11 chars */
      im->user[i] = (char *) ImageIO_alloc((strlen("  views: ") + 11 + 1));
      sprintf( im->user[i++], "  views: %d", analyzeHeader->hist.views );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  vols_added: ") + 11 + 1));
      sprintf( im->user[i++], "  vols_added: %d", analyzeHeader->hist.vols_added );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  start_field: ") + 11 + 1));
      sprintf( im->user[i++], "  start_field: %d", analyzeHeader->hist.start_field );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  field_skip: ") + 11 + 1));
      sprintf( im->user[i++], "  field_skip: %d", analyzeHeader->hist.field_skip );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  omax: ") + 11 + 1));
      sprintf( im->user[i++], "  omax: %d", analyzeHeader->hist.omax );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  omin: ") + 11 + 1));
      sprintf( im->user[i++], "  omin: %d", analyzeHeader->hist.omin );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  smax: ") + 11 + 1));
      sprintf( im->user[i++], "  smax: %d", analyzeHeader->hist.smax );
      
      im->user[i] = (char *) ImageIO_alloc((strlen("  smin: ") + 11 + 1));
      sprintf( im->user[i++], "  smin: %d", analyzeHeader->hist.smin );


      /* header is read. close header file and open data file. */
      if( name != NULL ) {

	int length = strlen(name) ;
	char* data_filename = (char *) ImageIO_alloc(length+4) ;
	
	if( strcmp( name+length-4, ".hdr" ) )
	  {
	    fprintf (stderr,
		     "_readAnalyzeHeader: error: file header extension must be .hdr\n");
	    ImageIO_free( data_filename );
	    return -1;
	  }
	 
	ImageIO_close(im);
	
	strcpy(data_filename,name);
	strcpy(data_filename+length-3, "img.gz");
	_openReadImage(im,data_filename);
	
	if(!im->fd) {

	  strcpy(data_filename,name);
	  strcpy(data_filename+length-3, "img");
	  _openReadImage(im,data_filename);
	  if(!im->fd) {
	    fprintf(stderr, "_readAnalyzeHeader: error: unable to open data file \'%s\'\n", data_filename);
	    ImageIO_free( data_filename );
	    return -1;
	    
	  }
	}
	ImageIO_free( data_filename );
      }
      
      /* check header validity */
      if(im->xdim > 0 && im->ydim > 0 && im->zdim > 0 && im->vdim > 0 &&
	 im->vx > 0.0 && im->vy > 0.0 && im->vz > 0.0 &&
	 (im->wordKind == WK_FLOAT || (im->wordKind == WK_FIXED &&
				       im->sign != SGN_UNKNOWN)) &&
	 im->endianness != END_UNKNOWN) {
	 return 0;
      }
      else return -1;
   }
   else return -1;
}







int readAnalyzeHeader( const char* name, _image* im)
{
  struct dsr analyzeHeader ;
  return( _readAnalyzeHeader( im, name, &analyzeHeader ) );
}



int
writeAnalyzeHeader( const _image* im )
{
  const char *proc = "writeAnalyzeHeader";
  struct dsr hdr;
  int i ;
  int imin = 0;
  int imax = 0;
 
   memset(&hdr,0, sizeof(struct dsr));
   
   for(i=0;i<8;i++) {
      hdr.dime.pixdim[i] = 0.0;
   }
   hdr.dime.vox_offset  = 0.0;
   hdr.dime.funused1    = 0.0;
   hdr.dime.funused2    = 0.0;
   hdr.dime.funused3    = 0.0;
   hdr.dime.cal_max     = 0.0;
   hdr.dime.cal_min     = 0.0;

   hdr.dime.dim[0] = 4;
   hdr.dime.dim[1] = im->xdim;
   hdr.dime.dim[2] = im->ydim;
   hdr.dime.dim[3] = im->zdim;
   hdr.dime.dim[4] = 1 ;

   if ( im->wordKind == WK_FIXED && im->sign == SGN_UNSIGNED ) {
     if( im->wdim == 1 ) {
       
	if ( im->vdim == 1 ) {
	  hdr.dime.datatype = DT_UNSIGNED_CHAR ;
	}
	else if ( im->vdim == 3 ) {
	  hdr.dime.datatype = DT_RGB ;
	}
	else {
	  fprintf( stderr, "%s: unsupported image type\n", proc );
	  return -1;
	}
	
	{
	  unsigned char *buf = (unsigned char *)im->data;
	  int size = im->xdim * im->ydim * im->zdim * im->vdim;
	  imin = imax = *buf;
	  for (i=0; i<size; i++, buf++) {
	    if ( imax < *buf ) imax = *buf;
	    if ( imin > *buf ) imin = *buf;
	  }
	}
	
      }
      else if ( im->wdim == 2 ) {
	if ( im->vdim == 1 ) {
	  unsigned short int *buf = (unsigned short int*)im->data;
	  int size = im->xdim * im->ydim *im->zdim;
	  int i;
	  imin = imax = *buf;
	  for (i=0; i<size; i++, buf++) {
	    if ( imax < *buf ) imax = *buf;
	    if ( imin > *buf ) imin = *buf;
	  }
	  if ( imax < 32768 ) {
	    hdr.dime.datatype = DT_SIGNED_SHORT ;
	  } 
	  else {
	    fprintf( stderr, "%s: conversion from unsigned short to short impossible, max=%d\n", proc, imax );
	    return -1;
	  }
	}
	else {
	  fprintf( stderr, "%s: unsupported image type\n", proc );
	  return -1;
	}
      }
      
      else {
	  fprintf( stderr, "%s: unsupported image type\n", proc );
	  return -1;
      }
      
   } /* if ( im->wordKind == WK_FIXED && im->sign == SGN_UNSIGNED ) */

   else if( im->wordKind == WK_FIXED && im->sign == SGN_SIGNED ) {
     if ( im->vdim != 1 ) {
       fprintf( stderr, "%s: unsupported image type\n", proc );
       return -1;
     }
     if( im->wdim == 2 ) {
       short int *buf = (short int*)im->data;
       int size = im->xdim * im->ydim *im->zdim;
       int i;
       imin = imax = *buf;
       for (i=0; i<size; i++, buf++) {
	 if ( imax < *buf ) imax = *buf;
	 if ( imin > *buf ) imin = *buf;
       }
       hdr.dime.datatype = DT_SIGNED_SHORT ;
     }
     else if( im->wdim == 4 ) {
       int *buf = (int*)im->data;
       int size = im->xdim * im->ydim *im->zdim;
       int i;
       imin = imax = *buf;
       for (i=0; i<size; i++, buf++) {
	 if ( imax < *buf ) imax = *buf;
	 if ( imin > *buf ) imin = *buf;
       }
       hdr.dime.datatype = DT_SIGNED_INT ;
     }
     else {
       fprintf( stderr, "%s: unsupported image type\n", proc );
       return -1;
     }
   }
   else if( im->wordKind == WK_FLOAT ) {
     if ( im->vdim != 1 ) {
       fprintf( stderr, "%s: unsupported image type\n", proc );
       return -1;
     }
     if( im->wdim == 4 ) {
       hdr.dime.datatype = DT_FLOAT ;
     }
     else if( im->wdim == 8 ) {
       hdr.dime.datatype = DT_DOUBLE ;
     }
     else {
       fprintf( stderr, "%s: unsupported image type\n", proc );
       return -1;
     }
   }
   else
   {
      fprintf( stderr, "%s: unsupported image type\n", proc );
      return -1;
   }
	 
   hdr.dime.bitpix = 8*im->wdim*im->vdim ;

   hdr.hk.regular = 'r';
   hdr.hk.sizeof_hdr = sizeof(struct dsr);

   /* this is probably bad and should be changed to the
      real values, but I'm too lazy to do it now. AG */
   hdr.dime.glmax  = 0 ;  /* maximum voxel value  */
   hdr.dime.glmin  = 0 ;  /* minimum voxel value */
   
   /* corrected GM
    */
   hdr.dime.glmax  = imax ;  /* maximum voxel value  */
   hdr.dime.glmin  = imin ;  /* minimum voxel value */
    
/*     Set the voxel dimension fields: 
       A value of 0.0 for these fields implies that the value is unknown.
         Change these values to what is appropriate for your data
         or pass additional command line arguments     */      
         
    hdr.dime.pixdim[1] = (float)im->vx;
    hdr.dime.pixdim[2] = (float)im->vy;
    hdr.dime.pixdim[3] = (float)im->vz;
    
/*   Assume zero offset in .img file, byte at which pixel
       data starts in the image file */

    hdr.dime.vox_offset = 0.0; 
    
/*   Planar Orientation;    */
/*   Movie flag OFF: 0 = transverse, 1 = coronal, 2 = sagittal
     Movie flag ON:  3 = transverse, 4 = coronal, 5 = sagittal  */  

    hdr.hist.orient     = 0;  
    
/*   up to 3 characters for the voxels units label; i.e. mm., um., cm. */

    strcpy(hdr.dime.vox_units,"mm.");
   
/*   up to 7 characters for the calibration units label; i.e. HU */

    strcpy(hdr.dime.cal_units," ");  
    
/*     Calibration maximum and minimum values;  
       values of 0.0 for both fields imply that no 
       calibration max and min values are used    */

    hdr.dime.cal_max = 0.0; 
    hdr.dime.cal_min = 0.0;

    if(ImageIO_write(im, &hdr, sizeof(struct dsr)) !=sizeof(struct dsr) )
       return -1;

    return 1 ;
}





/* Writes the given image body in an already opened file.*/
int writeAnalyzeData(const _image *im) {
  unsigned int lineSize = im->wdim * im->xdim * im->vdim ;
  unsigned long size = lineSize * im->ydim * im->zdim;
  unsigned int nwrt ;

  if(im->openMode != OM_CLOSE) {

#ifdef _REVERSE_LINES_IN_ANALYZE_
    char* data = (char *)im->data ;
    char* buf = data + size - lineSize ;
    
    while( buf >= data )
      {
	nwrt = ImageIO_write(im, buf, lineSize);
	if(nwrt != lineSize) return -1;
	buf -= lineSize ;
      }
#else
    nwrt = ImageIO_write(im, im->data, size);
    if(nwrt != size) return -1;
#endif    

    return 1 ;
  }
  else return -1;
}


/* Writes the given image body in an already opened file.*/
int printAnalyzeHeader( const char* name )
{
  _image *im;
  struct dsr analyzeHeader ;

  im = _initImage();
  _openReadImage(im, name);
  if(!im->fd) {
    fprintf(stderr, "printAnalyzeHeader: error: unable to open file \'%s\'\n", name);
    _freeImage(im);
    return -1;
  }

  if ( _readAnalyzeHeader( im, name, &analyzeHeader ) != 1 ) {
    fprintf(stderr, "printAnalyzeHeader: error: unable to read header in file \'%s\'\n", name);
    _freeImage(im);
    return -1;
  }


  


  ImageIO_close(im);
  im->fd = NULL;
  im->openMode = OM_CLOSE;
  _freeImage(im);
  return( 1 );


}
