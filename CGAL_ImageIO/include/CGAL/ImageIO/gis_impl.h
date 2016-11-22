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

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

#include <string>
#include <sstream>
#include <cstring>

#include <CGAL/ImageIO/inr.h>
#include <CGAL/ImageIO/fgetns.h>

#define _LGTH_STRING_ 1024

CGAL_INLINE_FUNCTION
PTRIMAGE_FORMAT createGisFormat() {
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));

  f->testImageFormat=&testGisHeader;
  f->readImageHeader=&readGisHeader;
  f->writeImage=&writeGis;
  strcpy(f->fileExtension,".dim,.dim.gz,.ima,.ima.gz");
  strcpy(f->realName,"Gis");
  return f;
}
CGAL_INLINE_FUNCTION
int writeGis( char *name, _image* im) {
  char *outputName;
  std::size_t length, extLength=0;
  int res;
  std::size_t done;

  length=strlen(name);
  outputName= (char *)ImageIO_alloc(length+8);
  
  if ( strncmp( name+length-4, ".dim", 4 ) == 0 ) {
    extLength = 4;
  }
  else if ( strncmp( name+length-4, ".ima", 4 ) == 0 ) {
    extLength = 4;
  }
  else if ( strncmp( name+length-7, ".ima.gz", 7 ) == 0 ) {
    extLength = 7;
  }
  else if ( strncmp( name+length-7, ".dim.gz", 7 ) == 0 ) {
    extLength = 7;
  }

  strncpy( outputName, name, length-extLength );
  if ( strncmp( name+length-7, ".dim.gz", 7 ) == 0 )
    strcpy( outputName+length-extLength, ".dim.gz" );
  else
    strcpy( outputName+length-extLength, ".dim" );

  _openWriteImage(im, outputName);
  if( !im->fd ) {
    fprintf(stderr, "writeGis: error: unable to open file \'%s\'\n", outputName);
    if ( outputName != NULL ) ImageIO_free( outputName );
    return ImageIO_OPENING;
  }

  if ( !writeGisHeader(im) ) {
    fprintf(stderr, "writeGis: error: unable to write header of \'%s\'\n",
	    outputName);
    if ( outputName != NULL ) ImageIO_free( outputName );
    ImageIO_close( im );
    im->fd = NULL;
    im->openMode = OM_CLOSE;
    return -1;
  }

  ImageIO_close(im);
  
  strncpy( outputName, name, length-extLength );
  if ( strncmp( name+length-3, ".gz", 3 ) == 0 ) {
    strcpy( outputName+length-extLength, ".ima.gz" );
  }
  else {
    strcpy( outputName+length-extLength, ".ima" );
  }

  _openWriteImage(im, outputName);

  if( !im->fd ) {
    fprintf(stderr, "writeGis: error: unable to open file \'%s\'\n", outputName);
    if ( outputName != NULL ) ImageIO_free( outputName );
    return ImageIO_OPENING;
  }

  if ( im->dataMode == DM_ASCII ) {
    std::size_t i, j, n, size;
    char *str = (char*)ImageIO_alloc( _LGTH_STRING_+1 );
    size = im->xdim * im->ydim * im->zdim * im->vdim;
    n = ( im->xdim < 16 ) ? im->xdim : 16;
    i = 0;

    switch( im->wordKind ) {
    default :
      fprintf(stderr, "writeGis: such word kind not handled in ascii mode for file \'%s\'\n", outputName);
      if ( outputName != NULL ) ImageIO_free( outputName );
      return( -3 );
    case WK_FIXED :
      switch ( im->wdim ) {
      default :
	fprintf(stderr, "writeGis: such word dim not handled in ascii mode for file \'%s\'\n", outputName);
	if ( outputName != NULL ) ImageIO_free( outputName );
	return( -3 );
      case 1 :
	switch ( im->sign ) {
	default :
	  fprintf(stderr, "writeGis: such sign not handled in ascii mode for file \'%s\'\n", outputName);
	  if ( outputName != NULL ) ImageIO_free( outputName );
	  return( -3 );
	case SGN_UNSIGNED :
	  {
	    unsigned char *theBuf = ( unsigned char * )im->data;
	    do {
	      memset( str, 0, _LGTH_STRING_ );
	      for ( j=0; j<n && i<size; j++, i++ ) {
		sprintf( str+strlen(str), "%d", theBuf[i] );
		if ( j<n && i<size ) sprintf( str+strlen(str), " " );
	      }
	      sprintf( str+strlen(str), "\n" );
              done = ImageIO_write( im, str, strlen( str )  );
              res = (done == strlen( str )) ? int(done) : -1;
	      if ( res  <= 0 ) {
		fprintf(stderr, "writeGis: error when writing data in \'%s\'\n", outputName);
		if ( outputName != NULL ) ImageIO_free( outputName );
		return( -3 );
	      }
	    } while ( i < size );
	  }
	  break;
	case SGN_SIGNED :
	  {
	    char *theBuf = ( char * )im->data;
	    do {
	      memset( str, 0, _LGTH_STRING_ );
	      for ( j=0; j<n && i<size; j++, i++ ) {
		sprintf( str+strlen(str), "%d", theBuf[i] );
		if ( j<n && i<size ) sprintf( str+strlen(str), " " );
	      }
	      sprintf( str+strlen(str), "\n" );
              done = ImageIO_write( im, str, strlen( str )  );
              res = (done == strlen( str )) ? int(done) : -1;
	      if ( res  <= 0 ) {
		fprintf(stderr, "writeGis: error when writing data in \'%s\'\n", outputName);
		if ( outputName != NULL ) ImageIO_free( outputName );
		return( -3 );
	      }
	    } while ( i < size );
	  }
	  break;
	} /* end of switch ( im->sign ) */
	break;
      case 2 :
	switch ( im->sign ) {
	default :
	  fprintf(stderr, "writeGis: such sign not handled in ascii mode for file \'%s\'\n", outputName);
	  if ( outputName != NULL ) ImageIO_free( outputName );
	  return( -3 );
	case SGN_UNSIGNED :
	  {
	    unsigned short int *theBuf = ( unsigned short int * )im->data;
	    do {
	      memset( str, 0, _LGTH_STRING_ );
	      for ( j=0; j<n && i<size; j++, i++ ) {
		sprintf( str+strlen(str), "%d", theBuf[i] );
		if ( j<n && i<size ) sprintf( str+strlen(str), " " );
	      }
	      sprintf( str+strlen(str), "\n" );
	      done = ImageIO_write( im, str, strlen( str )  );
              res = (done == strlen( str )) ? int(done) : -1;
	      if ( res  <= 0 ) {
		fprintf(stderr, "writeGis: error when writing data in \'%s\'\n", outputName);
		if ( outputName != NULL ) ImageIO_free( outputName );
		return( -3 );
	      }
	    } while ( i < size );
	  }
	  break;
	case SGN_SIGNED :
	  {
	    short int *theBuf = ( short int * )im->data;
	    do {
	      memset( str, 0, _LGTH_STRING_ );
	      for ( j=0; j<n && i<size; j++, i++ ) {
		sprintf( str+strlen(str), "%d", theBuf[i] );
		if ( j<n && i<size ) sprintf( str+strlen(str), " " );
	      }
	      sprintf( str+strlen(str), "\n" );
	      done = ImageIO_write( im, str, strlen( str )  );
              res = (done == strlen( str )) ? int(done) : -1;
	      if ( res  <= 0 ) {
		fprintf(stderr, "writeGis: error when writing data in \'%s\'\n", outputName);
		if ( outputName != NULL ) ImageIO_free( outputName );
		return( -3 );
	      }
	    } while ( i < size );
	  }
	  break;
	} /* end of switch ( im->sign ) */
	break;	
      } /* end of switch ( im->wdim ) */
    } /* end of switch( im->wordKind ) */

    ImageIO_free( str ); 
    if (outputName != NULL) ImageIO_free(outputName);
    return static_cast<int>(res);
  }
  else {
    bool ret = _writeInrimageData(im);
    if (outputName != NULL) ImageIO_free(outputName);
    return (ret ? 1 : -1);
  }
}

CGAL_INLINE_FUNCTION
int testGisHeader(char *,const char *name) {
   if (( !strncmp(name+strlen(name)-4, ".dim", 4)) || 
       ( !strncmp(name+strlen(name)-4, ".ima", 4)) || 
       ( !strncmp(name+strlen(name)-7, ".ima.gz", 7)) || 
       ( !strncmp(name+strlen(name)-7, ".dim.gz", 7)) )
     return 0;
  else 
    return -1;
}



CGAL_INLINE_FUNCTION
int readGisHeader( const char* name,_image* im)
{
  char *s, *str = NULL;
  int status;
  int n=0, nusermax = 20;

  str = (char*)ImageIO_alloc( _LGTH_STRING_+1 );
  
  if ( !fgetns(str, _LGTH_STRING_, im) ) 
    { ImageIO_free( str ); return -1; }

  std::istringstream iss;
  iss.str(str);
  iss >> im->xdim >> im->ydim >> im->zdim >> im->vdim;

  status = (int)iss.str().length();
  switch ( status ) {
  case 2 :    im->zdim = 1; CGAL_FALLTHROUGH;
  case 3 :    im->vdim = 1; CGAL_FALLTHROUGH;
  case 4 :    break;
  default :
    fprintf( stderr, "readGisHeader: unable to read dimensions in '%s'\n", name );
    ImageIO_free( str ); 
    return -1;
  }
  if ( im->vdim > 1 ) {
    im->vectMode = VM_INTERLACED;
  } 
  else {
    im->vectMode = VM_SCALAR;
  }
#define ADD_USER_STRING { \
    if ( n == 0 ) { \
      im->user = (char**)ImageIO_alloc( nusermax * sizeof( char*) ); \
      for ( n=0; n<nusermax; n++ ) im->user[n] = NULL; \
      n = 0; \
    } \
    im->user[n] = (char*)ImageIO_alloc( 1+strlen( s ) ); \
    strcpy( im->user[n++], s ); \
  }



  while( fgetns( str, _LGTH_STRING_, im ) != 0 ) {
    s = str;
    do {

      while ( *s == ' ' || *s == '\t' ) s++;

      if ( !strncmp( s, "-dx ", 4 ) ) {
	s += 4;
	status = sscanf( s, "%lf", &(im->vx) );
	if ( status != 1 ) {
	  fprintf( stderr, "readGisHeader: error while reading -dx in '%s'\n", s-4 );
	  *s = '\0';
	} 
	else {
	  while ( *s == '.' || (*s >= '0' && *s <= '9') ) s++;
	}
      }
      else if ( !strncmp( s, "-dy ", 4 ) ) {
	s += 4;
	status = sscanf( s, "%lf", &(im->vy) );
	if ( status != 1 ) {
	  fprintf( stderr, "readGisHeader: error while reading -dy in '%s'\n", s-4 );
	  *s = '\0';
	} 
	else {
	  while ( *s == '.' || (*s >= '0' && *s <= '9') ) s++;
	}
      }
      else if ( !strncmp( s, "-dz ", 4 ) ) {
	s += 4;
	status = sscanf( s, "%lf", &(im->vz) );
	if ( status != 1 ) {
	  fprintf( stderr, "readGisHeader: error while reading -dz in '%s'\n", s-4 );
	  *s = '\0';
	} 
	else {
	  while ( *s == '.' || (*s >= '0' && *s <= '9') ) s++;
	}
      }
      else if ( !strncmp( s, "-dt ", 4 ) ) {
	ADD_USER_STRING
	s += 4; 
	while ( *s == '.' || (*s >= '0' && *s <= '9') ) s++; 
      }

      else if ( !strncmp( s, "-type ", 6 ) ) {
	s += 6; 
	if ( !strncmp( s, "U8", 2 ) ) {
	  im->wdim     = 1;
	  im->wordKind = WK_FIXED;
	  im->sign     = SGN_UNSIGNED;
	  s += 2;
	}
	else if ( !strncmp( s, "S8", 2 ) ) {
	  im->wdim     = 1;
	  im->wordKind = WK_FIXED;
	  im->sign     = SGN_SIGNED;
	  s += 2;
	} 
	else if ( !strncmp( s, "U16", 3 ) ) {
	  im->wdim     = 2;
	  im->wordKind = WK_FIXED;
	  im->sign     = SGN_UNSIGNED;
	  s += 3;
	}
	else if ( !strncmp( s, "S16", 3 ) ) {
	  im->wdim     = 2;
	  im->wordKind = WK_FIXED;
	  im->sign     = SGN_SIGNED;
	  s += 3;
	}
	else if ( !strncmp( s, "U32", 3 ) ) {
	  im->wdim     = 4;
	  im->wordKind = WK_FIXED;
	  im->sign     = SGN_UNSIGNED;
	  s += 3;
	}
	else if ( !strncmp( s, "S32", 3 ) ) {
	  im->wdim     = 4;
	  im->wordKind = WK_FIXED;
	  im->sign     = SGN_SIGNED;
	  s += 3;
	}
	else if ( !strncmp( s, "FLOAT", 5 ) ) {
	  im->wdim     = sizeof( float );
	  im->wordKind = WK_FLOAT;
	  im->sign     = SGN_UNKNOWN;
	  s += 5;
	}
	else if ( !strncmp( s, "DOUBLE", 6 ) ) {
	  im->wdim     = sizeof( double );
	  im->wordKind = WK_FLOAT;
	  im->sign     = SGN_UNKNOWN;
	  s += 6;
	}
	else {
	  fprintf( stderr, "readGisHeader: unknown type '%s'\n", s-6 );
	  *s = '\0'; 
	}
      }

      else if ( !strncmp( s, "-bo ", 4 ) ) {
	s += 4;
	if ( !strncmp( s, "ABCD", 4 ) ) {
	  im->endianness = END_BIG;
	  s += 4;
	}
	else if ( !strncmp( s, "SUN", 3 ) ) {
	  im->endianness = END_BIG;
	  s += 3;
	}
	else if ( !strncmp( s, "DCBA", 4 ) ) {
	  im->endianness = END_LITTLE;
	  s += 4;
	}
	else if ( !strncmp( s, "ALPHA", 5 ) ) {
	  im->endianness = END_LITTLE;
	  s += 5;
	}
	else {
	  fprintf( stderr, "readGisHeader: unknown byte order '%s'\n", s-4 );
	  *s = '\0'; 
	}
      }

      else if ( !strncmp( s, "-ar ", 4 ) ) {
	s += 4;
	if ( !strncmp( s, "SUN", 3 ) ) {
	  im->endianness = END_BIG;
	  s += 3;
	}
	else if ( !strncmp( s, "ALPHA", 5 ) ) {
	  im->endianness = END_LITTLE;
	  s += 5;
	}
	else {
	  fprintf( stderr, "readGisHeader: unknown architecture '%s'\n", s-4 );
	  *s = '\0'; 
	}
      }

      else if ( !strncmp( s, "-om ", 4 ) ) {
	s += 4;
	if ( !strncmp( s, "binar", 5 ) ) {
	  im->dataMode = DM_BINARY;
	  s += 5;
	} 
	else if ( !strncmp( s, "ascii", 5 ) ) {
	  im->dataMode = DM_ASCII;
	  s += 5;
	} 
	else {
	  fprintf( stderr, "readGisHeader: unknown data type '%s'\n", s-4 );
	  ImageIO_free( str );
	  return -1;
	}
      }
      
      
      else {
	fprintf( stderr, "readGisHeader: unknown indentifier '%s'\n", s );
	ADD_USER_STRING
	*s = '\0'; 
      }

    } while( *s != '\0' && *s != '\n' );

  }
  ImageIO_free( str );
  

  if ( im->endianness == END_UNKNOWN ) {
    im->endianness = _getEndianness();
  }


  /* header is read. close header file and open data file. */
  if( name != NULL ) {
    
    std::size_t length = strlen(name) ;
    char* data_filename = (char *) ImageIO_alloc(length+4) ;
    
    if( strcmp( name+length-4, ".dim" ) ) {
      fprintf (stderr,
	       "readGisHeader: error: file header extension must be .dim\n");
      ImageIO_free( data_filename );
      return -1;
    }
    
    ImageIO_close(im);


    /* open data file
     */
    
    strcpy(data_filename,name);
    strcpy(data_filename+length-3, "ima.gz");
    _openReadImage(im,data_filename);
    
    if(!im->fd) {
      
      strcpy(data_filename,name);
      strcpy(data_filename+length-3, "ima");
      _openReadImage(im,data_filename);
      if(!im->fd) {
	fprintf(stderr, "readGisHeader: error: unable to open data file \'%s\'\n", data_filename);
	ImageIO_free( data_filename );
	return -1;
	
      }
    }
    ImageIO_free( data_filename );





    /* read data if ascii
       only U8 and S8
     */
    if ( im->dataMode == DM_ASCII ) {
      std::size_t size = std::size_t(im->xdim) * im->ydim * im->zdim * im->vdim * im->wdim;
      unsigned int n;
      char *tmp;
      int ret, iv=0;
      
      if ( im->wdim != 1 || im->wordKind != WK_FIXED ) {
	fprintf(stderr, "readGisHeader: error: unable to read such ascii type\n" );
	return -1;
      }

      n = 0 ;
      if ( size <= 0 ) return -1;
      if(!im->data) {
	im->data = ( void*) ImageIO_alloc(size);
	if(!im->data) return -1;
      }
      
      n = 0;
      str = (char*)ImageIO_alloc( _LGTH_STRING_+1 );
      while( fgetns( str, _LGTH_STRING_, im ) != 0 &&
	     n < im->xdim * im->ydim * im->zdim * im->vdim ) {
	tmp = str;
	while ( *tmp != '\n' && *tmp != '\0' && *tmp != EOF &&
		n < im->xdim * im->ydim * im->zdim * im->vdim ) {
	  /* skip trailing whitespace
	   */
	  while ( *tmp == ' ' || *tmp == '\t' )
	    tmp++;
	  if ( *tmp == '\0' || *tmp == '\n' || *tmp == EOF )
	    continue;
	  
	  /* read a number
	   */
	  switch ( im->wordKind ) {
	  case WK_FIXED :
	    ret = sscanf( tmp, "%d", &iv );
	  break;
	  default :
	    ImageIO_free( im->data ); im->data = NULL;
	    ImageIO_free( str );
	    return -1;
	  }
	  
	  if ( ret != 1 ) {
	    fprintf( stderr, "readGisHeader: error in reading ascii data\n" );
	    ImageIO_free( im->data ); im->data = NULL;
	    ImageIO_free( str );
	    return -1;	
	  }
	  
	  if ( im->wordKind == WK_FIXED 
	       && im->sign == SGN_UNSIGNED 
	       && im->wdim == 1 ) {
	    unsigned char *buf = (unsigned char *)im->data;
	    buf += n;
	    if ( iv < 0 )        *buf = (unsigned char)0;
	    else if ( iv > 255 ) *buf = (unsigned char)255;
	    else                 *buf = (unsigned char)iv;
	    n ++;
	  }
	  else if ( im->wordKind == WK_FIXED 
		    && im->sign == SGN_SIGNED 
		    && im->wdim == 1 ) {
	    char *buf = (char *)im->data;
	    buf += n;
	    if ( iv < -128 )     *buf = (char)-128;
	    else if ( iv > 127 ) *buf = (char)127;
	    else                 *buf = (char)iv;
	    n ++;
	  }
	  else if ( im->wordKind == WK_FIXED 
		    && im->sign == SGN_UNSIGNED 
		    && im->wdim == 2 ) {
	    unsigned short int *buf = (unsigned short int *)im->data;
	    buf += n;
	    if ( iv < 0 )          *buf = (unsigned short int)0;
	    else if ( iv > 65535 ) *buf = (unsigned short int)65535;
	    else                   *buf = (unsigned short int)iv;
	    n ++;
	  }
	  else if ( im->wordKind == WK_FIXED 
		    && im->sign == SGN_SIGNED 
		    && im->wdim == 2 ) {
	    short int *buf = (short int *)im->data;
	    buf += n;
	    if ( iv < -32768 )     *buf = (short int)-32768;
	    else if ( iv > 32767 ) *buf = (short int)32767;
	    else                   *buf = (short int)iv;
	    n ++;
	  }
	  else {
	    ImageIO_free( im->data ); im->data = NULL;
	    ImageIO_free( str );
	    return -1;
	  }
	  


	  /* skip a number 
	   */
	  while ( (*tmp >= '0' && *tmp <= '9') || *tmp == '.' || *tmp == '-' )
	    tmp++;
	}
      }
      ImageIO_free( str );
      ImageIO_close(im);
    }

  }

  
  /* check header validity */
  if ( im->xdim > 0 && im->ydim > 0 && im->zdim > 0 && im->vdim > 0 &&
       im->vx > 0.0 && im->vy > 0.0 && im->vz > 0.0 &&
       ( im->wordKind == WK_FLOAT ||
	 (im->wordKind == WK_FIXED && im->sign != SGN_UNKNOWN) ) &&
       im->endianness != END_UNKNOWN ) {
    return 0;
  }

  return -1;
}



CGAL_INLINE_FUNCTION
bool writeGisHeader( const _image* inr )
{
  const char *proc = "writeGisHeader";
  std::ostringstream oss;

  if ( inr->vectMode == VM_NON_INTERLACED ) {
    fprintf( stderr, "%s: can not write non interlaced data\n", proc );
    return false;
  }

  /* dimensions
   */
  oss << " " << inr->xdim << " " << inr->ydim;
  if ( inr->vdim > 1 ) {
    oss << " " << inr->zdim << " " << inr->vdim;
  }
  else if ( inr->zdim > 1 ) {
    oss << " " << inr->zdim;
  }
  oss << "\n";

  /* type
   */
  oss << "-type ";
  switch ( inr->wordKind ) {
  case WK_FIXED :
    switch( inr->sign ) {
    case SGN_UNSIGNED :
      oss << "U" << 8*inr->wdim;
      break;
    case SGN_SIGNED :
      oss << "S" << 8*inr->wdim;
      break;
    default :
      fprintf( stderr, "%s: unknown wordSign\n", proc );
      return false;    
    }
    break;
  case WK_FLOAT :
    if ( inr->wdim == sizeof( float ) ) {
      oss << "FLOAT";
    }
    else if ( inr->wdim  == sizeof( double ) ) {
      oss << "DOUBLE";
    }
    else {
      fprintf( stderr, "%s: unknown WK_FLOAT word dim\n", proc );
      return false;    
    }
    break;
  default :
    fprintf( stderr, "%s: unknown wordKind for image\n", proc );
    return false;  
  }
  oss << "\n";
  
  oss << "-dx "<< inr->vx <<"\n";
  oss << "-dy "<< inr->vy <<"\n";
  if ( inr->zdim > 1 )
    oss << "-dz " << inr->vz << "\n";

  if ( inr->wdim > 1 ) {
    oss << "-bo ";
    switch ( _getEndianness() ) {
    default :
    case END_LITTLE :
      oss << "DCBA"; break;
    case END_BIG :
      oss << "ABCD"; break;
    }
    oss << "\n";
  }
  switch ( inr->dataMode ) {
  default :
  case DM_BINARY :
    oss << "-om binar\n";
    break;
  case DM_ASCII :
    oss << "-om ascii\n";
  }
  if( ImageIO_write( inr, oss.str().data(), oss.str().length()) == 0) {
    return false;
  }
  return true;
}


CGAL_INLINE_FUNCTION
int writeGisData( const _image*) 
{
  return -1;
}
