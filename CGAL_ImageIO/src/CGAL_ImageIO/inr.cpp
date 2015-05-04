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

#include "inr.h"
#include "fgetns.h"

#include <string>
#include <sstream>

/* Magic header for inrimages v4 */
#define INR4_MAGIC "#INRIMAGE-4#{"


/** Magic header for inrimages */
#define INR_MAGIC "#INR"

typedef struct stringListElementStruct {
  char *string;
  struct stringListElementStruct *next;
} stringListElement;
/* list element with a pointer on a string */

typedef struct {
  stringListElement *begin, *end;
} stringListHead;
/* string list descriptor */

#include <clocale>
class Set_numeric_locale {
  const char * old_locale;
public:
  Set_numeric_locale(const char* locale) 
    : old_locale(std::setlocale(LC_NUMERIC, locale))
  {
  }

  ~Set_numeric_locale() {
    std::setlocale(LC_NUMERIC, old_locale);
  }
};

static void addStringElement(stringListHead *strhead,
			     const char *str);
/* add a string element at the tail of given list */

static void concatStringElement(const stringListHead *strhead,
				const char *str);
/* concat given string at the last element of given list */


/* Writes the given inrimage header in an already opened file.*/
int _writeInrimageHeader(const _image *im, ENDIANNESS end) {
  unsigned int pos, i;
  char type[30], endianness[5], buf[257], scale[20];
  std::ostringstream oss;

  Set_numeric_locale num_locale("C");

  if(im->openMode != OM_CLOSE) {
    /* fix word kind */
    switch(im->wordKind) {

    case WK_FLOAT:
      sprintf(type, "float");
      scale[0] = '\0';
      break;

    case WK_FIXED:
      switch(im->sign) {
      case SGN_SIGNED:
	sprintf(type, "signed fixed");
	break;
      case SGN_UNSIGNED:
	sprintf(type, "unsigned fixed");
	break;
      default:
	return -1;
      }
      sprintf(scale, "SCALE=2**0\n");
      break;
      
    default:
      return -1;
    }
    
    switch(end) {
    case END_LITTLE:
      sprintf(endianness, "decm");
      break;
    case END_BIG:
      sprintf(endianness, "sun");
      break;
    default:
      /* fix architecture endianness */
      if( _getEndianness() == END_LITTLE)
	sprintf(endianness, "decm");
      else
	sprintf(endianness, "sun");
      break;
    }

    /* write header information */
    oss << INR4_MAGIC << "\n";
    oss << "XDIM=" << im->xdim << "\n";
    oss << "YDIM=" << im->ydim << "\n";
    oss << "ZDIM=" << im->zdim << "\n";
    oss << "VDIM=" << im->vdim << "\n";
    oss << "TYPE=" << type << "\n";
    oss << "PIXSIZE=" << im->wdim*8 <<" bits\n";
    oss << scale << "CPU=" << endianness << "\n";
    oss << "VX=" << im->vx << "\n";
    oss << "VY=" << im->vy << "\n";
    oss << "VZ=" << im->vz << "\n";

    if ( im->cx != 0 ) 
      oss << "XO="<< im->cx << "\n";
    if ( im->cy != 0 ) 
      oss << "YO="<< im->cy << "\n";
    if ( im->cz != 0 ) 
      oss << "ZO="<< im->cz << "\n";
    if ( im->tx != 0.0 ) 
      oss << "TX="<< im->tx << "\n";
    if ( im->ty != 0.0 ) 
      oss << "TY="<< im->ty << "\n";
    if ( im->tz != 0.0 ) 
      oss << "TZ="<< im->tz << "\n";
    if ( im->rx != 0.0 ) 
      oss << "RX="<< im->rx <<"\n";
    if ( im->ry != 0.0 ) 
      oss << "RY="<< im->ry << "\n";
    if ( im->rz != 0.0 ) 
      oss << "RZ=" << im->rz <<"\n";

    pos = oss.str().length();
    
    if(ImageIO_write(im, oss.str().data(), oss.str().length()) == 0)
      return -1;
    
    
    /* write user strings */
    if ( im->user != NULL ) {
      for(i = 0; i < im->nuser; i++) {
	if ( im->user[i] == NULL ) continue;
	pos += strlen(im->user[i]) + 2;
	if(ImageIO_write(im, "#", 1) == 0) return -1;
	if(ImageIO_write(im, im->user[i], strlen(im->user[i])) == 0) return -1;
	if(ImageIO_write(im, "\n", 1) == 0) return -1;
      }
    }
    /* write end of header */
    pos = pos % 256;
    if(pos > 252) {
      for(i = pos; i < 256; i++)
	if(ImageIO_write(im, "\n", 1) != 1) return -1;
      pos = 0;
    }
    buf[0] = '\0';
    for(i = pos; i < 252; i++) strcat(buf, "\n");
    strcat(buf, "##}\n");
    
    if(ImageIO_write(im, buf, strlen(buf)) == 0) return -1;
    else return 1;
  }

  else return -1;
}



/* Writes the given image body in an already opened file.*/
int _writeInrimageData(const _image *im) {
  unsigned long size, nbv, nwrt, i;
  unsigned int v;
  unsigned char **vp;
  
  if(im->openMode != OM_CLOSE) {

    /* scalar or interlaced vectors */
    if(im->vectMode != VM_NON_INTERLACED) {
      size = im->xdim * im->ydim * im->zdim * im->vdim * im->wdim;
      nwrt = ImageIO_write(im, im->data, size);
      if(nwrt != size) return -1;
      else return 1;
    }

    /* non interlaced vectors: interlace for saving */
    else {
      nbv = im->xdim * im->ydim * im->zdim;
      size = im->xdim * im->ydim * im->zdim * im->wdim;
      vp = (unsigned char **) ImageIO_alloc(im->vdim * sizeof(unsigned char *));
      for(v = 0; v < im->vdim; v++)
	vp[v] = (unsigned char *) im->data + v * size;
      for(i = 0; i < nbv; i++)
	for(v = 0; v < im->vdim; v++) {
	  if(ImageIO_write(im, (const void *) vp[v], im->wdim) != im->wdim)
	    return -1;
	  vp[v] += im->wdim;
	}
      ImageIO_free(vp);
      return 1;
    }
  }
  else return -1;
}




/* read header of an opened inrimage */
int readInrimageHeader(const char *,_image *im) {
  char str[257];
  int n, nusr;
  stringListHead strl = { NULL, NULL };
  stringListElement *oel, *el;

  Set_numeric_locale num_locale("C");

  if(im->openMode != OM_CLOSE) {
    /* read image magic number */
    if(!fgetns(str, 257, im )) return -1;
    if(strcmp(str, INR4_MAGIC)) return -1;


    /* while read line does not begin with '#' or '\n', read line
       and decode field */
    if(!fgetns(str, 257, im)) return -1;

    while(str[0] != '#' && str[0] != '\0') {

      if(!strncmp(str, "XDIM=", 5)) {
        std::istringstream iss(str+5);
        if(!(iss >> im->xdim)) return -1;
      }
      else if(!strncmp(str, "YDIM=", 5)) {
        std::istringstream iss(str+5);
        if(!(iss >> im->ydim)) return -1;
      }
      else if(!strncmp(str, "ZDIM=", 5)) {
        std::istringstream iss(str+5);
        if(!(iss >> im->zdim)) return -1;
      }
      else if(!strncmp(str, "VDIM=", 5)) {
        std::istringstream iss(str+5);
        if(!(iss >> im->vdim)) return -1;
	if(im->vdim == 1) im->vectMode = VM_SCALAR;
	else im->vectMode = VM_INTERLACED;
      }
      else if(!strncmp(str, "VX=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->vx)) return -1;
      }
      else if(!strncmp(str, "VY=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->vy)) return -1;
      }
      else if(!strncmp(str, "VZ=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->vz)) return -1;
      }
      else if(!strncmp(str, "TYPE=", 5)) {
	if(!strncmp(str+5, "float", 5)) im->wordKind = WK_FLOAT;
	else {
	  if(!strncmp(str+5, "signed fixed", 12)) {
	    im->wordKind = WK_FIXED;
	    im->sign = SGN_SIGNED;
	  }
	  else if(!strncmp(str+5, "unsigned fixed", 14)) {
	    im->wordKind = WK_FIXED;
	    im->sign = SGN_UNSIGNED;
	  }
	  else return -1;
	}
      }
      /* before "sscanf(str+8, "%i %n", &im->wdim, &n) != 1"
	 was used. 
	 However the man said 
         ...
	 n      Nothing is expected; instead, the number of charac­
              ters consumed thus far from  the  input  is  stored
              through  the  next pointer, which must be a pointer
              to int.  This is not a conversion, although it  can
              be  suppressed  with  the  *  flag.  The C standard
              says: `Execution of a %n directive does not  incre­
              ment  the  assignment count returned at the comple­
              tion of execution' but  the  Corrigendum  seems  to
              contradict  this.  Probably  it is wise not to make
              any assumptions on the effect of %n conversions  on
              the return value.
	 ...
	 Thus I change it. It was yielding a RETURN_FAILURE with 
	 insight (GM).
      */
      else if(!strncmp(str, "PIXSIZE=", 8)) {
        std::istringstream iss(str+8);
        if(!(iss >> im->wdim)) return -1;
	if(im->wdim != 8 && im->wdim != 16 && im->wdim != 32 &&
	   im->wdim != 64) return -1;
	
	if ( im->wdim <= 9 ) {
	  if(strncmp(str+8+1, " bits", 5)) return -1;
	}
	else if ( im->wdim <= 99 ) {
	  if(strncmp(str+8+2, " bits", 5)) return -1;
	}
	else {
	  return -1;
	}

	im->wdim >>= 3;
      }
      else if(!strncmp(str, "SCALE=", 6)) ;
      else if(!strncmp(str, "CPU=", 4)) {
	if(!strncmp(str+4, "decm", 4)) im->endianness = END_LITTLE;
	else if(!strncmp(str+4, "alpha", 5)) im->endianness = END_LITTLE;
	else if(!strncmp(str+4, "pc", 2)) im->endianness = END_LITTLE;
	else if(!strncmp(str+4, "sun", 3)) im->endianness = END_BIG;
	else if(!strncmp(str+4, "sgi", 3)) im->endianness = END_BIG;
	else return -1;
      }

      else if(!strncmp(str, "XO=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->cx)) return -1;
      }
      else if(!strncmp(str, "YO=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->cy)) return -1;
      }
      else if(!strncmp(str, "ZO=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->cz)) return -1;
      }

      else if(!strncmp(str, "TX=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->tx)) return -1;
      }
      else if(!strncmp(str, "TY=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->ty)) return -1;
      }
      else if(!strncmp(str, "TZ=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->tz)) return -1;
      }
      else if(!strncmp(str, "RX=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->rx)) return -1;
      }
      else if(!strncmp(str, "RY=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->ry)) return -1;
      }
      else if(!strncmp(str, "RZ=", 3)) {
        std::istringstream iss(str+3);
        if(!(iss >> im->rz)) return -1;
      }

      if(!fgetns(str, 257, im)) return -1;
    }

    /* parse user strings */
    im->nuser = nusr = 0;
    while(str[0] == '#' && strncmp(str, "##}", 3)) {
      addStringElement(&strl, str + 1);
      while(strlen(str) == 256) {
	if(!fgetns(str, 257, im)) return -1;
	concatStringElement(&strl, str);
      }
      nusr++;
      if(!fgetns(str, 257, im)) return -1;      
    }
    
    /* go to end of header */
    while(strncmp(str, "##}", 3)) {
      if(!fgetns(str, 257, im)) return -1;
    }
    

    /* check header validity */
    if(im->xdim > 0 && im->ydim > 0 && im->zdim > 0 && im->vdim > 0 &&
       im->vx > 0.0 && im->vy > 0.0 && im->vz > 0.0 &&
       (im->wordKind == WK_FLOAT || (im->wordKind == WK_FIXED &&
				     im->sign != SGN_UNKNOWN)) &&
       im->endianness != END_UNKNOWN) {
      if(nusr > 0) {
	im->nuser = nusr;
	im->user = (char **) ImageIO_alloc(im->nuser * sizeof(char *));
	oel = NULL;
	for(el = strl.begin, n = 0; el != NULL; el = oel, n++) {
	  im->user[n] = el->string;
	  oel = el->next;
	  ImageIO_free(el);
	}
      }
      return 0;
    }
    else return -1;

  }
  else return -1;
}



/* add a string element at the tail of given list */
static void addStringElement(stringListHead *strhead, const char *str) {
  stringListElement *el;

  el = (stringListElement *) ImageIO_alloc(sizeof(stringListElement));
  /* was strdup(str); */
  el->string = (char*)ImageIO_alloc( strlen(str)+1);
  memcpy(el->string, str,  strlen(str)+1);
  el->next = NULL;
  if(strhead->begin == NULL)
    strhead->begin = strhead->end = el;
  else {
    strhead->end->next = el;
    strhead->end = el;
  }
}


/* concat given string at the last element of given list */
static void concatStringElement(const stringListHead *strhead,
				const char *str) {
  stringListElement *el;

  el = strhead->end;
  el->string = (char *) realloc(el->string,
				strlen(el->string) + strlen(str) + 1);
  strcat(el->string, str);
}

int testInrimageHeader(char *magic,const char *) {
  if (!strcmp(magic, INR_MAGIC))
    return 0;
  else 
    return -1;
}

int writeInrimage(char *name,_image *im) {
  int res;

  _openWriteImage( im, name );

  if(!im->fd) {
    fprintf(stderr, "writeInrimage: error: unable to open file \'%s\'\n", name );
    return ImageIO_OPENING;
  }

  res = _writeInrimageHeader(im, END_UNKNOWN);
  if (res < 0) {
    fprintf(stderr, "writeInrimage: error: unable to write header of \'%s\'\n",
	    name);
    ImageIO_close( im );
    im->fd = NULL;
    im->openMode = OM_CLOSE;
    return( res );
  }
  
  res = _writeInrimageData( im );
  if (res < 0) {
    fprintf(stderr, "writeInrimage: error: unable to write data of \'%s\'\n",
	    name);
    ImageIO_close( im );
    im->fd = NULL;
    im->openMode = OM_CLOSE;
    return( res );
  }

  ImageIO_close( im );
  im->fd = NULL;
  im->openMode = OM_CLOSE;

  return ( res );  
}

PTRIMAGE_FORMAT createInrimageFormat() {
  PTRIMAGE_FORMAT f=(PTRIMAGE_FORMAT) ImageIO_alloc(sizeof(IMAGE_FORMAT));

  f->testImageFormat=&testInrimageHeader;
  f->readImageHeader=&readInrimageHeader;
  f->writeImage=&writeInrimage;
  strcpy(f->fileExtension,".inr,.inr.gz,.gradient,.gradient.gz,.gradient_direction,.gradient_direction.gz");
  strcpy(f->realName,"Inrimage");
  return f;
}
