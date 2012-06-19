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

#include "gis.h" 
#include "inr.h"

/* get a string from a file and discard the ending newline character
   if any */
char *fgetns(char *str, int n,  _image *im ) {
  char *ret;
  int l;

  memset( str, 0, n );
  ret = ImageIO_gets( im, str, n );

  if(!ret) return NULL;

  l = strlen(str);
  if(l > 0 && str[l-1] == '\n') str[l-1] = '\0';
  return ret;
}
