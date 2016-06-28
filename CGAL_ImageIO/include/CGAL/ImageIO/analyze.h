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

#ifndef ANALYZE_H
#define ANALYZE_H

#ifdef _MSC_VER
#pragma warning ( disable : 4068 4786 4081 4267 )
#endif



#include <CGAL/ImageIO.h>

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

#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/analyze_impl.h>
#endif // CGAL_HEADER_ONLY

#endif
