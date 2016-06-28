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

/*************************************************************************
 * iobmp.h - I procedures for BMP raw images
 *
 * $Id$
 *
 * Copyright INRIA
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * 
 * CREATION DATE: 
 * Wed Oct  6 17:03:48 MET DST 1999
 *
 * ADDITIONS, CHANGES
 *
 *
 */

#ifndef _bmp_h_
#define _bmp_h_



#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h> /* open, close */
#include <sys/stat.h> /* open, close */
#include <sys/types.h> /* open, close */
#include <string.h>
#include <CGAL/ImageIO.h>
extern int readBmpImage(const char *name, _image *im);
extern void *_readBmpImage( const char *name, int *dimx, int *dimy, int *dimz );
int testBmpHeader(char *magic,const char *name);
/** creates an return the file format structure associated with the BMP file format */
PTRIMAGE_FORMAT createBMPFormat();

extern void IoBmp_verbose ( );
extern void IoBmp_noverbose ( );

#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/bmp_impl.h>
#endif // CGAL_HEADER_ONLY

#endif /* _bmp_h_ */
