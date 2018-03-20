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

#ifndef PNM_H
#define PNM_H

#include <stdio.h>

#include <CGAL/ImageIO.h>


int readPpmImage(const char *name,_image *im);
int writePpmImage(char *name, _image *im);
int readPgmAsciiImage(const char *name,_image *im);
int readPgmImage(const char *name,_image *im);
int writePgmImage(char *name,  _image *im);
int testPgmAsciiHeader(char *magic,const char *name);
int testPgmHeader(char *magic,const char *name);
int testPpmHeader(char *magic,const char *name);
PTRIMAGE_FORMAT createPgmFormat();
PTRIMAGE_FORMAT createPgmAscIIFormat();
PTRIMAGE_FORMAT createPpmFormat();

#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/pnm_impl.h>
#endif // CGAL_HEADER_ONLY

#endif
