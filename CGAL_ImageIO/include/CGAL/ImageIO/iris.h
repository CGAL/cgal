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

#ifndef IRIS_H
#define IRIS_H

#include <CGAL/ImageIO.h>



PTRIMAGE_FORMAT createIrisFormat();
int readIrisImage( const char *name, _image *im );
int testIrisHeader(char *magic,const char *name);

#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/iris_impl.h>
#endif // CGAL_HEADER_ONLY

#endif
