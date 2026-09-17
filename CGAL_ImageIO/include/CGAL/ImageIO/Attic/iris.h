// Copyright (c) 2005-2008 ASCLEPIOS Project, INRIA Sophia-Antipolis (France)
// All rights reserved.
//
// This file is part of the ImageIO Library, and as been adapted for CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later
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
