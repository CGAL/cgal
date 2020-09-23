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

#ifndef GIF_H
#define GIF_H


#include <stdio.h>
#include <CGAL/ImageIO.h>
int readGifImage(const char *name,_image *im);

int testGifHeader(char *magic,const char *name);
/** creates an return the file format structure associated with the Gif file format */
PTRIMAGE_FORMAT createGifFormat();

#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/gif_impl.h>
#endif // CGAL_HEADER_ONLY

#endif
