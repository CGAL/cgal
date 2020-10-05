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


#ifndef INR_H
#define INR_H


#include <CGAL/ImageIO.h>



/** Writes the given inrimage header in an already opened file.
    @param im image descriptor
    @param end image data endianness (END_UNKNOWN to use architecture endianness) */
int _writeInrimageHeader(const _image *im,
                         ENDIANNESS end);

/** Writes the given image body in an already opened file.
    @param im image descriptor */
bool _writeInrimageData(const _image *im);

/** read header from an opened inrimage file
    @param im opened inrmage descriptor */
int readInrimageHeader(const char *,_image *im);

/** test if an image file has the Inrimage format, given the magic string
    and filename
    @param magic this string is of size 5 and corresponds to the first 4 characters of the file (terminated with a null character)
    @param filename the file name
    @return -1 if it is not an Inrimage file and 0 otherwise */
int testInrimageHeader(char *magic,const char *filename);

/** calls _writeInrimageHeader and _writeInrimageData
    @param basename the file name without any extension
    @param im structure
    @return same as  _writeInrimageHeader*/
int writeInrimage(char *basename,_image *im);
/** creates an return the file format structure associated with the Inrimage file format */
PTRIMAGE_FORMAT createInrimageFormat();

#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/inr_impl.h>
#endif // CGAL_HEADER_ONLY


#endif
