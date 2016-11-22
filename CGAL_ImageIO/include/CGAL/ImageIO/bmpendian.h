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

/*
 * from bmp.zip, see the url http://www.ddj.com/ftp/1995/1995.03/
 * author Dr. Dobb's
 */

/*
 * This is the header for endian.c - functions to read/write our
 * CGAL_INT8, CGAL_INT16 and CGAL_INT32 types from/to a little-endian file.
 */

#ifndef __ENDIAN_H_INCLUDED__
#define __ENDIAN_H_INCLUDED__

#include <cstdio>
#include <CGAL/ImageIO/bmptypes.h>
#include <cstdio>

/*
 * Read the basic types as little-endian values.  The return value will be
 * zero if successful, EOF, otherwise.
 */
extern int readINT8little(FILE *f, CGAL_INT8 *i);
extern int readINT16little(FILE *f, CGAL_INT16 *i);
extern int readINT32little(FILE *f, CGAL_INT32 *i);
extern int readUINT8little(FILE *f, CGAL_UINT8 *i);
extern int readUINT16little(FILE *f, CGAL_UINT16 *i);
extern int readUINT32little(FILE *f, CGAL_UINT32 *i);

/*
 * Write the basic types as little-endian values.  The return value will be
 * zero if successful, EOF otherwise.
 */
int writeINT8little(FILE *f, CGAL_INT8 i);
int writeINT16little(FILE *f, CGAL_INT16 i);
int writeINT32little(FILE *f, CGAL_INT32 i);
int writeUINT8little(FILE *f, CGAL_UINT8 i);
int writeUINT16little(FILE *f, CGAL_UINT16 i);
int writeUINT32little(FILE *f, CGAL_UINT32 i);

#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/bmpendian_impl.h>
#endif // CGAL_HEADER_ONLY

#endif  /* __ENDIAN_H_INCLUDED__ */


/*
 * Formatting information for emacs in c-mode
 *
 * Local Variables:
 * c-indent-level:4
 * c-continued-statement-offset:4
 * c-brace-offset:-4
 * c-brace-imaginary-offset:0
 * c-argdecl-indent:4
 * c-label-offset:-4
 * End:
 */

