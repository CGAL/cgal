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

/*************************************************************************
 * typedefs.h -
 *
 * $Id$
 *
 * Copyright©INRIA 1998
 *
 * AUTHOR:
 * Gregoire Malandain (greg@sophia.inria.fr)
 * http://www.inria.fr/epidaure/personnel/malandain/
 *
 * CREATION DATE:
 * June, 9 1998
 *
 *
 *
 */

#ifndef _typedefs_h_
#define _typedefs_h_





/* Differents type coding for images and buffers.
 */
typedef enum {
  TYPE_UNKNOWN /* unknown type */,
  CGAL_UCHAR  /* unsigned char */,
  CGAL_SCHAR  /* signed char */,
  CGAL_USHORT /* unsigned short int */,
  CGAL_SSHORT /* signed short int */,
  CGAL_INT    /* signed int */,
  CGAL_ULINT  /* unsigned long int */,
  CGAL_FLOAT  /* float */,
  CGAL_DOUBLE  /* double */
} ImageType, bufferType;

typedef char               s8;
typedef unsigned char      u8;
typedef short int          s16;
typedef unsigned short int u16;
typedef int                i32;
typedef int                s32;
typedef unsigned long int  u64;
typedef float              r32;
typedef double             r64;





/* Typedef Booleen
 */
typedef enum {
  False = 0,
  True = 1
} typeBoolean;



#endif
