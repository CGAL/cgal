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
 * convert.h - conversion between types
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
 * ADDITIONS, CHANGES
 *
 */

#ifndef _convert_h_
#define _convert_h_







#include <stdio.h> /* memcpy */
#include <string.h> /* memcpy */
#include <CGAL/ImageIO/typedefs.h>






/* Converts buffers from one type into another.
 *
 * DESCRIPTION:
 *
 * Conversions from types DOUBLE or FLOAT to
 * others types are all implemented. I hope
 * I do not make mistakes.
 *
 * - When a real value is converted into an integer
 *   one, we take the nearest integer, i.e.
 *   (int)(real value - 0.5) for negative
 *   numbers, and (int)(real value + 0.5) for
 *   positive ones.
 *
 * - Some conversions are not yet implemented.
 *
 * - When types and pointers are equal, nothing is
 *   done.
 */
extern void ConvertBuffer( void *bufferIn, /* buffer to be converted */
                           bufferType typeIn, /* type of this buffer */
                           void *bufferOut, /* result buffer */
                           bufferType typeOut, /* type of this buffer */
                           int bufferLength /* buffers' length */
                           );

extern void Convert_r32_to_s8( r32 *theBuf,
                               s8 *resBuf,
                               int size );

extern void Convert_r32_to_u8( r32 *theBuf,
                               u8 *resBuf,
                               int size );

extern void Convert_r32_to_s16( r32 *theBuf,
                                s16 *resBuf,
                                int size );

extern void Convert_r32_to_u16( r32 *theBuf,
                                u16 *resBuf,
                                int size );





#ifdef CGAL_HEADER_ONLY
#include <CGAL/ImageIO/convert_impl.h>
#endif // CGAL_HEADER_ONLY

#endif
