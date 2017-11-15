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

/*
 * from bmp.zip, see the url http://www.ddj.com/ftp/1995/1995.03/
 * author Dr. Dobb's
 */

/*
 * Data types used in bitmap files.
 */

#ifndef __BMPTYPES_H_INCLUDED__
#define __BMPTYPES_H_INCLUDED__

/*****************************************************************************
*
* Data types.
*
* CGAL_INT8 is an integer of at least 8 bits wide.
* CGAL_INT16 is an integer of at least 16 bits wide.
* CGAL_INT32 is an integer of at least 32 bits wide.
*
* CGAL_UINT8 is an unsigned CGAL_INT8
* CGAL_UINT16 is an unsigned CGAL_INT16
* CGAL_UINT32 is an unsigned CGAL_INT32
*/

#include <boost/cstdint.hpp>

typedef char            CGAL_INT8;
typedef short           CGAL_INT16;
typedef boost::int32_t  CGAL_INT32;
typedef unsigned char   CGAL_UINT8;
typedef unsigned short  CGAL_UINT16;
typedef boost::uint32_t CGAL_UINT32;

/*****************************************************************************
*
* Constants.  Each set corresponds to a field in a structure.  Note that some
* sets only have one value.  Default values for all fields are the value
* corresponding to 0. 
*/

/*
 * Constants used in the "type" field of Bitmapfileheader and
 * BITMAPARRAYHEADER structures.  Note that these are all two-character
 * mnemonics as well as integer constants.
 */
#define TYPE_ICO        (0x4349)   /* 'IC' */
#define TYPE_BMP        (0x4d42)   /* 'BM' */
#define TYPE_PTR        (0x5450)   /* 'PT' */
#define TYPE_ICO_COLOR  (0x4943)   /* 'CI' */
#define TYPE_PTR_COLOR  (0x5043)   /* 'CP' */
#define TYPE_ARRAY      (0x4142)   /* 'BA' */
   
/*
 * Compression schemes.  Note that BITFIELDS (from NT) uses the same number as
 * HUFFMAN1D (from OS/2)
 */
#define COMPRESSION_NONE       (0)
#define COMPRESSION_RLE_8      (1)
#define COMPRESSION_RLE_4      (2)
#define COMPRESSION_HUFFMAN1D  (3)
#define COMPRESSION_BITFIELDS  (3)
#define COMPRESSION_RLE_24     (4)
#define COMPRESSION_LAST       (4)
   
/*
 * units of resolution
 */
#define UNITS_PELS_PER_METER (0)
#define UNITS_LAST           (0)

/*
 * origin of coordinate space
 */   
#define ORIGIN_LOWER_LEFT  (0)
#define ORIGIN_LAST        (0)

/*
 * halftoning algorithms
 */   
#define HALFTONING_NONE             (0)
#define HALFTONING_ERROR_DIFFUSION  (1)
#define HALFTONING_PANDA            (2)
#define HALFTONING_SUPER_CIRCLE     (3)
#define HALFTONING_LAST             (3)
   
/*
 * color table encoding
 */
#define COLOR_ENCODING_RGB   (0)
#define COLOR_ENCODING_LAST  (0)

/*****************************************************************************
*
* Structures.
*/
   
/*
 * Bitmapfileheader defines a single bitmap image.  Its analogue in the
 * Windows SDK is the Bitmapfileheader.  Its analogues in the OS/2 Toolkit are
 * the Bitmapfileheader and Bitmapfileheader2 structures.
 *
 * A BITMAPHEADER structure is always concatenated to the end of a
 * Bitmapfileheader structure.
 */
typedef struct Bitmapfileheader
{
    CGAL_UINT16    type;
    CGAL_UINT32    size;
    CGAL_INT16     xHotspot;
    CGAL_INT16     yHotspot;
    CGAL_UINT32    offsetToBits;
} Bitmapfileheader;


/*
 * BITMAPARRAYHEADER is used to establish a linked list of Bitmapfileheader
 * structures for a bitmap file with multiple images in it.  There is no
 * equivalent structure in the Windows SDK.  Its analogues in the OS/2 toolkit
 * are the BITMAPARRAYFILEHEADER and BITMAPARRAYFILEHEADER2 strucutres.
 *
 * A Bitmapfileheader structure is always concatenated to the end of a
 * BITMAPARRAYHEADER structure.
 */
typedef struct BITMAPARRAYHEADER
{
    CGAL_UINT16    type;
    CGAL_UINT32    size;
    CGAL_UINT32    next;
    CGAL_UINT16    screenWidth;
    CGAL_UINT16    screenHeight;
} BITMAPARRAYHEADER;
   

/*
 * BITMAPHEADER defines the properties of a bitmap.  Its analogues in the
 * Windows SDK are the BITMAPCOREINFOHEADER and BITMAPINFOHEADER structures.
 * Its analogues in the OS/2 Toolkit are the BITMAPINFOHEADER and
 * BITMAPINFOHEADER2 structures.
 *
 * A color table is concatenated to this structure.  The number of elements in
 * the color table determined by the bit-depth of the image.
 *
 * Note, that if the field "size" is 12 or less, then the width and height
 * fields should be read as CGAL_UINT16's instead of CGAL_UINT32's.
 *
 * Also note that if the field "size" is greater than 12, then the color table
 * will have an extra byte of padding between each structures (to longword
 * align it)
 *
 * The different sizes for the width, height, and color table are the only
 * differences between the "old" and "new" bitmap file formats.
 */
typedef struct BITMAPHEADER
{
    CGAL_UINT32 size;
    CGAL_INT32  width;
    CGAL_INT32  height;
    CGAL_UINT16 numBitPlanes;
    CGAL_UINT16 numBitsPerPlane;
    CGAL_UINT32 compressionScheme;
    CGAL_UINT32 sizeOfImageData;
    CGAL_UINT32 xResolution;
    CGAL_UINT32 yResolution;
    CGAL_UINT32 numColorsUsed;
    CGAL_UINT32 numImportantColors;
    CGAL_UINT16 resolutionUnits;
    CGAL_UINT16 padding;
    CGAL_UINT16 origin;
    CGAL_UINT16 halftoning;
    CGAL_UINT32 halftoningParam1;
    CGAL_UINT32 halftoningParam2;
    CGAL_UINT32 colorEncoding;
    CGAL_UINT32 identifier;
} BITMAPHEADER;


/*
 * RGB defines a single color palette entry.  Its analogues in the Windows SDK
 * are the RGBTRIPLE and RGBQUAD structures.  Its analogues in the OS/2
 * Toolkit are the RGB and RGB2 structure. 
 */
typedef struct RGB
{
    CGAL_UINT8 blue;
    CGAL_UINT8 green;
    CGAL_UINT8 red;
} RGB;

#endif   /* __BMPTYPES_H_INCLUDED__ */

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
