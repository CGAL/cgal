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
 * This file contains mid-level functions for reading bitmap structures and
 * high-level functions that read bitmap files.
 */

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

/*****************************************************************************
 *
 * Mid-level functions.
 *
 * These functions read in the various structures defined in bmptypes.h.  Each
 * function assumes that the file pointer is positioned at the start of the
 * given structure.  Upon completion, each function will leave the file
 * pointer positioned on the byte immediately following the structure.  Return
 * values will be 0 for success and non-zero for failure.  In all cases, a
 * return value of non-zero means that the file position has been left in an
 * indeterminate state and further reading should not be attempted.
 */

/*
 * Read a Bitmapfileheader structure.
 */
CGAL_INLINE_FUNCTION
int readBitmapFileHeader(FILE *fp, Bitmapfileheader *bfh)
{
    int rc;
  
    rc = readUINT16little(fp, &(bfh->type));
    if (rc != 0)
	return rc;
    
    rc = readUINT32little(fp, &(bfh->size));
    if (rc != 0)
	return rc;
    
    rc = readINT16little(fp, &(bfh->xHotspot));
    if (rc != 0)
	return rc;
    
    rc = readINT16little(fp, &(bfh->yHotspot));
    if (rc != 0)
	return rc;
    
    rc = readUINT32little(fp, &(bfh->offsetToBits));
    return rc;
}

/*
 * Read a BITMAPARRAYHEADER
 */
CGAL_INLINE_FUNCTION
int readBitmapArrayHeader(FILE *fp, BITMAPARRAYHEADER *bah)
{
    int rc;
    
    rc = readUINT16little(fp, &(bah->type));
    if (rc != 0)
	return rc;
    rc = readUINT32little(fp, &(bah->size));
    if (rc != 0)
	return rc;
    rc = readUINT32little(fp, &(bah->next));
    if (rc != 0)
	return rc;
    rc = readUINT16little(fp, &(bah->screenWidth));
    if (rc != 0)
	return rc;
    rc = readUINT16little(fp, &(bah->screenHeight));
    return rc;
}

/*
 * Read the BITMAPHEADER structure.  This one requires a bit of sanity
 * checking.  The length of the structure on the disk is specified in the
 * first field.  We must stop reading after that many bytes, and if that value
 * is longer than sizeof(BITMAPHEADER), we must skip over any leftover bytes.
 * Finally, if size is 12, then width an height are really 16-bit values, and
 * we have to read them differently so they'll be properly stored in the
 * 32-bit fields BITMAPHEADER uses.
 */
CGAL_INLINE_FUNCTION
int readBitmapHeader(FILE *fp, BITMAPHEADER *bh)
{
    int    rc, oldFormat;
	unsigned int bytesRead;
    CGAL_UINT16 tempVal;
    
    /*
     * Clear the structure.  Default values for all fields are zeros.  This
     * will prevent garbage from being returned if the structure is truncated
     * on disk.
     */
    memset(bh, 0, sizeof(BITMAPHEADER));
    
    /*
     * Read the size of the structure.  From here on in, we'll have to be sure
     * we don't read more bytes than this value.
     */
    rc = readUINT32little(fp, &(bh->size));
    if (rc != 0)
	return rc;
    bytesRead = 4;
    
    /*
     * If the size is 12 bytes or less, than this is an "old format"
     * structure.  So the width and height fields will have to be read
     * differently.
     */
    if (bh->size <= 12)
	oldFormat = 1;
    else
	oldFormat = 0;
    
    /*
     * Width and height are read differently for old and new format files.  In
     * the old format, they're 16-bit values.  In the new format, they're
     * 32-bits long.
     */
    if (oldFormat)
    {
	rc = readUINT16little(fp, &tempVal);
	if (rc != 0)
	    return rc;
	bh->width = tempVal;
	bytesRead += 2;
    }
    else
    {
	rc = readINT32little(fp, &(bh->width));
	if (rc != 0)
	    return rc;
	bytesRead += 4;
    }
    if (bytesRead >= bh->size)
	return 0;
    
    if (oldFormat)
    {
	rc = readUINT16little(fp, &tempVal);
	if (rc != 0)
	    return rc;
	bh->height = tempVal;
	bytesRead += 2;
    }
    else
    {
	rc = readINT32little(fp, &(bh->height));
	if (rc != 0)
	    return rc;
	bytesRead += 4;
    }
    if (bytesRead >= bh->size)
	return 0;
    
    /*
     * From this point on, old and new formats are identical to each other,
     * and we can proceed as if there was no difference.  For each field, we
     * read it in and increment the count of bytes read.  If at any time we
     * have read the amount we got earlier (in the size field), then stop and
     * leave the rest of the fields as zeros.
     */
    rc = readUINT16little(fp, &(bh->numBitPlanes));
    if (rc != 0)
	return rc;
    bytesRead += 2;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT16little(fp, &(bh->numBitsPerPlane));
    if (rc != 0)
	return rc;
    bytesRead += 2;
    if (bytesRead >= bh->size)
	return 0;
  
    /*
     * Old format stop here.  But we don't have to check, because in that
     * format, 12 bytes have been read and the function will have exited 
     * without any extra checking.
     */
    rc = readUINT32little(fp, &(bh->compressionScheme));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT32little(fp, &(bh->sizeOfImageData));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT32little(fp, &(bh->xResolution));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT32little(fp, &(bh->yResolution));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT32little(fp, &(bh->numColorsUsed));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT32little(fp, &(bh->numImportantColors));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT16little(fp, &(bh->resolutionUnits));
    if (rc != 0)
	return rc;
    bytesRead += 2;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT16little(fp, &(bh->padding));
    if (rc != 0)
	return rc;
    bytesRead += 2;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT16little(fp, &(bh->origin));
    if (rc != 0)
	return rc;
    bytesRead += 2;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT16little(fp, &(bh->halftoning));
    if (rc != 0)
	return rc;
    bytesRead += 2;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT32little(fp, &(bh->halftoningParam1));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT32little(fp, &(bh->halftoningParam2));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
  
    rc = readUINT32little(fp, &(bh->colorEncoding));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    rc = readUINT32little(fp, &(bh->identifier));
    if (rc != 0)
	return rc;
    bytesRead += 4;
    if (bytesRead >= bh->size)
	return 0;
    
    /*
     * If there are more bytes in the file than this, then the file is using a
     * future format that doesn't exist yet.  Skip over the bytes.  Assuming
     * this future format somewhat resembles what we know now, ignoring the
     * fields will be safe.  We _MUST_ skip them, though, since the color
     * table begins on the byte after this structure, and we have to position
     * the file pointer there.
     */
    return fseek(fp, (bh->size - bytesRead), SEEK_CUR);
}


/*
 * readRgb reads in a single RGB structure from the disk.  The numBytes field
 * indicates how many bytes the field occupies on the disk.  It assumes that
 * each component is one byte on disk and the rest is padding.  This will
 * compensate for the old/new differences in color tables.  (Old format
 * bitmaps use 3 bytes per entry, while new format bitmaps use 4.)  Note how
 * it will never read more than the number of bytes requested.
 */
CGAL_INLINE_FUNCTION
int readRgb(FILE *fp, RGB *rgb, int numBytes)
{
    int rc;
    
    if (numBytes == 0)
	return 0;
    rc = readUINT8little(fp, &(rgb->blue));
    if (rc != 0)
	return rc;
    
    if (numBytes == 1)
	return 0;
    rc = readUINT8little(fp, &(rgb->green));
    if (rc != 0)
	return rc;
    
    if (numBytes == 2)
	return 0;
    rc = readUINT8little(fp, &(rgb->red));
    if (rc != 0)
	return rc;
    
    if (numBytes == 3)
	return 0;
    
    /* Skip over extra bytes if more than three were requested */
    return fseek(fp, (numBytes - 3), SEEK_CUR);
}

/*
 * A color table is a block of RGB structures, all the same size.  Read it by
 * calling readRgb repeatedly.
 */
CGAL_INLINE_FUNCTION
int readColorTable(FILE *fp, RGB *rgb, int numEntries, int numBytesPerEntry)
{
    int i, rc;
    
    for (i=0; i<numEntries; i++)
    {
	rc = readRgb(fp, &(rgb[i]), numBytesPerEntry);
	if (rc != 0)
	    return rc;
    }
    return 0;
}


/*
 * ReadBitsUncompressed.  Reads pixel values into an array of RGB
 * values.  It assmes that there is no compression.  Note that there we're
 * only handling bit depths of 1, 4, 8, 16, and 24. Note that OS/2 bitmaps
 * can (in theory) have any number of bits per pixel, so you might find a
 * strange bitmap file that this can't handle, but the chances of finding such
 * a file this are nearly zero.
 *
 * Each row of pixels is always padded to a 4-byte boundary.
 */
CGAL_INLINE_FUNCTION
int readBitsUncompressed(FILE *fp, RGB *image, int width, int height,
			 int depth, RGB *colorTable)
{
    CGAL_UINT8 temp;
    int   rc, padBytes, i;
    long  row, column, pixel, value;
    
    switch (depth) {
    case 1:
	/*
	 * For 1 bit per pixel, each byte is 8 pixels.  Each one is an index
	 * into the color table (1 or 0).  Most significant byte first.  All
	 * is padded to 32-bit boundaries as well.
	 */
	pixel = 0;
	if (((width % 32) == 0) || ((width % 32) > 24))
	    padBytes = 0;
	else if ((width % 32) <= 8)
	    padBytes = 3;
	else if ((width % 32) <= 16)
	    padBytes = 2;
	else
	    padBytes = 1;

	for (row = height; row > 0; row--)
	{
	    for (column = width; column > 0; column -= 8)
	    {
		rc = readUINT8little(fp, &temp);
		if (rc != 0)
		    return rc;
		for (i=0; i < ((column < 8) ? column : 8); i++)
		{
		    /*
		     * For each byte read, bit-decompose it.  Note that the
		     * last byte on a row could have less than 8 bits used.
		     * Most significant bits come first.
		     */
		    value = ((temp & (1 << (7-i))) == 0) ? 0 : 1;
		    image[pixel].red   = colorTable[value].red;
		    image[pixel].green = colorTable[value].green;
		    image[pixel].blue  = colorTable[value].blue;
		    pixel++;
		}
	    }
	    if (padBytes != 0)
	    {
		rc = fseek(fp, padBytes, SEEK_CUR);
		if (rc != 0)
		    return rc;
	    }
	}
	break;

    case 4:
	/*
	 * For 4 bits per pixel, each byte is two pixels.  The upper half go to
	 * the first pixel, and the lower half to the second.
	 */
	pixel = 0;
	if (((width % 8) == 0) || ((width % 8) > 6))
	    padBytes = 0;
	else if ((width % 8) <= 2)
	    padBytes = 3;
	else if ((width % 8) <= 4)
	    padBytes = 2;
	else
	    padBytes = 1;
	
	for (row = height; row > 0; row--)
	{
	    for (column = width; column > 0; column -= 2)
	    {
		/*
		 * Each byte here is two pixels.  Note that the last byte on a
		 * row may only contain one pixel.
		 */
		rc = readUINT8little(fp, &temp);
		if (rc != 0)
		    return rc;
		/*
		 * First pixel is the upper 4 bits
		 */
		value = temp >> 4;
		image[pixel].red   = colorTable[value].red;
		image[pixel].green = colorTable[value].green;
		image[pixel].blue  = colorTable[value].blue;
		pixel++;

		/*
		 * Second pixel is lower 4 bits.  If this is the last byte in
		 * the row, and there are an odd number of pixels per row, then
		 * this is not valid data.
		 */
		if (column == 1)
		{
		    value = temp & 0x0f;
		    image[pixel].red   = colorTable[value].red;
		    image[pixel].green = colorTable[value].green;
		    image[pixel].blue  = colorTable[value].blue;
		    pixel++;
		}
	    }
	    if (padBytes != 0)
	    {
		rc = fseek(fp, padBytes, SEEK_CUR);
		if (rc != 0)
		    return rc;
	    }
	}
	break;

    case 8:
	/*
	 * For 8 bits per pixel, each byte is one pixel.
	 */
	pixel = 0;
	padBytes = ((width % 4) == 0) ? 0 : (4 - (width % 4));

	for (row=height; row > 0; row--)
	{
	    for (column=width; column > 0; column--)
	    {
		rc = readUINT8little(fp, &temp);
		if (rc != 0)
		    return rc;
		image[pixel].red   = colorTable[temp].red;
		image[pixel].green = colorTable[temp].green;
		image[pixel].blue  = colorTable[temp].blue;
		pixel++;
	    }
	    if (padBytes != 0)
	    {
		rc = fseek(fp, padBytes, SEEK_CUR);
		if (rc != 0)
		    return rc;
	    }
	}
	break;

    case 16:
	/*
	 * For 16 bits per pixel, you must read two bytes per pixel.  But
	 * there's a catch. The data is big endian!  This is because all pixel
	 * data (for all formats, actually) is stored as a packed array,
	 * stored in pixel order.
	 */
	pixel = 0;
	padBytes = ((width % 2) == 0) ? 0 : 2;
	for (row=height; row > 0; row--)
	{
	    for (column=width; column > 0; column--)
	    {
		/*
		 * Read a 16-bit integer as big endian.  Do this by reading
		 * two bytes and mathematically combine them.  After that,
		 * proceed as usual.
		 */
		rc = readUINT8little(fp, &temp);
		if (rc != 0)
		    return rc;
		value = ((long)temp) << 8;
		rc = readUINT8little(fp, &temp);
		if (rc != 0)
		    return rc;
		value |= temp;

		image[pixel].red   = colorTable[value].red;
		image[pixel].green = colorTable[value].green;
		image[pixel].blue  = colorTable[value].blue;
		pixel++;
	    }
	    if (padBytes != 0)
	    {
		rc = fseek(fp, padBytes, SEEK_CUR);
		if (rc != 0)
		    return rc;
	    }
	}
	break;

    case 24:
	/*
	 * For 24 bits per pixel, it's an RGB structure.  Note that the color
	 * table is ignore for bit depths greater than 24 bits.
	 */
	pixel = 0;
	padBytes = width % 4;

	for (row=height; row > 0; row--)
	{
	    for (column=width; column > 0; column--)
	    {
		rc = readRgb(fp, image+pixel, 3);
		pixel++;
	    }
	    if (padBytes != 0)
	    {
		rc = fseek(fp, padBytes, SEEK_CUR);
		if (rc != 0)
		    return rc;
	    }
	}          
	break;
    }
    
    return 0;
}


/*
 * ReadMaskBitsUncompressed.  Reads a monochrome mask into an array of
 * characters.  It assmes that there is no compression.  This is very similar
 * (internally) to the readBitsUncompressed function.  Note that if the data
 * read isn't really one-bit-deep data, you'll probably get garbage back.
 */
CGAL_INLINE_FUNCTION
int readMaskBitsUncompressed(FILE *fp, char *image, int width, int height)
{
    CGAL_UINT8 temp;
    int   rc, padBytes, i;
    long  row, column, pixel;
    char value;
    
    /*
     * see the one-bit-depth part of readBitsUncompressed for comments
     */
    pixel = 0;
    if (((width % 32) == 0) || ((width % 32) > 24))
	padBytes = 0;
    else if ((width % 32) <= 8)
	padBytes = 3;
    else if ((width % 32) <= 16)
	padBytes = 2;
    else
	padBytes = 1;

    for (row = height; row > 0; row--)
    {
	for (column = width; column > 0; column -= 8)
	{
	    rc = readUINT8little(fp, &temp);
	    if (rc != 0)
		return rc;
	    for (i=0; i < ((column < 8) ? column : 8); i++)
	    {
		value = ((temp & (1 << (7-i))) == 0) ? 0 : 1;
		image[pixel] = value;
		pixel++;
	    }
	}
	if (padBytes != 0)
	{
	    rc = fseek(fp, padBytes, SEEK_CUR);
	    if (rc != 0)
		return rc;
	}
    }

    return 0;
}


/*
 * reflectYRGB takes an array of RGB vales and the dimensions they represent
 * and flips it vertically.  This will convert a bottom-left origin array to a
 * top-left origin array.
 */
CGAL_INLINE_FUNCTION
void reflectYRGB(RGB *image, int width, int height)
{
    int row, col;
    RGB temp;
    
    for (row = 0; row < (height / 2); row++)
    {
	for (col = 0; col < width; col++)
	{
	    /* Swap pixels at (x,y) with (x,height-y) */
	    memcpy(&temp, image+(row * width + col), sizeof(RGB));
	    memcpy(image+(row * width + col),
		   image+((height - row - 1) * width + col), sizeof(RGB));
	    memcpy(image+((height - row - 1) * width + col), &temp,
		   sizeof(RGB));
	}
    }
}
/*
 * reflectYchar takes an array of char values and the dimensions they
 * represent and flips it vertically.  This will convert a bottom-left origin
 * array to a top-left origin array.
 */
CGAL_INLINE_FUNCTION
void reflectYchar(char *image, int width, int height)
{
    int row, col;
    char temp;
    
    for (row = 0; row < (height / 2); row++)
    {
	for (col = 0; col < width; col++)
	{
	    /* Swap values at (x,y) with (x,height-y) */
	    temp = image[row * width + col];
	    image[row * width + col]=image[(height - row - 1) * width + col];
	    image[(height - row - 1) * width + col] = temp;
	}
    }
}

/*****************************************************************************
 *
 * High-level functions
 *
 * These functions read in specific types of bitmap files.  Each assumes that
 * the file pointer is positioned at the appropriate place in a bitmap file.
 * (At the start of a Bitmapfileheader for all functions except
 * readMultipleImages, which assumes the file pointer to be positioned on the
 * start of a BITMAPARRAYHEADER.  These functions will leave the file pointer
 * on the byte after the image's color table.
 *
 * The coordinate speaces in the returned arrays will have an upper-left
 * origin.   As before, a non-zero return value indicates that something went
 * wrong.
 *
 * Note that the BMP and mono-ICO functions will not return 1000 if the image
 * is of type color-icon.  This is because a color icon consists of a bitmap
 * and a monochrome icon.
 *
 * return values:
 *	   0 - success
 *	1000 - incorrect file type for the routine called
 *	1001 - image data out of range or damaged file
 *	1002 - good data, but the routine called can't handle it (yet)
 *	1003 - out of memory allocating color table
 *	1004 - out of memory allocating image
 *	1005 - out of memory allocating image arrays
 *	1006 - Illegal image type in a multi-image array
 *
 *     other - I/O error of some kind
 */


/*
 * readSingleImageBMP will read a single BMP image and generate an array of RGB
 * triples that contain the RGB values for every pixel.  It will also return
 * the dimensions of the image.
 */
CGAL_INLINE_FUNCTION
int readSingleImageBMP(FILE *fp, RGB **argb, CGAL_UINT32 *width, CGAL_UINT32 *height)
{
    Bitmapfileheader  bfh;
    BITMAPHEADER      bh;
    RGB              *colorTable = (RGB*)NULL;
    RGB              *image = (RGB*)NULL;
    int               rc, depth, inverted;
    int               numColors;
    long              numPixels, endPos;
    
    /*
     * First, get the file header and sanity check it.  The type field must be
     * TYPE_BMP or we're aborting.
     */
    rc = readBitmapFileHeader(fp, &bfh);
    if (rc != 0)
	return rc;
    if ((bfh.type != TYPE_BMP) &&
	(bfh.type != TYPE_ICO_COLOR) &&
	(bfh.type != TYPE_PTR_COLOR))
	return 1000;
    
    /*
     * Immediately following a file header is always the bitmap header.  Read
     * it.  Sanity check any values we might need.  Specifically, less than
     * 32-bit depth, known compression scheme, known origin, and known color
     * encoding, and valid height/width.  Note that negative height is OK,
     * that indicates an upper-left origin for a Windows bitmap.
     */
    rc = readBitmapHeader(fp, &bh);
    if (rc != 0)
	return rc;
    depth = bh.numBitPlanes * bh.numBitsPerPlane;
    if ((depth > 32) ||
	(bh.compressionScheme > COMPRESSION_LAST) ||
	(bh.origin > ORIGIN_LAST) ||
	(bh.colorEncoding > COLOR_ENCODING_LAST) ||
	(bh.width < 1) ||
	(bh.height == 0))
	return 1001;
    
    /*
     * If the height is negative, then this is a Windows bitmap whose origin
     * is the upper-left corner and not the lower-left.  The inverted flag
     * indicates a lower-left origin.  Our code will be outputting an
     * upper-left origin pixel array.
     */
    if (bh.height < 0)
    {
	inverted = 0;
	bh.height = -bh.height;
    }
    else
	inverted = 1;
    
    /*
     * Now, sanity check a few fields that are valid, but I don't have code to
     * deal with them yet.  This includes: more than one bit plane, any
     * compression scheme, and bit depths that are not 1, 4, 8, 16, or 24.
     */
    if ((bh.numBitPlanes > 1) ||
	((bh.numBitsPerPlane != 1) &&
	 (bh.numBitsPerPlane != 4) &&
	 (bh.numBitsPerPlane != 8) &&
	 (bh.numBitsPerPlane != 16) &&
	 (bh.numBitsPerPlane != 24)) ||
	(bh.compressionScheme != COMPRESSION_NONE))
	return 1002;
    
    /*
     * Allocate and read the color table.  The file pointer has been
     * positioned in the right place by the readBitmapHeader function.  Note
     * that images with 24-bits or more color depth have no color table.  They
     * are  already RGB.  When reading the color table, be sure to check for
     * old/new format bitmaps. 
     */
    if (depth < 24)
    {
	numColors = 1 << depth;
	colorTable = (RGB *)calloc(numColors, sizeof(RGB));
	if (colorTable == NULL)
	    return 1003;
	if (bh.size <= 12)
	    rc = readColorTable(fp, colorTable, numColors, 3);
	else
	    rc = readColorTable(fp, colorTable, numColors, 4);
	if (rc != 0)
	{
	    free(colorTable);
	    return rc;
	}
    }
    
    /*
     * We're at the end of the color table.  Preserve this position.  We'll
     * need to leave the file pointer there before returning from this
     * function. 
     */
    endPos = ftell(fp);

    /*
     * Allocate the array of pixels and fill it in.
     */
    numPixels = bh.width * bh.height;
    image = (RGB *)calloc(numPixels, sizeof(RGB));
    if (image == NULL)
    {
	free (colorTable);
	return 1004;
    }
    
    /*
     * Seek to the bits
     */
    rc = fseek(fp, bfh.offsetToBits, SEEK_SET);
    if (rc != 0)
    {
	free (colorTable);
	free (image);
	return rc;
    }
    
    /*
     * Read the bits.  If code for decompressing bits should be written,
     * insert the call here.
     */
    switch ((int)bh.compressionScheme) {
    case COMPRESSION_NONE:
	rc = readBitsUncompressed(fp, image, bh.width, bh.height, depth,
				  colorTable);
	break;
    }
    
    if (rc != 0)
    {
	free(image);
	return rc;
    }
    
    /*
     * If the origin is lower-left, flip the image upside down
     */
    if (inverted)
	reflectYRGB(image, bh.width, bh.height);

    /*
     * Return the output values.  Set the file pointer to the byte after the
     * color table.
     */
    *argb = image;
    *width = bh.width;
    *height = bh.height;
    fseek(fp, endPos, SEEK_SET);

    /*
     * Clean up and return.  Note that we don't return the color table.  This
     * is because we're returning an array of RGB values for the image - such
     * a table would be redundant.
     */
    if (colorTable != NULL)
	free(colorTable);

    return 0;
}


/*
 * Read in a monochrome OS/2 icon/pointer.  return two arrays of bytes
 * (interpreted as booleans) for the XOR and AND masks.
 */
CGAL_INLINE_FUNCTION
int readSingleImageICOPTR(FILE *fp, char **xorMask, char **andMask,
		          CGAL_UINT32 *width, CGAL_UINT32 *height) 
{
    Bitmapfileheader  bfh;
    BITMAPHEADER      bh;
    char             *mask1, *mask2;
    int               rc;
    long              numPixels, endPos;

    /*
     * Read and sanity check the header.  Monochrom OS/2 icons are TYPE_ICO.
     * Monochrome pointers are TYPE_PTR.  Color ICO and PTR is also allowed,
     * because monochrome masks are part of those images.
     */
    rc = readBitmapFileHeader(fp, &bfh);
    if (rc != 0)
	return rc;
    if ((bfh.type != TYPE_ICO) &&
	(bfh.type != TYPE_PTR) &&
	(bfh.type != TYPE_ICO_COLOR) &&
	(bfh.type != TYPE_PTR_COLOR))
	return 1000;

    /*
     * Now read the bitmap data and sanity check it.  Since this is a
     * monochrome icon, bit depth must be 1.  Additionally, a known
     * compression scheme, known origin, known color encoding, and valid
     * height/width.  Height can't be less than 0, as it can with color
     * images, since this code is only used on for OS/2-type images.
     */
    rc = readBitmapHeader(fp, &bh);
    if (rc != 0)
	return rc;
    if ((bh.numBitPlanes != 1) ||
	(bh.numBitsPerPlane != 1) ||
	(bh.compressionScheme > COMPRESSION_LAST) ||
	(bh.origin > ORIGIN_LAST) ||
	(bh.colorEncoding > COLOR_ENCODING_LAST) ||
	(bh.width < 1) ||
	(bh.height < 1))
	return 1001;

    /*
     * Sanity check some valid fields that I can't deal with yet.
     */
    if (bh.compressionScheme != COMPRESSION_NONE)
	return 1002;

    /*
     * Skip over the color table, since this is a monochrome mask.  Note that
     * the size is already known - 2 entries - which is 6 or 8 bytes.
     * this isn't, and we don't.
     */
    if (bh.size <= 12)
	rc = fseek(fp, 6, SEEK_CUR);
    else
	rc = fseek(fp, 8, SEEK_CUR);
    if (rc != 0)
    {
	return rc;
    }

    /*
     * Save this file position.  we'll have to seek back to it after reading
     * in the image data.
     */
    endPos = ftell(fp);

    /*
     * The image is actually two images.  The top half is an AND mask and the
     * bottom half is an XOR mask.  Allocate the images.
     */
    numPixels = bh.width * bh.height / 2;
    mask1 = (char *)malloc(numPixels);
    if (mask1 == NULL)
	return 1004;
    mask2 = (char *)malloc(numPixels);
    if (mask2 == NULL)
    {
	free(mask1);
	return 1004;
    }

    /*
     * Seek to the bit data
     */
    rc = fseek(fp, bfh.offsetToBits, SEEK_SET);
    if (rc != 0)
    {
	free (mask1);
	free (mask2);
	return rc;
    }

    /*
     * Read in the bits.  Note: since these are really two images, two calls
     * to readMaskBitsUncompressed are made, and the height used is 1/2 the
     * height mentioned in the header.
     */
    switch ((int) bh.compressionScheme) {
    case COMPRESSION_NONE:
	rc = readMaskBitsUncompressed(fp, mask1, bh.width, bh.height/2);
	if (rc == 0)
	    rc = readMaskBitsUncompressed(fp, mask2, bh.width, bh.height/2);
	break;
    }

    if (rc != 0)
    {
	free (mask1);
	free (mask2);
	return rc;
    }

    /*
     * A mask will never have an upper-left origin, since Windows will never
     * produce one in a bitmap file.
     */
    reflectYchar(mask1, bh.width, bh.height / 2);
    reflectYchar(mask2, bh.width, bh.height / 2);

    /*
     * Return everything we've read.
     */
    *xorMask = mask1;
    *andMask = mask2;
    *width = bh.width;
    *height = bh.height / 2;
    fseek(fp, endPos, SEEK_SET);
    
    return 0;
}


/*
 * Read in a color OS/2 icon.  return an array of RGBs for the colors.
 * and two arrays of bytes (interpreted as booleans) for the XOR and AND
 * masks.
 */
CGAL_INLINE_FUNCTION
int readSingleImageColorICOPTR(FILE *fp, RGB **argb, char **xorMask,
			       char **andMask, CGAL_UINT32 *width, CGAL_UINT32 *height)
{
    CGAL_UINT32 width1, height1, width2, height2;
    int rc;

    /*
     * Color icons consist of a monochrome icon followed by a bitmap.  This
     * makes reading them easy - first do a monochrome mask read, and then do
     * a color bitmap read.  We should probably add some more descriptive
     * error codes here.
     *
     * First read the mask.
     */
    rc = readSingleImageICOPTR(fp, xorMask, andMask, &width2, &height2);
    if (rc != 0)
    {
	return rc;
    }

    /*
     * Next, read the color bitmap part
     */
    rc = readSingleImageBMP(fp, argb, &width1, &height1);
    if (rc != 0)
    {
	return rc;
    }

    /*
     * Now, just sanity check the image - the dimensions for the image should
     * match the dimensions of the masks.
     */
    if ((width1 != width2) || (height1 != height2))
	return 1001;

    *width = width1;
    *height = height1;
    return 0;
}


/*
 * readMultipleImage runs down the list of images in a file and returns them
 * all.  ImageCount is the number of images in the file.  The other returned
 * values are parallel arrays.  If an element in aargb, axorMask, or aandMask
 * is NULL, then that image has no such array.  (Bitmaps have no xor or and
 * masks, monochrome icons have no color arrays.
 *
 * Note that on errors other than 1000 and 1005, the arrays will contain good
 * data - the images that have been read properly will be in the arrays.
 * Images that have not yet been read will consist of NULL pointers in the
 * arrays. 
 */
CGAL_INLINE_FUNCTION
int readMultipleImage(FILE *fp, RGB ***argbs, char ***xorMasks,
		      char ***andMasks, CGAL_UINT32 **widths, CGAL_UINT32 **heights,
		      int *imageCount)
{
    int rc;
    long filePos;
    BITMAPARRAYHEADER bah;
    CGAL_UINT16 imageType;
    int count;
    
    /*
     * First count the images.  Preserve the file position for later.  If some
     * structure in the list isn't an array header, return 1000.
     */
    filePos = ftell(fp);
    count = 0;
    do
    {
	rc = readBitmapArrayHeader(fp, &bah);
	if (rc != 0)
	    return rc;
	if (bah.type != TYPE_ARRAY)
	    return 1000;
	fseek(fp, bah.next, SEEK_SET);
	count++;
    } while (bah.next != 0);
    fseek(fp, filePos, SEEK_SET);
    
    /*
     * Allocate the arrays.  Return 1005 on any failures
     */
    *argbs = (RGB **)calloc(count, sizeof(RGB *));
    if (*argbs == NULL)
	return 1005;
    *xorMasks = (char **)calloc(count, sizeof(char *));
    if (*xorMasks == NULL)
    {
	free(*argbs);
	return 1005;
    }
    *andMasks = (char **)calloc(count, sizeof(char *));
    if (*andMasks == NULL)
    {
	free(*argbs);
	free(*xorMasks);
	return 1005;
    }
    *widths = (CGAL_UINT32 *)calloc(count, sizeof(CGAL_UINT32));
    if (*widths == NULL)
    {
	free(*argbs);
	free(*xorMasks);
	free(*andMasks);
	return 1005;
    }
    *heights = (CGAL_UINT32 *)calloc(count, sizeof(CGAL_UINT32));
    if (*heights == NULL)
    {
	free(*argbs);
	free(*xorMasks);
	free(*andMasks);
	free(*widths);
	return 1005;
    }
    *imageCount = count;

    /*
     * Loop through the images again.  This time, read each image
     */
    count = 0;
    do
    {
	rc = readBitmapArrayHeader(fp, &bah);
	if (rc != 0)
	    return rc;
	/*
	 * Get the image type.  Preserve the position, since we're reading
	 * into the next structure.
	 */
	filePos = ftell(fp);
	rc = readUINT16little(fp, &imageType);
	if (rc != 0)
	    return rc;
	rc = fseek(fp, filePos, SEEK_SET);

	/*
	 * Now that we know what kind of image we're about to read, read it. 
	 */
	switch(imageType) {
	case TYPE_BMP:
	    rc = readSingleImageBMP(fp, (*argbs)+count, (*widths)+count,
				    (*heights)+count);
	    break;
	case TYPE_ICO:
	case TYPE_PTR:
	    rc = readSingleImageICOPTR(fp, (*xorMasks)+count,
				       (*andMasks)+count, (*widths)+count,
				       (*heights)+count);
	    break;
	case TYPE_ICO_COLOR:
	case TYPE_PTR_COLOR:
	    rc = readSingleImageColorICOPTR(fp, (*argbs)+count,
					    (*xorMasks)+count,
					    (*andMasks)+count,
					    (*widths)+count,
					    (*heights)+count);
	    break;
	default:
	    return 1006;
	}
	if (rc != 0)
	    return rc;
	
	fseek(fp, bah.next, SEEK_SET);
	count++;
    } while (bah.next != 0);

    return 0;
}


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

