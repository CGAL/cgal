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

#ifdef CGAL_HEADER_ONLY
#define CGAL_INLINE_FUNCTION inline
#else
#define CGAL_INLINE_FUNCTION
#endif

/*
 * These functions read and write our basic integer types from a little-endian
 * file.  The endian and word-size of the host machine will not affect this
 * code.  The only assumption made is that the C data type (char) is one byte
 * long.  This should be a safe assumption.
 */

#include <stdio.h>
#include <CGAL/ImageIO/bmptypes.h>

/*****************************************************************************
*
* Read functions.  All read functions take an open file pointer as the first
* parameter and a pointer to data as the second parameter.  The return value
* will be 0 on success, and EOF on failure.  If successful, the second
* parameter will point to the data read.
*/

/*
 * The CGAL_INT8 and CGAL_UINT8 types are stored as a single byte on disk.  The INT8
 * type is a signed integer with range (-128..127).  The CGAL_UINT8 type is an
 * unsigned integer with range (0..255).
 */
CGAL_INLINE_FUNCTION
int readINT8little(FILE *f, CGAL_INT8 *i)
{
    int rc;

    rc = fgetc(f);
    if (rc == EOF)
        return EOF;

    *i = CGAL_INT8(rc & 0xff);
    return 0;
}

CGAL_INLINE_FUNCTION
int readUINT8little(FILE *f, CGAL_UINT8 *i)
{
    int  rc;

    rc = fgetc(f);
    if (rc == EOF)
        return EOF;

    *i = CGAL_UINT8(rc & 0xff);
    return 0;
}


/*
 * The CGAL_INT16 and CGAL_UINT16 types are stored as two bytes on disk.  The INT16 type
 * is a signed integer with range (-32768..32767).  The CGAL_UINT16 type is an
 * unisgned integer with range (0..65535).
 */
CGAL_INLINE_FUNCTION
int readINT16little(FILE *f, CGAL_INT16 *i)
{
    int rc;
    CGAL_INT16 temp = 0;

    temp = (fgetc(f) & 0xff);

    rc = fgetc(f);
    if (rc == EOF)
        return EOF;

    temp = temp | CGAL_INT16((rc & 0xff) << 8);
    *i = temp;
    return 0;
}

CGAL_INLINE_FUNCTION
int readUINT16little(FILE *f, CGAL_UINT16 *i)
{
    int rc;
    CGAL_UINT16 temp = 0;

    temp = (fgetc(f) & 0xff);

    rc = fgetc(f);
    if (rc == EOF)
        return EOF;

    temp = CGAL_INT16(temp | ((rc & 0xff) << 8));
    *i = temp;
    return 0;
}

/*
 * The CGAL_INT32 and CGAL_UINT32 types are stored as four bytes on disk.  The INT32
 * type is a signed integer with range (-2147483648..2147483647).  The CGAL_UINT32
 * type is an unisgned integer with range (0..4294967295).
 */
CGAL_INLINE_FUNCTION
int readINT32little(FILE *f, CGAL_INT32 *i)
{
    int rc;
    CGAL_INT32 temp = 0;

    temp = ((long)fgetc(f) & 0xff);
    temp = CGAL_INT32(temp | (((long)fgetc(f) & 0xff) << 8));
    temp = CGAL_INT32(temp | (((long)fgetc(f) & 0xff) << 16));

    rc = fgetc(f);
    if (rc == EOF)
        return EOF;

    temp = CGAL_INT32(temp | (((long)rc & 0xff) << 24));
    *i = temp;
    return 0;
}

CGAL_INLINE_FUNCTION
int readUINT32little(FILE *f, CGAL_UINT32 *i)
{
    int rc;
    CGAL_UINT32 temp = 0;

    temp = ((long)fgetc(f) & 0xff);
    temp = CGAL_UINT32(temp | (((long)fgetc(f) & 0xff) << 8));
    temp = CGAL_UINT32(temp | (((long)fgetc(f) & 0xff) << 16));

    rc = fgetc(f);
    if (rc == EOF)
        return EOF;

    temp = CGAL_UINT32(temp | (((long)rc & 0xff) << 24));
    *i = temp;
    return 0;
}

/*****************************************************************************
*
* Write functions.  All write functions take an open file pointer as the first
* parameter and a data as the second parameter.  The return value will be 0 on
* success, and EOF on failure.  If successful, the second parameter will have
* been written to the open file.
*/

CGAL_INLINE_FUNCTION
int writeINT8little(FILE *f, CGAL_INT8 i)
{
    return fputc(i, f);
}

CGAL_INLINE_FUNCTION
int writeUINT8little(FILE *f, CGAL_UINT8 i)
{
    return fputc(i, f);
}

CGAL_INLINE_FUNCTION
int writeINT16little(FILE *f, CGAL_INT16 i)
{
    int rc;

    rc = fputc((i & 0xff), f);
    if (rc == EOF)
        return EOF;

    return fputc(((i >> 8) & 0xff), f);
}

CGAL_INLINE_FUNCTION
int writeUINT16little(FILE *f, CGAL_UINT16 i)
{
    int rc;

    rc = fputc((i & 0xff), f);
    if (rc == EOF)
        return EOF;

    return fputc(((i >> 8) & 0xff), f);
}

CGAL_INLINE_FUNCTION
int writeINT32little(FILE *f, CGAL_INT32 i)
{
    int rc;

    rc = fputc((i & 0xff), f);
    if (rc == EOF)
        return EOF;

    rc = fputc(((i >> 8) & 0xff), f);
    if (rc == EOF)
        return EOF;

    rc = fputc(((i >> 16) & 0xff), f);
    if (rc == EOF)
        return EOF;

    return fputc(((i >> 24) & 0xff), f);
}


CGAL_INLINE_FUNCTION
int writeUINT32little(FILE *f, CGAL_UINT32 i)
{
    int rc;

    rc = fputc((i & 0xff), f);
    if (rc == EOF)
        return EOF;

    rc = fputc(((i >> 8) & 0xff), f);
    if (rc == EOF)
        return EOF;

    rc = fputc(((i >> 16) & 0xff), f);
    if (rc == EOF)
        return EOF;

    return fputc(((i >> 24) & 0xff), f);
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

