
/*
 * These functions read and write our basic integer types from a little-endian
 * file.  The endian and word-size of the host machine will not affect this
 * code.  The only assumption made is that the C data type (char) is one byte
 * long.  This should be a safe assumption.
 */

#include <stdio.h>
#include "bmptypes.h"
#include "bmpendian.h"

/*****************************************************************************
*
* Read functions.  All read functions take an open file pointer as the first
* parameter and a pointer to data as the second parameter.  The return value
* will be 0 on success, and EOF on failure.  If successful, the second
* parameter will point to the data read.
*/

/*
 * The INT8 and UINT8 types are stored as a single byte on disk.  The INT8
 * type is a signed integer with range (-128..127).  The UINT8 type is an
 * unsigned integer with range (0..255).
 */
int readINT8little(FILE *f, INT8 *i)
{
    int rc;
    
    rc = fgetc(f);
    if (rc == EOF)
	return rc;
    
    *i = (rc & 0xff);
    return 0;
}

int readUINT8little(FILE *f, UINT8 *i)
{
    int  rc;
    
    rc = fgetc(f);
    if (rc == EOF)
	return rc;
    
    *i = (rc & 0xff);
    return 0;
}


/*
 * The INT16 and UINT16 types are stored as two bytes on disk.  The INT16 type
 * is a signed integer with range (-32768..32767).  The UINT16 type is an
 * unisgned integer with range (0..65535).
 */
int readINT16little(FILE *f, INT16 *i)
{
    int rc;
    INT16 temp = 0;
    
    temp = (fgetc(f) & 0xff);
    
    rc = fgetc(f);
    if (rc == EOF)
	return rc;
    
    temp |= ((rc & 0xff) << 8);
    *i = temp;
    return 0;
}

int readUINT16little(FILE *f, UINT16 *i)
{
    int rc;
    UINT16 temp = 0;
    
    temp = (fgetc(f) & 0xff);
    
    rc = fgetc(f);
    if (rc == EOF)
	return rc;
    
    temp |= ((rc & 0xff) << 8);
    *i = temp;
    return 0;
}

/*
 * The INT32 and UINT32 types are stored as four bytes on disk.  The INT32
 * type is a signed integer with range (-2147483648..2147483647).  The UINT32
 * type is an unisgned integer with range (0..4294967295).
 */
int readINT32little(FILE *f, INT32 *i)
{
    int rc;
    INT32 temp = 0;
    
    temp = ((long)fgetc(f) & 0xff);
    temp |= (((long)fgetc(f) & 0xff) << 8);
    temp |= (((long)fgetc(f) & 0xff) << 16);
    
    rc = fgetc(f);
    if (rc == EOF)
	return rc;
    
    temp |= (((long)rc & 0xff) << 24);
    *i = temp;
    return 0;
}

int readUINT32little(FILE *f, UINT32 *i)
{
    int rc;
    UINT32 temp = 0;
    
    temp = ((long)fgetc(f) & 0xff);
    temp |= (((long)fgetc(f) & 0xff) << 8);
    temp |= (((long)fgetc(f) & 0xff) << 16);
    
    rc = fgetc(f);
    if (rc == EOF)
	return rc;
    
    temp |= (((long)rc & 0xff) << 24);
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

int writeINT8little(FILE *f, INT8 i)
{
    return fputc(i, f);
}

int writeUINT8little(FILE *f, UINT8 i)
{
    return fputc(i, f);
}

int writeINT16little(FILE *f, INT16 i)
{
    int rc;
    
    rc = fputc((i & 0xff), f);
    if (rc == EOF)
	return rc;
    
    return fputc(((i >> 8) & 0xff), f);
}

int writeUINT16little(FILE *f, UINT16 i)
{
    int rc;
    
    rc = fputc((i & 0xff), f);
    if (rc == EOF)
	return rc;
    
    return fputc(((i >> 8) & 0xff), f);
}

int writeINT32little(FILE *f, INT32 i)
{
    int rc;
    
    rc = fputc((i & 0xff), f);
    if (rc == EOF)
	return rc;
    
    rc = fputc(((i >> 8) & 0xff), f);
    if (rc == EOF)
	return rc;
    
    rc = fputc(((i >> 16) & 0xff), f);
    if (rc == EOF)
	return rc;
    
    return fputc(((i >> 24) & 0xff), f);
}


int writeUINT32little(FILE *f, UINT32 i)
{
    int rc;
    
    rc = fputc((i & 0xff), f);
    if (rc == EOF)
	return rc;
    
    rc = fputc(((i >> 8) & 0xff), f);
    if (rc == EOF)
	return rc;
    
    rc = fputc(((i >> 16) & 0xff), f);
    if (rc == EOF)
	return rc;
    
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

