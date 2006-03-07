/*
 * from bmp.zip, see the url http://www.ddj.com/ftp/1995/1995.03/
 * author Dr. Dobb's
 */

/*
 * This is the header for endian.c - functions to read/write our
 * INT8, INT16 and INT32 types from/to a little-endian file.
 */

#ifndef __ENDIAN_H_INCLUDED__
#define __ENDIAN_H_INCLUDED__
 
/*
 * Read the basic types as little-endian values.  The return value will be
 * zero if successful, EOF, otherwise.
 */
extern int readINT8little(FILE *f, INT8 *i);
extern int readINT16little(FILE *f, INT16 *i);
extern int readINT32little(FILE *f, INT32 *i);
extern int readUINT8little(FILE *f, UINT8 *i);
extern int readUINT16little(FILE *f, UINT16 *i);
extern int readUINT32little(FILE *f, UINT32 *i);

/*
 * Write the basic types as little-endian values.  The return value will be
 * zero if successful, EOF otherwise.
 */
int writeINT8little(FILE *f, INT8 i);
int writeINT16little(FILE *f, INT16 i);
int writeINT32little(FILE *f, INT32 i);
int writeUINT8little(FILE *f, UINT8 i);
int writeUINT16little(FILE *f, UINT16 i);
int writeUINT32little(FILE *f, UINT32 i);

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

