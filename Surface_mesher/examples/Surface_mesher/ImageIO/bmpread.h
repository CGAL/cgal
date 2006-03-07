/*
 * from bmp.zip, see the url http://www.ddj.com/ftp/1995/1995.03/
 * author Dr. Dobb's
 */

/*
 * This is the header for readbmp.c - functions to read the bitmap file
 * structures.  See readbmp.c for details.
 */

#ifndef __READBMP_H_INCLUDED__
#define __READBMP_H_INCLUDED__

/*
 * Mid-level functions
 */
int readBitmapFileHeader(FILE *fp, BITMAPFILEHEADER *bfh);
int readBitmapArrayHeader(FILE *fp, BITMAPARRAYHEADER *bah);
int readBitmapHeader(FILE *fp, BITMAPHEADER *bh);
int readRgb(FILE *fp, RGB *rgb, int numBytes);
int readColorTable(FILE *fp, RGB *rgb, int numEntries, int numBytesPerEntry);

int readBitsUncompressed(FILE *fp, RGB *image, int width, int height,
			 int depth, RGB* colorTable);
int readMaskBitsUncompressed(FILE *fp, char *image, int width, int height);

void reflectYRGB(RGB *image, int width, int height);
void reflectYchar(char *image, int width, int height);

/*
 * High level functions.
 */
int readSingleImageBMP(FILE *fp, RGB **argb, UINT32 *width, UINT32 *height);
int readSingleImageICOPTR(FILE *fp, char **xorMask, char **andMask,
		          UINT32 *width, UINT32 *height);
int readSingleImageColorICOPTR(FILE *fp, RGB **argb, char **xorMask,
			       char **andMask, UINT32 *width, UINT32 *height);
int readMultipleImage(FILE *fp, RGB ***argbs, char ***xorMasks,
		      char ***andMasks, UINT32 **widths, UINT32 **heights,
		      int *imageCount);

#endif

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

