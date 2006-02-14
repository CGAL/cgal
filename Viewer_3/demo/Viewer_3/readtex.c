/* readtex.c */

/*
 * Read an SGI .rgb image file and generate a mipmap texture set.
 * Much of this code was borrowed from SGI's tk OpenGL toolkit.
 */



#include <GL/gl.h>
#include <GL/glu.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h>


#ifndef SEEK_SET
#  define SEEK_SET 0
#endif


/*
** RGB Image Structure
*/

typedef struct _TK_RGBImageRec {
   GLint sizeX, sizeY;
   GLint components;
   unsigned char *data;
} TK_RGBImageRec;



/******************************************************************************/

typedef struct _rawImageRec {
    unsigned short imagic;
    unsigned short type;
    unsigned short dim;
    unsigned short sizeX, sizeY, sizeZ;
    unsigned long min, max;
    unsigned long wasteBytes;
    char name[80];
    unsigned long colorMap;
    FILE *file;
    unsigned char *tmp, *tmpR, *tmpG, *tmpB, *tmpA;
    unsigned long rleEnd;
    GLuint *rowStart;
    GLint *rowSize;
} rawImageRec;

/******************************************************************************/

static void ConvertShort(unsigned short *array, long length)
{
   unsigned long b1, b2;
   unsigned char *ptr;

   ptr = (unsigned char *)array;
   while (length--) {
      b1 = *ptr++;
      b2 = *ptr++;
      *array++ = (unsigned short) ((b1 << 8) | (b2));
   }
}

static void ConvertLong(GLuint *array, long length)
{
   unsigned long b1, b2, b3, b4;
   unsigned char *ptr;

   ptr = (unsigned char *)array;
   while (length--) {
      b1 = *ptr++;
      b2 = *ptr++;
      b3 = *ptr++;
      b4 = *ptr++;
      *array++ = (b1 << 24) | (b2 << 16) | (b3 << 8) | (b4);
   }
}

static rawImageRec *RawImageOpen(const char *fileName)
{
   union {
      int testWord;
      char testByte[4];
   } endianTest;
   rawImageRec *raw;
   GLenum swapFlag;
   int x;

   endianTest.testWord = 1;
   if (endianTest.testByte[0] == 1) {
      swapFlag = GL_TRUE;
   } else {
      swapFlag = GL_FALSE;
   }

   raw = (rawImageRec *)malloc(sizeof(rawImageRec));
   if (raw == NULL) {
      fprintf(stderr, "Out of memory!\n");
      return NULL;
   }
   if ((raw->file = fopen(fileName, "rb")) == NULL) {
      perror(fileName);
      return NULL;
   }

   fread(raw, 1, 12, raw->file);

   if (swapFlag) {
      ConvertShort(&raw->imagic, 6);
   }

   raw->tmp = (unsigned char *)malloc(raw->sizeX*256);
   raw->tmpR = (unsigned char *)malloc(raw->sizeX*256);
   raw->tmpG = (unsigned char *)malloc(raw->sizeX*256);
   raw->tmpB = (unsigned char *)malloc(raw->sizeX*256);
   if (raw->sizeZ==4) {
      raw->tmpA = (unsigned char *)malloc(raw->sizeX*256);
   }
   if (raw->tmp == NULL || raw->tmpR == NULL || raw->tmpG == NULL ||
       raw->tmpB == NULL) {
      fprintf(stderr, "Out of memory!\n");
      return NULL;
   }

   if ((raw->type & 0xFF00) == 0x0100) {
      x = raw->sizeY * raw->sizeZ * sizeof(GLuint);
      raw->rowStart = (GLuint *)malloc(x);
      raw->rowSize = (GLint *)malloc(x);
      if (raw->rowStart == NULL || raw->rowSize == NULL) {
         fprintf(stderr, "Out of memory!\n");
         return NULL;
      }
      raw->rleEnd = 512 + (2 * x);
      fseek(raw->file, 512, SEEK_SET);
      fread(raw->rowStart, 1, x, raw->file);
      fread(raw->rowSize, 1, x, raw->file);
      if (swapFlag) {
         ConvertLong(raw->rowStart, (long) (x/sizeof(GLuint)));
         ConvertLong((GLuint *)raw->rowSize, (long) (x/sizeof(GLint)));
      }
   }
   return raw;
}

static void RawImageClose(rawImageRec *raw)
{

   fclose(raw->file);
   free(raw->tmp);
   free(raw->tmpR);
   free(raw->tmpG);
   free(raw->tmpB);
   if (raw->sizeZ>3) {
      free(raw->tmpA);
   }
   free(raw);
}

static void RawImageGetRow(rawImageRec *raw, unsigned char *buf, int y, int z)
{
   unsigned char *iPtr, *oPtr, pixel;
   int count, done = 0;

   if ((raw->type & 0xFF00) == 0x0100) {
      fseek(raw->file, (long) raw->rowStart[y+z*raw->sizeY], SEEK_SET);
      fread(raw->tmp, 1, (unsigned int)raw->rowSize[y+z*raw->sizeY],
            raw->file);
      
      iPtr = raw->tmp;
      oPtr = buf;
      while (!done) {
         pixel = *iPtr++;
         count = (int)(pixel & 0x7F);
         if (!count) {
			 done = 1;
            return;
         }
         if (pixel & 0x80) {
            while (count--) {
               *oPtr++ = *iPtr++;
            }
         } else {
            pixel = *iPtr++;
            while (count--) {
               *oPtr++ = pixel;
            }
         }
      }
   } else {
      fseek(raw->file, 512+(y*raw->sizeX)+(z*raw->sizeX*raw->sizeY),
            SEEK_SET);
      fread(buf, 1, raw->sizeX, raw->file);
   }
}


static void RawImageGetData(rawImageRec *raw, TK_RGBImageRec *final)
{
   unsigned char *ptr;
   int i, j;

   final->data = (unsigned char *)malloc((raw->sizeX+1)*(raw->sizeY+1)*4);
   if (final->data == NULL) {
      fprintf(stderr, "Out of memory!\n");
   }

   ptr = final->data;
   for (i = 0; i < (int)(raw->sizeY); i++) {
      RawImageGetRow(raw, raw->tmpR, i, 0);
      RawImageGetRow(raw, raw->tmpG, i, 1);
      RawImageGetRow(raw, raw->tmpB, i, 2);
      if (raw->sizeZ>3) {
         RawImageGetRow(raw, raw->tmpA, i, 3);
      }
      for (j = 0; j < (int)(raw->sizeX); j++) {
         *ptr++ = *(raw->tmpR + j);
         *ptr++ = *(raw->tmpG + j);
         *ptr++ = *(raw->tmpB + j);
         if (raw->sizeZ>3) {
            *ptr++ = *(raw->tmpA + j);
         }
      }
   }
}


static TK_RGBImageRec *tkRGBImageLoad(const char *fileName)
{
   rawImageRec *raw;
   TK_RGBImageRec *final;

   raw = RawImageOpen(fileName);
   if (!raw) {
      fprintf(stderr, "File not found\n");
      return NULL;
   }
   final = (TK_RGBImageRec *)malloc(sizeof(TK_RGBImageRec));
   if (final == NULL) {
      fprintf(stderr, "Out of memory!\n");
      return NULL;
   }
   final->sizeX = raw->sizeX;
   final->sizeY = raw->sizeY;
   final->components = raw->sizeZ;
   RawImageGetData(raw, final);
   RawImageClose(raw);
   return final;
}


static void FreeImage( TK_RGBImageRec *image )
{
   free(image->data);
   free(image);
}


/*
 * Load an SGI .rgb file and generate a set of 2-D mipmaps from it.
 * Input:  imageFile - name of .rgb to read
 *         intFormat - internal texture format to use, or number of components
 * Return:  GL_TRUE if success, GL_FALSE if error.
 */
GLboolean LoadRGBMipmaps( const char *imageFile, GLint intFormat )
{
   GLint error;
   GLenum format;
   TK_RGBImageRec *image;

   image = tkRGBImageLoad( imageFile );
   if (!image) {
      return GL_FALSE;
   }

   if (image->components==3) {
      format = GL_RGB;
   }
   else if (image->components==4) {
      format = GL_RGBA;
   }
   else {
      /* not implemented */
      fprintf(stderr,
              "Error in LoadRGBMipmaps %d-component images not implemented\n",
              image->components );
      return GL_FALSE;
   }

   error = gluBuild2DMipmaps( GL_TEXTURE_2D,
                              intFormat,
                              image->sizeX, image->sizeY,
                              format,
                              GL_UNSIGNED_BYTE,
                              image->data );

   FreeImage(image);
   return error ? GL_FALSE : GL_TRUE;
}



/*
 * Load an SGI .rgb file and return a pointer to the image data.
 * Input:  imageFile - name of .rgb to read
 * Output:  width - width of image
 *          height - height of image
 *          format - format of image (GL_RGB or GL_RGBA)
 * Return:  pointer to image data or NULL if error
 */
GLubyte *LoadRGBImage( const char *imageFile, GLint *width, GLint *height,
                       GLenum *format )
{
   TK_RGBImageRec *image;
   GLint bytes;
   GLubyte *buffer;

   image = tkRGBImageLoad( imageFile );
   if (!image) {
      return NULL;
   }

   if (image->components==3) {
      *format = GL_RGB;
   }
   else if (image->components==4) {
      *format = GL_RGBA;
   }
   else {
      /* not implemented */
      fprintf(stderr,
              "Error in LoadRGBImage %d-component images not implemented\n",
              image->components );
      return NULL;
   }

   *width = image->sizeX;
   *height = image->sizeY;

   bytes = image->sizeX * image->sizeY * image->components;
   buffer = (GLubyte *) malloc(bytes);
   if (!buffer)
      return NULL;

   memcpy( (void *) buffer, (void *) image->data, bytes );

   FreeImage(image);

   return buffer;
}

