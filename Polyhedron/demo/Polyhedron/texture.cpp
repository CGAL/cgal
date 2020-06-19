/***************************************************************************
texture.cpp
----------------------------------------------------------------------------
begin                : june 2003
copyright            : (C) 2003 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
***************************************************************************/

#include <math.h>
#include <string.h>
#include "texture.h"


//////////////////////////////////////////////
// CONSTRUCTORS
//////////////////////////////////////////////

//********************************************
// Constructor
//********************************************
Texture::Texture()
{
  m_pData = NULL;
  m_Width = 0;
  m_WidthByte32 = 0;
  m_Height = 0;
  m_Depth = 0;
}

//********************************************
// Destructor
//********************************************
Texture::~Texture()
{
  Free();
}


//////////////////////////////////////////////
// DATA
//////////////////////////////////////////////

//********************************************
// Alloc
//********************************************
int Texture::Alloc(unsigned int width, unsigned int height, unsigned int depth)
{
  Free();

  unsigned int Width32 = WidthByte32(width,depth);

  m_pData = new unsigned char [Width32 * height];
  if(m_pData == NULL)
  {
    return 0;
  }

  // Set members variables
  m_Width       = width;
  m_WidthByte32 = Width32;
  m_Height      = height;
  m_Depth       = depth;
  UpdateHeader();

  return 1;
}

//********************************************
// Free
//********************************************
void Texture::Free()
{
  if(m_pData != NULL)
  {
    delete [] m_pData;
    m_pData = NULL;
  }
  m_Width = 0;
  m_Height = 0;
  m_Depth = 0;
}



//********************************************
// UpdateWidthByte32
//********************************************
void Texture::UpdateWidthByte32()
{
  m_WidthByte32 = WidthByte32(m_Width,m_Depth);
}

//********************************************
// WidthByte32
//********************************************
unsigned int Texture::WidthByte32(unsigned int width, unsigned int depth)
{
  // 32 bits alignment (4 bytes)
  int rest=(width*depth/8)%4;
  if(rest != 0)
    return (width*depth/8 + 4-rest);
  else
    return (width*depth/8);
}

//********************************************
// UpdateHeader
//********************************************
void Texture::UpdateHeader()
{
  UpdateWidthByte32();
}



//////////////////////////////////////////////
//////////////////////////////////////////////
// CHECKING
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// IsValid
//********************************************
int Texture::IsValid()
{
  int success = 0;
  success = (m_Depth == 24) || (m_Depth == 32);
  success &= (m_Width != 0);
  success &= (m_Height != 0);
  success &= (m_pData != NULL);
  return success;
}

//********************************************
// HigherPowerOfTwo
//********************************************
int Texture::HigherPowerOfTwo(int value)
{
  if(value <= 0)
    return value;

  int power = 1;
  int x = 0;

  while(1)
  {
    x = (int)pow(2.0,(double)power);
    if(x >= value)
      return x;
    power++;
  }
}

//********************************************
// LowerPowerOfTwo
//********************************************
int Texture::LowerPowerOfTwo(int value)
{
  if(value <= 0)
    return value;

  int power = 1;
  int x = 0;

  while(1)
  {
    x = (int)pow(2.0,(double)power);
    if(x >= value)
      return (int)pow(2.0,(double)power-1);
    power++;
  }
}

//********************************************
// SameSize
//********************************************
int Texture::SameSize(Texture *pTexture)
{
  int success = (m_Width == pTexture->GetWidth());
  success &= (m_Height == pTexture->GetHeight());
  return success;
}


//********************************************
// Flip BGR to RGB
//********************************************
int Texture::BGRtoRGB()
{
  if(!IsValid())
    return 0;

  unsigned char pixel;
  int BytePerPixel = m_Depth/8;
  for(unsigned int j=0;j<m_Height;j++)
    for(unsigned int i=0;i<m_Width;i++)
    {
      pixel = m_pData[m_WidthByte32*j+i*BytePerPixel+2];
      m_pData[m_WidthByte32*j+i*BytePerPixel+2] = m_pData[m_WidthByte32*j+i*BytePerPixel];
      m_pData[m_WidthByte32*j+i*BytePerPixel] = pixel;
    }
  return 1;
}

//////////////////////////////////////////////
//////////////////////////////////////////////
// DUPLICATE
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// Extract
//********************************************
int Texture::Extract(int left, int top, int right, int bottom)
{
  // Saturate
  if(right == -1)
    right = m_Width-1;
  if(bottom == -1)
    bottom = m_Height-1;

  // Check
  if(left >= right || top >= bottom)
    return 0;
  if(left < 0  || left >= (int)m_Width ||
    right < 0 || right >= (int)m_Width)
    return 0;
  if(top < 0  || top >= (int)m_Height ||
    bottom < 0 || bottom >= (int)m_Height)
    return 0;

  int NewWidth = right-left+1;
  int NewWidthByte32 = WidthByte32(NewWidth,m_Depth);
  int NewHeight = bottom-top+1;
  int BytePerPixel = m_Depth / 8;
  int i,j,k;

  //TRACE("Start extracting...\n");
  //TRACE("New width : %d\n",NewWidth);
  //TRACE("New height : %d\n",NewHeight);

  // Alloc
  unsigned char *pData = new unsigned char[NewWidthByte32*NewHeight];
  if(pData == NULL)
  {
    //TRACE("Insufficiant memory");
    return 0;
  }

  for(j=0;j<NewHeight;j++)
    for(i=0;i<NewWidth;i++)
      for(k=0;k<BytePerPixel;k++)
        pData[NewWidthByte32*j+i*BytePerPixel+k] = m_pData[m_WidthByte32*(m_Height-1-(j+top))+(i+left)*BytePerPixel+k];

  // Replace datas
  delete [] m_pData;
  m_pData = pData;
  m_Width = NewWidth;
  m_WidthByte32 = NewWidthByte32;
  m_Height = NewHeight;

  UpdateHeader();

  return 1;
}


//********************************************
// DuplicateMirror
//********************************************
int Texture::DuplicateMirror(int left, int top, int right, int bottom)
{

  if(!Extract(left,top,right,bottom))
    return 0;

  left = 0;
  right = m_Width-1;
  top = 0;
  bottom = m_Height-1;

  int NewWidth = 2*m_Width;
  int NewWidthByte32 = WidthByte32(NewWidth,m_Depth);
  int NewHeight = 2*m_Height;
  int BytePerPixel = m_Depth / 8;
  int i,j,k;

  //TRACE("Start duplicate mirror...\n");
  //TRACE("New width : %d\n",NewWidth);
  //TRACE("New widthbyte32 : %d\n",NewWidthByte32);
  //TRACE("New height : %d\n",NewHeight);

  // Alloc
  unsigned char *pData = new unsigned char[NewWidthByte32*NewHeight];
  if(pData == NULL)
  {
    //TRACE("Insufficiant memory");
    return 0;
  }

  // o o
  // x o
  for(j=0;j<NewHeight/2;j++)
    for(i=0;i<NewWidth/2;i++)
      for(k=0;k<BytePerPixel;k++)
        pData[NewWidthByte32*j+i*BytePerPixel+k] = m_pData[m_WidthByte32*(bottom-(j+top))+(i+left)*BytePerPixel+k];
  // o o
  // o x
  for(j=0;j<NewHeight/2;j++)
    for(i=NewWidth/2;i<NewWidth;i++)
      for(k=0;k<BytePerPixel;k++)
        pData[NewWidthByte32*j+i*BytePerPixel+k] = m_pData[m_WidthByte32*(bottom-(j+top))+(right-(i-NewWidth/2+left))*BytePerPixel+k];
  // x o
  // o o
  for(j=NewHeight/2;j<NewHeight;j++)
    for(i=0;i<NewWidth/2;i++)
      for(k=0;k<BytePerPixel;k++)
        pData[NewWidthByte32*j+i*BytePerPixel+k] = m_pData[m_WidthByte32*(j-NewHeight/2+top)+(i+left)*BytePerPixel+k];
  // o x
  // o o
  for(j=NewHeight/2;j<NewHeight;j++)
    for(i=NewWidth/2;i<NewWidth;i++)
      for(k=0;k<BytePerPixel;k++)
        pData[NewWidthByte32*j+i*BytePerPixel+k] = m_pData[m_WidthByte32*(j-NewHeight/2+top)+(right-(i-NewWidth/2+left))*BytePerPixel+k];

  // Replace datas
  delete [] m_pData;
  m_pData = pData;
  m_Width = NewWidth;
  m_WidthByte32 = NewWidthByte32;
  m_Height = NewHeight;

  UpdateHeader();

  return 1;
}




//********************************************
// DuplicateRepeatWidth
//********************************************
int Texture::DuplicateRepeatWidth(int left, int top, int right, int bottom)
{
  if(!Extract(left,top,right,bottom))
    return 0;

  left = 0;
  right = m_Width-1;
  top = 0;
  bottom = m_Height-1;

  int NewWidth = 2*m_Width;
  int NewWidthByte32 = WidthByte32(NewWidth,m_Depth);
  int NewHeight = m_Height;
  int BytePerPixel = m_Depth / 8;
  int i,j,k;

  ////TRACE("Start duplicate repeat width...\n");
  ////TRACE("New width : %d\n",NewWidth);
  ////TRACE("New widthbyte32 : %d\n",NewWidthByte32);
  ////TRACE("New height : %d\n",NewHeight);

  // Alloc
  unsigned char *pData = new unsigned char[NewWidthByte32*NewHeight];
  if(pData == NULL)
  {
    ////TRACE("Insufficiant memory");
    return 0;
  }

  // x o
  for(j=0;j<NewHeight;j++)
    for(i=0;i<NewWidth/2;i++)
      for(k=0;k<BytePerPixel;k++)
        pData[NewWidthByte32*j+i*BytePerPixel+k] = m_pData[m_WidthByte32*(bottom-(j+top))+(i+left)*BytePerPixel+k];
  // o x
  for(j=0;j<NewHeight;j++)
    for(i=NewWidth/2;i<NewWidth;i++)
      for(k=0;k<BytePerPixel;k++)
        pData[NewWidthByte32*j+i*BytePerPixel+k] = m_pData[m_WidthByte32*(bottom-(j+top))+(i-NewWidth/2+left)*BytePerPixel+k];

  // Replace datas
  delete [] m_pData;
  m_pData = pData;
  m_Width = NewWidth;
  m_WidthByte32 = NewWidthByte32;
  m_Height = NewHeight;

  UpdateHeader();

  return 1;
}


//********************************************
// Fill
//********************************************
void Texture::Fill(unsigned char r,
                   unsigned char g,
                   unsigned char b)
{
  if(!IsValid()) return;
  if(m_Depth != 24) return;
  for(unsigned int j=0;j<m_Height;j++)
    for(unsigned int i=0;i<m_Width;i++)
    {
      m_pData[m_WidthByte32*j+i*3] = b;
      m_pData[m_WidthByte32*j+i*3+1] = g;
      m_pData[m_WidthByte32*j+i*3+2] = r;
    }
}


//***************************************
// GreyToColor
//***************************************
void Texture::GreyToColor(unsigned char grey,
                          unsigned char r,
                          unsigned char g,
                          unsigned char b)
{
  if(!IsValid()) return;
  if(m_Depth != 24) return;
  for(unsigned int j=0;j<m_Height;j++)
    for(unsigned int i=0;i<m_Width;i++)
    {
      if(m_pData[m_WidthByte32*j+i*3] == grey)
      {
        m_pData[m_WidthByte32*j+i*3] = b;
        m_pData[m_WidthByte32*j+i*3+1] = g;
        m_pData[m_WidthByte32*j+i*3+2] = r;
      }
    }
}

//***************************************
// ColorToColor
//***************************************
void Texture::ColorToColor(unsigned char r1,
                           unsigned char g1,
                           unsigned char b1,
                           unsigned char r2,
                           unsigned char g2,
                           unsigned char b2)
{
  if(!IsValid()) return;
  if(m_Depth != 24) return;
  for(unsigned int j=0;j<m_Height;j++)
    for(unsigned int i=0;i<m_Width;i++)
    {
      if(m_pData[m_WidthByte32*j+i*3] == b1 &&
        m_pData[m_WidthByte32*j+i*3+1] == g1 &&
        m_pData[m_WidthByte32*j+i*3+2] == r1)
      {
        m_pData[m_WidthByte32*j+i*3]   = b2;
        m_pData[m_WidthByte32*j+i*3+1] = g2;
        m_pData[m_WidthByte32*j+i*3+2] = r2;
      }
    }
}

//////////////////////////////////////////////
//////////////////////////////////////////////
// ALPHA
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// SetAlphaLayer
//********************************************
int Texture::SetAlphaLayer(unsigned char alpha) // 0 - 255
{
  // Check
  if(!IsValid())
    return 0;

  if(m_Depth != 32)
    return 0;

  // Fill alpha layer
  int size = m_Width * m_Height;
  for(int i=0;i<4*size;i+=4)
    m_pData[i+3] = alpha;

  return 1;
}

//********************************************
// AddAlphaLayer
//********************************************
int Texture::AddAlphaLayer(unsigned char alpha) // 0 - 255
{
  // Check
  if(!IsValid())
    return 0;

  // Has soon alpha
  if(HasAlpha())
    return SetAlphaLayer(alpha);

  // Alloc memory
  unsigned char *pData = new unsigned char[4*m_Width*m_Height];
  if(pData == NULL)
  {
    //TRACE("Texture::AddAlphaLayer : insufficiant memory");
    return 0;
  }

  // Fill new data
  int size = m_Width * m_Height;
  int BytePerPixel = m_Depth / 8;
  for(int i=0;i<size;i++)
  {
    pData[4*i+0] = m_pData[BytePerPixel*i+0];
    pData[4*i+1] = m_pData[BytePerPixel*i+1];
    pData[4*i+2] = m_pData[BytePerPixel*i+2];
    pData[4*i+3] = alpha;
  }

  // Set new depth
  m_Depth = 32;

  // Replace datas
  delete [] m_pData;
  m_pData = pData;

  return 1;
}


//********************************************
// SetAlpha
// From RGB to grey scales, then alpha layer
//********************************************
int Texture::SetAlphaLayer(Texture *pTexture)
{
  // Check
  if(!IsValid())
    return 0;
  if(!pTexture->IsValid())
    return 0;

  if(!SameSize(pTexture))
    return 0;

  if(!AddAlphaLayer(0))
    return 0;

  // Fill new data
  unsigned char *pData = pTexture->GetData();
  int size = m_Width * m_Height;
  int BytePerPixel = pTexture->GetDepth() / 8;
  for(int i=0;i<size;i++)
    m_pData[4*i+3] = (unsigned char)((int)pData[BytePerPixel*i+0]+
    (int)pData[BytePerPixel*i+1]+
    (int)pData[BytePerPixel*i+2])/3;

  return 1;
}

//////////////////////////////////////////////
//////////////////////////////////////////////
// DISPLAY
//////////////////////////////////////////////
//////////////////////////////////////////////


//********************************************
// ReadBuffer
//********************************************
int Texture::ReadBuffer(unsigned char *buffer,
                        int width,
                        int height,
                        int depth)
{
  if(buffer == NULL)
    return 0;

  if(!Alloc(width,height,depth))
    return 0;

  int BytePerPixel = depth / 8;

  for(int j=0;j<height;j++)
    for(int i=0;i<width;i++)
      for(int k=0;k<BytePerPixel;k++)
        m_pData[m_WidthByte32*j + i*BytePerPixel+k] =
        buffer[(width*j+i)*BytePerPixel+k];

  return 1;
}


//********************************************
// ReadBufferByte32
//********************************************
int Texture::ReadBufferByte32(unsigned char *pData,
                              int width,
                              int height)
{
  // alloc 32 bits buffer
  if(!Alloc(width,height,32))
    return 0;

  if(pData == NULL)
    return 0;

  memcpy(m_pData,pData,height*m_WidthByte32);
  return 1;
}


//***************************************
// Copy
//***************************************
void Texture::Copy(Texture *pTexture)
{
  unsigned char *pBuffer = pTexture->GetData();
  if(pBuffer == NULL)
    return;

  unsigned int width = pTexture->GetWidth();
  unsigned int height = pTexture->GetHeight();
  unsigned int depth = pTexture->GetDepth();
  if(!Alloc(width,height,depth))
    return;

  unsigned int BytePerPixel = depth / 8;

  for(unsigned int j=0;j<height;j++)
    for(unsigned int i=0;i<width;i++)
      for(unsigned int k=0;k<BytePerPixel;k++)
        m_pData[m_WidthByte32*j + i*BytePerPixel+k] =
        pBuffer[(width*j+i)*BytePerPixel+k];
}


//********************************************
// ReadBuffer
//********************************************
int Texture::ReadBuffer(float *buffer,
                        int width,
                        int height,
                        int depth)
{
  if(buffer == NULL)
    return 0;

  if(!Alloc(width,height,depth))
    return 0;

  int BytePerPixel = depth / 8;

  for(int j=0;j<height;j++)
    for(int i=0;i<width;i++)
      for(int k=0;k<BytePerPixel;k++)
        m_pData[m_WidthByte32*j + i*BytePerPixel+k] =
        (unsigned char)(255.0f * buffer[(width*j+i)*BytePerPixel+k]);

  return 1;
}

//********************************************
// ReadBuffer
//********************************************
int Texture::ReadBuffer(float **ppBuffer,
                        int width,
                        int height,
                        float ratio)
{
  if(ppBuffer == NULL)
    return 0;

  if(!Alloc(width,height,24))
    return 0;

  for(int j=0;j<height;j++)
    for(int i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = (unsigned char)(ratio*ppBuffer[j][i]);
      m_pData[m_WidthByte32*j + i*3+1] = (unsigned char)(ratio*ppBuffer[j][i]);
      m_pData[m_WidthByte32*j + i*3+2] = (unsigned char)(ratio*ppBuffer[j][i]);
    }
  return 1;
}

//********************************************
// WriteBuffer
//********************************************
int Texture::WriteBuffer(float **ppBuffer,
                         int width,
                         int height)
{
  if(ppBuffer == NULL)
    return 0;

  for(int j=0;j<height;j++)
    for(int i=0;i<width;i++) // only first channel
      ppBuffer[j][i] = m_pData[m_WidthByte32*j + i*3];
  return 1;
}

//********************************************
// WriteBuffer32
//********************************************
int Texture::WriteBuffer32(float **ppBuffer,
                           int width,
                           int height)
{
  if(ppBuffer == NULL)
    return 0;
  //ASSERT(m_Depth == 32);
  unsigned int r,g,b;
  int tmp;
  for(int j=0;j<height;j++)
  {
    tmp = m_WidthByte32*j;
    for(int i=0;i<width;i++) // from 3 channels
    {
      b = m_pData[tmp + i*4];
      g = m_pData[tmp + i*4+1];
      r = m_pData[tmp + i*4+2];
      unsigned int value = (r << 16) + (g << 8) + b;
      ppBuffer[j][i] = (float)value;
    }
  }
  return 1;
}

//********************************************
// ReadBuffer
//********************************************
int Texture::ReadBuffer(double *buffer,
                        int width,
                        int height,
                        int depth)
{
  if(buffer == NULL)
    return 0;

  if(!Alloc(width,height,depth))
    return 0;

  int BytePerPixel = depth / 8;

  for(int j=0;j<height;j++)
    for(int i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*BytePerPixel]   = (unsigned char)buffer[width*j+i];
      m_pData[m_WidthByte32*j + i*BytePerPixel+1] = (unsigned char)buffer[width*j+i];
      m_pData[m_WidthByte32*j + i*BytePerPixel+2] = (unsigned char)buffer[width*j+i];
    }

  return 1;
}


//********************************************
// Grey
//********************************************
int Texture::Grey(unsigned int x, unsigned int y)
{
  int BytePerPixel = m_Depth / 8;
  // Grey scale
  if(BytePerPixel == 1)
    return (int)m_pData[m_WidthByte32*y + x];
  else // 24 or 32 bits (alpha layer)
    return (int)((int)m_pData[m_WidthByte32*y + x*BytePerPixel+0]+
    (int)m_pData[m_WidthByte32*y + x*BytePerPixel+1]+
    (int)m_pData[m_WidthByte32*y + x*BytePerPixel+2])/3;
}

//********************************************
// Color
//********************************************
void Texture::Color(unsigned int x, unsigned int y,
                    unsigned char *pRed, unsigned char *pGreen, unsigned char *pBlue)
{
  int BytePerPixel = m_Depth / 8;
  // Grey scale
  if(BytePerPixel == 1)
  {
    *pRed = m_pData[m_WidthByte32*y + x];
    *pGreen = m_pData[m_WidthByte32*y + x];
    *pBlue = m_pData[m_WidthByte32*y + x];
  }
  else // 24 or 32 bits (alpha layer)
  {
    *pRed = m_pData[m_WidthByte32*y + x*BytePerPixel];
    *pGreen = m_pData[m_WidthByte32*y + x*BytePerPixel+1];
    *pBlue = m_pData[m_WidthByte32*y + x*BytePerPixel+2];
  }
}





void Texture::GenerateMirrorV(unsigned int width,
                              unsigned int height,
                              unsigned char r,
                              unsigned char g,
                              unsigned char b,
                              unsigned char rb,
                              unsigned char gb,
                              unsigned char bb)
{
  Alloc(width,height,24);

  for(unsigned j=0;j<height;j++)
    for(unsigned int i=0;i<width/4;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = r;
      m_pData[m_WidthByte32*j + i*3+1] = g;
      m_pData[m_WidthByte32*j + i*3+2] = b;
    }
  for(unsigned j=0;j<height;j++)
    for(unsigned int i=width/4;i<3*width/4;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = rb;
      m_pData[m_WidthByte32*j + i*3+1] = gb;
      m_pData[m_WidthByte32*j + i*3+2] = bb;
    }
  for(unsigned j=0;j<height;j++)
    for(unsigned int i=3*width/4;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = r;
      m_pData[m_WidthByte32*j + i*3+1] = g;
      m_pData[m_WidthByte32*j + i*3+2] = b;
    }
}

void Texture::GenerateMirrorH(unsigned int width,
                              unsigned int height,
                              unsigned char r,
                              unsigned char g,
                              unsigned char b,
                              unsigned char rb,
                              unsigned char gb,
                              unsigned char bb)
{
  Alloc(width,height,24);

  for(unsigned j=0;j<height/4;j++)
    for(unsigned int i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = r;
      m_pData[m_WidthByte32*j + i*3+1] = g;
      m_pData[m_WidthByte32*j + i*3+2] = b;
    }
  for(unsigned j=height/4;j<3*height/4;j++)
    for(unsigned int i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = rb;
      m_pData[m_WidthByte32*j + i*3+1] = gb;
      m_pData[m_WidthByte32*j + i*3+2] = bb;
    }
  for(unsigned j=3*height/4;j<height;j++)
    for(unsigned int i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = r;
      m_pData[m_WidthByte32*j + i*3+1] = g;
      m_pData[m_WidthByte32*j + i*3+2] = b;
  }
}


void Texture::GenerateCheckerBoard(unsigned int width,
                                   unsigned int height,
                                   int size,
                                   unsigned char r,
                                   unsigned char g,
                                   unsigned char b,
                                   unsigned char rb,
                                   unsigned char gb,
                                   unsigned char bb)
{
  Alloc(width,height,24);

  for(unsigned j=0;j<height;j++)
    for(unsigned int i=0;i<width;i++)
    {
      int h = (int)((double)i/(double)size);
      int v = (int)((double)j/(double)size);
      bool h_odd = (h%2) == 0;
      bool v_odd = (v%2) == 0;
      if((h_odd        + v_odd)%2 == 0)
      {
        m_pData[m_WidthByte32*j + i*3]   = rb;
        m_pData[m_WidthByte32*j + i*3+1] = gb;
        m_pData[m_WidthByte32*j + i*3+2] = bb;
      }
      else
      {
        m_pData[m_WidthByte32*j + i*3]   = r;
        m_pData[m_WidthByte32*j + i*3+1] = g;
        m_pData[m_WidthByte32*j + i*3+2] = b;
      }
    }
}

void Texture::GenerateVStripes(unsigned int width,
                               unsigned int height,
                               int size,
                               unsigned char r,
                               unsigned char g,
                               unsigned char b,
                               unsigned char rb,
                               unsigned char gb,
                               unsigned char bb)
{
  Alloc(width,height,24);

  for(unsigned j=0;j<height;j++)
    for(unsigned int i=0;i<width;i++)
    {
      int h = (int)((double)i/(double)size);
      bool h_odd = (h%2) == 0;
      if(h_odd)
      {
        m_pData[m_WidthByte32*j + i*3]   = rb;
        m_pData[m_WidthByte32*j + i*3+1] = gb;
        m_pData[m_WidthByte32*j + i*3+2] = bb;
      }
      else
      {
        m_pData[m_WidthByte32*j + i*3]   = r;
        m_pData[m_WidthByte32*j + i*3+1] = g;
        m_pData[m_WidthByte32*j + i*3+2] = b;
      }
    }
}

void Texture::GenerateHStripes(unsigned int width,
                               unsigned int height,
                               int size,
                               unsigned char r,
                               unsigned char g,
                               unsigned char b,
                               unsigned char rb,
                               unsigned char gb,
                               unsigned char bb)
{
  Alloc(width,height,24);

  for(unsigned j=0;j<height;j++)
    for(unsigned int i=0;i<width;i++)
    {
      int v = (int)((double)j/(double)size);
      bool v_odd = (v%2) == 0;
      if(v_odd)
      {
        m_pData[m_WidthByte32*j + i*3]   = rb;
        m_pData[m_WidthByte32*j + i*3+1] = gb;
        m_pData[m_WidthByte32*j + i*3+2] = bb;
      }
      else
      {
        m_pData[m_WidthByte32*j + i*3]   = r;
        m_pData[m_WidthByte32*j + i*3+1] = g;
        m_pData[m_WidthByte32*j + i*3+2] = b;
      }
    }
}


//***************************************
// GenerateGrid
//***************************************
void Texture::GenerateGrid(unsigned int width,
                           unsigned int height,
                           int /*size*/,
                           unsigned int thickness,
                           unsigned char r,
                           unsigned char g,
                           unsigned char b,
                           unsigned char rb,
                           unsigned char gb,
                           unsigned char bb)
{
  Alloc(width,height,24);

  // fill background
  unsigned int i,j;
  for(j=0;j<height;j++)
    for(i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = rb;
      m_pData[m_WidthByte32*j + i*3+1] = gb;
      m_pData[m_WidthByte32*j + i*3+2] = bb;
    }

  // horizontal
  for(j=height/2-thickness/2;j<height/2+thickness/2;j++)
  {
    for(i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = r;
      m_pData[m_WidthByte32*j + i*3+1] = g;
      m_pData[m_WidthByte32*j + i*3+2] = b;
    }
  }

  // vertical
  for(unsigned int i=width/2-thickness/2;i<width/2+thickness/2;i++)
    for(unsigned int j=0;j<height;j++)
    {
      m_pData[m_WidthByte32*j + i*3]   = r;
      m_pData[m_WidthByte32*j + i*3+1] = g;
      m_pData[m_WidthByte32*j + i*3+2] = b;
    }
}

void Texture::GenerateGradientH(unsigned int width,
                                unsigned int height,
                                int /*size*/)
{
  Alloc(width,height,24);

  // fill background
  for(unsigned int i=0;i<width;i++)
  {
    unsigned char g = (unsigned char)(255.0 * (double)i / (double)(width-1));
    for(unsigned int j=0;j<height;j++)
    {
      m_pData[m_WidthByte32*j + i*3]   = g;
      m_pData[m_WidthByte32*j + i*3+1] = g;
      m_pData[m_WidthByte32*j + i*3+2] = g;
    }
  }
}
void Texture::GenerateGradientV(unsigned int width,
                                unsigned int height,
                                int /*size*/)
{
  Alloc(width,height,24);

  // fill background
  for(unsigned int j=0;j<height;j++)
  {
    unsigned char g = (unsigned char)(255.0 * (double)j / (double)(height-1));
    for(unsigned int i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = g;
      m_pData[m_WidthByte32*j + i*3+1] = g;
      m_pData[m_WidthByte32*j + i*3+2] = g;
    }
  }
}

//***************************************
// GenerateLines
//***************************************
void Texture::GenerateLines(unsigned int width,
                            unsigned int height,
                            int size,
                            unsigned char r,
                            unsigned char g,
                            unsigned char b,
                            unsigned char rb,
                            unsigned char gb,
                            unsigned char bb)
{
  Alloc(width,height,24);

  // fill background
  unsigned int i,j;
  for(j=0;j<height;j++)
    for(i=0;i<width;i++)
    {
      m_pData[m_WidthByte32*j + i*3]   = rb;
      m_pData[m_WidthByte32*j + i*3+1] = gb;
      m_pData[m_WidthByte32*j + i*3+2] = bb;
    }

  // horizontal zebra
  for(j=0;j<height;j++)
  {
    if((j/size)%2 == 0)
      for(i=0;i<width;i++)
      {
        m_pData[m_WidthByte32*j + i*3]   = r;
        m_pData[m_WidthByte32*j + i*3+1] = g;
        m_pData[m_WidthByte32*j + i*3+2] = b;
      }
  }
}






// ** EOF **



