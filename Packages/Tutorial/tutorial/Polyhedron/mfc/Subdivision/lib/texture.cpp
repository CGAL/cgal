/***************************************************************************
texture.cpp
----------------------------------------------------------------------------
begin                : june 2003
copyright            : (C) 2003 by Pierre Alliez - INRIA
email                : pierre.alliez@sophia.inria.fr
***************************************************************************/

#include "stdafx.h"
#include <math.h>
#include <string.h>
#include "texture.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#endif


//////////////////////////////////////////////
// CONSTRUCTORS
//////////////////////////////////////////////

//********************************************
// Constructor
//********************************************
CTexture::CTexture()
{
  m_pData = NULL;
  m_Width = 0;
  m_WidthByte32 = 0;
  m_Height = 0;
  m_Depth = 0;
  m_pFileName = new char[MAX_PATH];
  strcpy(m_pFileName,"");
}

//********************************************
// Destructor
//********************************************
CTexture::~CTexture()
{
  Free();
  delete [] m_pFileName;
}


//////////////////////////////////////////////
// DATA
//////////////////////////////////////////////

//********************************************
// Alloc
//********************************************
int CTexture::Alloc(unsigned int width, unsigned int height, unsigned int depth)
{
  Free();

  unsigned int Width32 = WidthByte32(width,depth);

  m_pData = new unsigned char [Width32 * height];
  if(m_pData == NULL)
  {
    TRACE("CTexture::Alloc : Insufficiant memory\n");
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
void CTexture::Free()
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


//////////////////////////////////////////////
//////////////////////////////////////////////
// FILE READING
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// ReadFile (dispatch function)
//********************************************
int CTexture::ReadFile(char *filename, unsigned int width, unsigned int height, unsigned int depth)
{
  // Cleanup
  Free();

  // Storage
  strcpy(m_pFileName,filename);

  // Extension
  TRACE("CTexture::ReadFile : file : %s\n",filename);
  int len = strlen(filename);
  char extension[10];
  strcpy(extension,&(filename[len-4]));

  if(extension == ".bmp")
    return ReadFileBMP(filename);
  if(extension == ".raw")
    return ReadFileRAW(filename,width,height,depth);

  return 0;
}


//********************************************
// ReadFileBMP (*.bmp)
//********************************************
// Read windows bmp files
// Accept only 24 bits
// Size : 2^n x 2^m
//********************************************
int CTexture::ReadFileBMP(char *filename)
{
  
  // Check for valid bmp file
  CFile file;
  CFileException ex;

  // Try to open file
  if(!file.Open(filename, CFile::modeRead | CFile::typeBinary,&ex))
  {
    TRACE("Unable to open file for reading");
    return 0;
  }

  // File header
  BITMAPFILEHEADER FileHeader;
  file.Read(&FileHeader,sizeof(BITMAPFILEHEADER));
  TRACE("FileHeader.bfType : %d\n",FileHeader.bfType);
  TRACE("FileHeader.bfSize : %d\n",FileHeader.bfSize);
  TRACE("FileHeader.bfReserved1 : %d\n",FileHeader.bfReserved1);
  TRACE("FileHeader.bfReserved2 : %d\n",FileHeader.bfReserved2);
  TRACE("FileHeader.bfOffBits : %d\n",FileHeader.bfOffBits);

  // Is it a Windows BMP file ? (BM)
  WORD sign = ((WORD) ('M' << 8) | 'B');
  if(FileHeader.bfType != sign)
  {
    TRACE("Invalid BMP file");
    file.Close();
    return 0;
  }

  file.Read(&m_Header,sizeof(BITMAPINFOHEADER));
  TRACE("\n");
  TRACE("IMAGE HEADER :\n");
  TRACE("biSize : %d\n",m_Header.biSize);
  TRACE("biWidth : %d\n",m_Header.biWidth);
  TRACE("biHeight : %d\n",m_Header.biHeight);
  TRACE("biPlanes : %d\n",m_Header.biPlanes);
  TRACE("biBitCount : %d\n",m_Header.biBitCount);
  TRACE("biCompression : %d\n",m_Header.biCompression);
  TRACE("biSizeImage : %d\n",m_Header.biSizeImage);
  TRACE("biXPelsPerMeter : %d\n",m_Header.biXPelsPerMeter);
  TRACE("biYPelsPerMeter : %d\n",m_Header.biYPelsPerMeter);
  TRACE("biClrUsed : %d\n",m_Header.biClrUsed);
  TRACE("biClrImportant : %d\n",m_Header.biClrImportant);

  // 24 bits ?
  if(m_Header.biPlanes != 1 ||
     m_Header.biBitCount != 24)
  {
    TRACE("Texture file must have 24 bits depth");
    file.Close();
    return 0;
  }

  // Alloc (does call Free before)
  Free();
  m_pData = new unsigned char[m_Header.biSizeImage];
  if(m_pData == NULL)
  {
    TRACE("Insufficiant memory");
    file.Close();
    return 0;
  }

  // Update datas
  m_Width = m_Header.biWidth;
  m_Height = m_Header.biHeight;
  m_Depth = m_Header.biBitCount;

  // Image reading
  file.Read(m_pData,m_Header.biSizeImage);

  // Close file
  file.Close();

  UpdateWidthByte32();
 
  return 1;
}

//********************************************
// UpdateWidthByte32
//********************************************
void CTexture::UpdateWidthByte32()
{
  m_WidthByte32 = WidthByte32(m_Width,m_Depth);
}

//********************************************
// WidthByte32
//********************************************
unsigned int CTexture::WidthByte32(unsigned int width, unsigned int depth)
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
void CTexture::UpdateHeader()
{
  UpdateWidthByte32();

  m_Header.biWidth = m_Width;
  m_Header.biHeight = m_Height;
  m_Header.biSizeImage = m_WidthByte32 * m_Height;

  m_Header.biSize = 40;
  m_Header.biPlanes = 1;
  m_Header.biBitCount = m_Depth;
  m_Header.biCompression = (WORD)0;
  m_Header.biXPelsPerMeter = 0;
  m_Header.biYPelsPerMeter = 0;
  m_Header.biClrUsed = 0;
  m_Header.biClrImportant = 0;
}


//********************************************
// ReadFileRAW (*.raw)
//********************************************
// Read raw files
// Accept only 24 or 32 bits
// Size : 2^n x 2^m
//********************************************
int CTexture::ReadFileRAW(char *filename, unsigned int width, unsigned int height, unsigned int depth)
{
  // Check for valid file
  FILE *fp = fopen(filename,"rb");

  // Try to open file
  if(!fp)
  {
    TRACE("Unable to open file for reading");
    return 0;
  }

  // Alloc (does call Free before)
  if(!Alloc(width,height,depth))
  {
    TRACE("Insufficiant memory");
    fclose(fp);
    return 0;
  }

  fread(m_pData,sizeof(unsigned char),m_Width*m_Height*depth/8,fp);

  // Close file
  fclose(fp);

  // Success, also set FileName
  strcpy(m_pFileName,filename);

  return 1;
}

//////////////////////////////////////////////
//////////////////////////////////////////////
// FILE SAVING
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// SaveFile (dispatch function)
//********************************************
int CTexture::SaveFile(char *filename)
{
  TRACE("CTexture::SaveFile : file : %s\n",filename);
  int len = strlen(filename);
  char extension[10];
  strcpy(extension,&(filename[len-4]));

  if(extension == ".raw")
    return SaveFileRAW(filename);
  if(extension == ".bmp")
    return SaveFileBMP(filename);

  return 0;
}


//********************************************
// SaveFileRAW
//********************************************
int CTexture::SaveFileRAW(char *filename)
{
  // Check for valid image
  if((m_Width * m_Height * m_Depth) == 0)
  {
    TRACE("CTexture::SaveFileRAW : invalid image");
    return 0;
  }

  // Check for valid file
  FILE *fp = fopen(filename,"wb");

  // Try to open file
  if(!fp)
  {
    TRACE("Unable to open file for writing");
    return 0;
  }

  // Image writing
  fwrite(m_pData,sizeof(unsigned char),m_Width*m_Height*m_Depth/8,fp);

  // Close file
  fclose(fp);

  return 1;
}


//********************************************
// SaveFileBMP (*.bmp)
//********************************************
// Save windows bmp files
// Accept only 24 bits
//********************************************
int CTexture::SaveFileBMP(char *filename)
{
  if(!IsValid())
    return 0;

  if(m_Depth != 24)
    return 0;

  // Check for valid bmp file
  CFile file;
  CFileException ex;

  // Try to open file
  if(!file.Open(filename,CFile::modeCreate | CFile::modeWrite | CFile::typeBinary,&ex))
  {
    TRACE("Unable to open file for writing");
    return 0;
  }

  // File header
  BITMAPFILEHEADER FileHeader;
  WORD sign = ((WORD) ('M' << 8) | 'B');
  FileHeader.bfType = sign;
  FileHeader.bfSize = 14 + 40 + m_WidthByte32*m_Height; 
  FileHeader.bfReserved1 = 0; 
  FileHeader.bfReserved2 = 0; 
  FileHeader.bfOffBits = 54; 

  TRACE("\nSave BMP File...\n");
  TRACE("FileHeader.bfType : %d\n",FileHeader.bfType);
  TRACE("FileHeader.bfSize : %d\n",FileHeader.bfSize);
  TRACE("FileHeader.bfReserved1 : %d\n",FileHeader.bfReserved1);
  TRACE("FileHeader.bfReserved2 : %d\n",FileHeader.bfReserved2);
  TRACE("FileHeader.bfOffBits : %d\n",FileHeader.bfOffBits);

  file.Write(&FileHeader,sizeof(BITMAPFILEHEADER));

  file.Write(&m_Header,sizeof(BITMAPINFOHEADER));

  // DEBUG
  TRACE("\n");
  TRACE("IMAGE HEADER :\n");
  TRACE("biSize : %d\n",m_Header.biSize);
  TRACE("biWidth : %d\n",m_Header.biWidth);
  TRACE("biHeight : %d\n",m_Header.biHeight);
  TRACE("biPlanes : %d\n",m_Header.biPlanes);
  TRACE("biBitCount : %d\n",m_Header.biBitCount);
  TRACE("biCompression : %d\n",m_Header.biCompression);
  TRACE("biSizeImage : %d\n",m_Header.biSizeImage);
  TRACE("biXPelsPerMeter : %d\n",m_Header.biXPelsPerMeter);
  TRACE("biYPelsPerMeter : %d\n",m_Header.biYPelsPerMeter);
  TRACE("biClrUsed : %d\n",m_Header.biClrUsed);
  TRACE("biClrImportant : %d\n",m_Header.biClrImportant);

  // Image writing
  file.Write(m_pData,m_Header.biSizeImage);

  // Close file
  file.Close();

  return 1;
}


//////////////////////////////////////////////
//////////////////////////////////////////////
// CHECKING
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// IsValid
//********************************************
int CTexture::IsValid()
{
  int success = 0;
  success = (m_Depth == 24) || (m_Depth == 32);
  success &= (m_Width != 0);
  success &= (m_Height != 0);
  success &= (m_pData != NULL);
  if(!success)
  {
    TRACE("\n");
    TRACE("Invalid texture\n");
    TRACE("Width  : %d\n",m_Width);
    TRACE("Height : %d\n",m_Height);
    TRACE("Depth  : %d\n",m_Depth);
    TRACE("Data   : %x\n",m_pData);
  }
  return success;
}


//********************************************
// HigherPowerOfTwo
//********************************************
int CTexture::HigherPowerOfTwo(int value)
{
  if(value <= 0)
    return value;

  int power = 1;
  int x = 0;

  while(1)
  {
    x = (int)pow(2,power);
    if(x >= value)
      return x;
    power++;
  }
}

//********************************************
// LowerPowerOfTwo
//********************************************
int CTexture::LowerPowerOfTwo(int value)
{
  if(value <= 0)
    return value;

  int power = 1;
  int x = 0;

  while(1)
  {
    x = (int)pow(2,power);
    if(x >= value)
      return (int)pow(2,power-1);
    power++;
  }
}

//********************************************
// SameSize
//********************************************
int CTexture::SameSize(CTexture *pTexture)
{
  int success = (m_Width == pTexture->GetWidth());
  success &= (m_Height == pTexture->GetHeight());
  return success;
}


//********************************************
// Flip BGR to RGB
//********************************************
int CTexture::BGRtoRGB()
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
int CTexture::Extract(int left, int top, int right, int bottom)
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

  TRACE("Start extracting...\n");
  TRACE("New width : %d\n",NewWidth);
  TRACE("New height : %d\n",NewHeight);

  // Alloc
  unsigned char *pData = new unsigned char[NewWidthByte32*NewHeight];
  if(pData == NULL)
  {
    TRACE("Insufficiant memory");
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
int CTexture::DuplicateMirror(int left, int top, int right, int bottom)
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

  TRACE("Start duplicate mirror...\n");
  TRACE("New width : %d\n",NewWidth);
  TRACE("New widthbyte32 : %d\n",NewWidthByte32);
  TRACE("New height : %d\n",NewHeight);

  // Alloc
  unsigned char *pData = new unsigned char[NewWidthByte32*NewHeight];
  if(pData == NULL)
  {
    TRACE("Insufficiant memory");
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
int CTexture::DuplicateRepeatWidth(int left, int top, int right, int bottom)
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

  TRACE("Start duplicate repeat width...\n");
  TRACE("New width : %d\n",NewWidth);
  TRACE("New widthbyte32 : %d\n",NewWidthByte32);
  TRACE("New height : %d\n",NewHeight);

  // Alloc
  unsigned char *pData = new unsigned char[NewWidthByte32*NewHeight];
  if(pData == NULL)
  {
    TRACE("Insufficiant memory");
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
void CTexture::Fill(unsigned char r,
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
void CTexture::GreyToColor(unsigned char grey,
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
void CTexture::ColorToColor(unsigned char r1,
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
int CTexture::SetAlphaLayer(unsigned char alpha) // 0 - 255
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
int CTexture::AddAlphaLayer(unsigned char alpha) // 0 - 255
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
    TRACE("CTexture::AddAlphaLayer : insufficiant memory");
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
int CTexture::SetAlphaLayer(CTexture *pTexture) 
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
// Draw
//********************************************
int CTexture::Draw(CDC *pDC, 
									 int xOffset, 
									 int yOffset, 
									 int width, 
									 int height)
{
  // Checking
  if(!IsValid())
    return 0;

  // Flags
  if(width == -1)
    width = m_Width;
  if(height == -1)
    height = m_Height;

  // Painting
  return SetDIBitsToDevice(pDC->m_hDC, xOffset, yOffset, width, height, 0, 0, 0, 
			   m_Height, GetData(), (CONST BITMAPINFO *)&m_Header, DIB_RGB_COLORS);
}

//********************************************
// Stretch
//********************************************
int CTexture::Stretch(CDC *pDC,
											CRect *pRect)
{
	// Checking
	if(!IsValid())
		return 0;
	
	SetStretchBltMode(pDC->m_hDC,COLORONCOLOR);
	
	// Painting
	return StretchDIBits(pDC->m_hDC,
		pRect->left,
		pRect->top,
		pRect->Width(),
		pRect->Height(),
		0,
		0,
		m_Width,
		m_Height,
		GetData(),
		(CONST BITMAPINFO *)&m_Header,
		DIB_RGB_COLORS,SRCCOPY);
}

//********************************************
// ReadBuffer
//********************************************
int CTexture::ReadBuffer(unsigned char *buffer, 
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
int CTexture::ReadBufferByte32(unsigned char *pData, 
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
void CTexture::Copy(CTexture *pTexture)
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
int CTexture::ReadBuffer(float *buffer, 
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
				(BYTE)(255.0f * buffer[(width*j+i)*BytePerPixel+k]);

  return 1;
}

//********************************************
// ReadBuffer
//********************************************
int CTexture::ReadBuffer(float **ppBuffer, 
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
int CTexture::WriteBuffer(float **ppBuffer, 
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
int CTexture::WriteBuffer32(float **ppBuffer, 
														int width, 
														int height)
{
  if(ppBuffer == NULL)
    return 0;
	ASSERT(m_Depth == 32);
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
int CTexture::ReadBuffer(double *buffer, 
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
int CTexture::Grey(unsigned int x, unsigned int y)
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
void CTexture::Color(unsigned int x, unsigned int y, 
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



//***************************************
// ExportHandle
//***************************************
HANDLE CTexture::ExportHandle()
{ 
  HANDLE handle;

  TRACE("Export handle...");
  // Process source handle size
  int size = sizeof(BITMAPINFOHEADER) + m_WidthByte32 * m_Height;
  // Alloc memory
  TRACE("alloc...");
  handle = (HANDLE)::GlobalAlloc (GHND,size);
  if(handle != NULL)
  {
    char *pData = (char *) ::GlobalLock((HGLOBAL)handle);
    TRACE("lock...");
    // Copy header
    TRACE("header...");
    memcpy(pData,&m_Header,sizeof(BITMAPINFOHEADER));
    // Copy datas
    TRACE("datas...");
    memcpy(pData+sizeof(BITMAPINFOHEADER),m_pData,m_WidthByte32*m_Height);
    // Unlock
    TRACE("unlock...");
    ::GlobalUnlock((HGLOBAL)handle);
  }
  TRACE("ok\n");
	ASSERT(handle);
  return handle;
}

//********************************************
// ImportHandle
//********************************************
int CTexture::ImportHandle(HANDLE handle)
{
	TRACE("Import handle...");
  ASSERT(handle != NULL);
	char *pData = (char *) ::GlobalLock((HGLOBAL)handle);
	
	// Header
	memcpy(&m_Header,pData,sizeof(BITMAPINFOHEADER));
	
	TRACE("\n");
	TRACE("IMAGE HEADER :\n");
	TRACE("biSize : %d\n",m_Header.biSize);
	TRACE("biWidth : %d\n",m_Header.biWidth);
	TRACE("biHeight : %d\n",m_Header.biHeight);
	TRACE("biPlanes : %d\n",m_Header.biPlanes);
	TRACE("biBitCount : %d\n",m_Header.biBitCount);
	TRACE("biCompression : %d\n",m_Header.biCompression);
	TRACE("biSizeImage : %d\n",m_Header.biSizeImage);
	TRACE("biXPelsPerMeter : %d\n",m_Header.biXPelsPerMeter);
	TRACE("biYPelsPerMeter : %d\n",m_Header.biYPelsPerMeter);
	TRACE("biClrUsed : %d\n",m_Header.biClrUsed);
	TRACE("biClrImportant : %d\n",m_Header.biClrImportant);
	
	// 24 bits ?
	if(m_Header.biPlanes != 1 ||
		m_Header.biBitCount != 24)
	{
		AfxMessageBox("Texture file must have 24 bits depth");
		return 0;
	}
	
	// Alloc (does call Free before)
	Alloc(m_Header.biWidth,m_Header.biHeight,m_Header.biBitCount);
	memcpy(m_pData,pData+sizeof(BITMAPINFOHEADER),m_WidthByte32*m_Height);
	
	::GlobalUnlock((HGLOBAL)handle);
	return 1;
}


//***************************************
// GenerateGrid
//***************************************
void CTexture::GenerateGrid(unsigned int width,
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
  for(unsigned int j=0;j<height;j++)
    for(unsigned int i=0;i<width;i++)
		{
			m_pData[m_WidthByte32*j + i*3]   = rb;
			m_pData[m_WidthByte32*j + i*3+1] = gb;
			m_pData[m_WidthByte32*j + i*3+2] = bb;
		}

	// horizontal
  for(j=0;j<height;j+=size)
    for(unsigned int i=0;i<width;i++)
		{
			m_pData[m_WidthByte32*j + i*3]   = r;
			m_pData[m_WidthByte32*j + i*3+1] = g;
			m_pData[m_WidthByte32*j + i*3+2] = b;
		}

	// vertical
  for(j=0;j<height;j++)
	  for(unsigned int i=0;i<width;i+=size)
		{
			m_pData[m_WidthByte32*j + i*3]   = r;
			m_pData[m_WidthByte32*j + i*3+1] = g;
			m_pData[m_WidthByte32*j + i*3+2] = b;
		}
}

//***************************************
// GenerateLines
//***************************************
void CTexture::GenerateLines(unsigned int width,
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
  for(unsigned int j=0;j<height;j++)
    for(unsigned int i=0;i<width;i++)
		{
			m_pData[m_WidthByte32*j + i*3]   = rb;
			m_pData[m_WidthByte32*j + i*3+1] = gb;
			m_pData[m_WidthByte32*j + i*3+2] = bb;
		}

  // horizontal zebra
  for(j=0;j<height;j++)
  {
    if((j/size)%2 == 0)
      for(unsigned int i=0;i<width;i++)
      {
        m_pData[m_WidthByte32*j + i*3]   = r;
        m_pData[m_WidthByte32*j + i*3+1] = g;
        m_pData[m_WidthByte32*j + i*3+2] = b;
      }
  }
}





// ** EOF **



