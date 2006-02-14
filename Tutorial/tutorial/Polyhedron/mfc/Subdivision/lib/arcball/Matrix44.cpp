//********************************************
// Matrix44.cpp
//********************************************
// class CMatrix44
//********************************************
// mmeyer@gg.caltech.edu
// Created  : 05/06/00
// Modified : 05/06/00
//********************************************

#include "stdafx.h"
#include <math.h>

#include "Matrix44.h"

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
CMatrix44::CMatrix44(const CMat44& m)
{
  m_data[0][0] = m[0][0];
  m_data[0][1] = m[0][1];
  m_data[0][2] = m[0][2];
  m_data[0][3] = m[0][3];
  m_data[1][0] = m[1][0];
  m_data[1][1] = m[1][1];
  m_data[1][2] = m[1][2];
  m_data[1][3] = m[1][3];
  m_data[2][0] = m[2][0];
  m_data[2][1] = m[2][1];
  m_data[2][2] = m[2][2];
  m_data[2][3] = m[2][3];
  m_data[3][0] = m[3][0];
  m_data[3][1] = m[3][1];
  m_data[3][2] = m[3][2];
  m_data[3][3] = m[3][3];
}

//********************************************
// Constructor
//********************************************
CMatrix44::CMatrix44(const double *data)
{
  Set(data);
}

//********************************************
// Constructor
//********************************************
CMatrix44::CMatrix44(double a11, double a12, double a13, double a14,
		     double a21, double a22, double a23, double a24,
		     double a31, double a32, double a33, double a34,
		     double a41, double a42, double a43, double a44)
{
  m_data[0][0] = a11;
  m_data[0][1] = a12;
  m_data[0][2] = a13;
  m_data[0][3] = a14;

  m_data[1][0] = a21;
  m_data[1][1] = a22;
  m_data[1][2] = a23;
  m_data[1][3] = a24;

  m_data[2][0] = a31;
  m_data[2][1] = a32;
  m_data[2][2] = a33;
  m_data[2][3] = a34;

  m_data[3][0] = a41;
  m_data[3][1] = a42;
  m_data[3][2] = a43;
  m_data[3][3] = a44;
}

//********************************************
// Constructor
//********************************************
CMatrix44::CMatrix44(const CMatrix44 &rMatrix)
{
  Set(&rMatrix);
}

//********************************************
// Constructor
//********************************************
CMatrix44::CMatrix44(const CMatrix44 *pMatrix)
{
  Set(pMatrix);
}


//////////////////////////////////////////////
// DATA
//////////////////////////////////////////////


//********************************************
// MakeIdentity
//  Sets the matrix to be the identity matrix.
//********************************************
void
CMatrix44::MakeIdentity()
{
  m_data[0][0] = 1.f;
  m_data[0][1] = 0.f;
  m_data[0][2] = 0.f;
  m_data[0][3] = 0.f;
  m_data[1][0] = 0.f;
  m_data[1][1] = 1.f;
  m_data[1][2] = 0.f;
  m_data[1][3] = 0.f;
  m_data[2][0] = 0.f;
  m_data[2][1] = 0.f;
  m_data[2][2] = 1.f;
  m_data[2][3] = 0.f;
  m_data[3][0] = 0.f;
  m_data[3][1] = 0.f;
  m_data[3][2] = 0.f;
  m_data[3][3] = 1.f;
}

//********************************************
// Clear
//********************************************
void
CMatrix44::Clear()
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] = 0.f;
}

//********************************************
// Set
//********************************************
void
CMatrix44::Set(const CMat44& m)
{
  m_data[0][0] = m[0][0];
  m_data[0][1] = m[0][1];
  m_data[0][2] = m[0][2];
  m_data[0][3] = m[0][3];
  m_data[1][0] = m[1][0];
  m_data[1][1] = m[1][1];
  m_data[1][2] = m[1][2];
  m_data[1][3] = m[1][3];
  m_data[2][0] = m[2][0];
  m_data[2][1] = m[2][1];
  m_data[2][2] = m[2][2];
  m_data[2][3] = m[2][3];
  m_data[3][0] = m[3][0];
  m_data[3][1] = m[3][1];
  m_data[3][2] = m[3][2];
  m_data[3][3] = m[3][3];
}

//********************************************
// Set
//********************************************
void
CMatrix44::Set(const double *data)
{
  int k = 0;
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] = data[k++];
}

//********************************************
// Set
//********************************************
void
CMatrix44::Set(const CMatrix44 &rMatrix)
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] = rMatrix.Get(i,j);
}

//********************************************
// Set
//********************************************
void
CMatrix44::Set(const CMatrix44 *pMatrix)
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] = pMatrix->Get(i,j);
}

//********************************************
// SetTranslate
//********************************************
void
CMatrix44::SetTranslate(double tx, double ty, double tz)
{
  m_data[0][0] = 1.f;
  m_data[0][1] = 0.f;
  m_data[0][2] = 0.f;
  m_data[0][3] = tx;
  m_data[1][0] = 0.f;
  m_data[1][1] = 1.f;
  m_data[1][2] = 0.f;
  m_data[1][3] = ty;
  m_data[2][0] = 0.f;
  m_data[2][1] = 0.f;
  m_data[2][2] = 1.f;
  m_data[2][3] = tz;
  m_data[3][0] = 0.f;
  m_data[3][1] = 0.f;
  m_data[3][2] = 0.f;
  m_data[3][3] = 1.f;
}

//********************************************
// SetRotate
//********************************************
void
CMatrix44::SetRotate(double ax, double ay, double az, double radAngle)
{
  double q[4];
  double c = cos(radAngle/2);
  double s = sin(radAngle/2);
  q[0] = c;
  q[1] = s * ax;
  q[2] = s * ay;
  q[3] = s * az;

  m_data[0][0] = 1.f - 2.f * (q[2] * q[2] + q[3] * q[3]);
  m_data[0][1] =     2.f * (q[1] * q[2] + q[3] * q[0]);
  m_data[0][2] =     2.f * (q[3] * q[1] - q[2] * q[0]);
  m_data[0][3] = 0.f;

  m_data[1][0] =     2.f * (q[1] * q[2] - q[3] * q[0]);
  m_data[1][1] = 1.f - 2.f * (q[3] * q[3] + q[1] * q[1]);
  m_data[1][2] =     2.f * (q[2] * q[3] + q[1] * q[0]);
  m_data[1][3] = 0.f;

  m_data[2][0] =     2.f * (q[3] * q[1] + q[2] * q[0]);
  m_data[2][1] =     2.f * (q[2] * q[3] - q[1] * q[0]);
  m_data[2][2] = 1.f - 2.f * (q[2] * q[2] + q[1] * q[1]);
  m_data[2][3] = 0.f;

  m_data[3][0] = 0.f;
  m_data[3][1] = 0.f;
  m_data[3][2] = 0.f;
  m_data[3][3] = 1.f;
}

//********************************************
// SetScale
//********************************************
void
CMatrix44::SetScale(double sx, double sy, double sz)
{
  m_data[0][0] = sx;
  m_data[0][1] = 0.f;
  m_data[0][2] = 0.f;
  m_data[0][3] = 0.f;
  m_data[1][0] = 0.f;
  m_data[1][1] = sy;
  m_data[1][2] = 0.f;
  m_data[1][3] = 0.f;
  m_data[2][0] = 0.f;
  m_data[2][1] = 0.f;
  m_data[2][2] = sz;
  m_data[2][3] = 0.f;
  m_data[3][0] = 0.f;
  m_data[3][1] = 0.f;
  m_data[3][2] = 0.f;
  m_data[3][3] = 1.f;
}


//********************************************
// Trace
//********************************************
void
CMatrix44::Trace() const
{
  TRACE("\n");
  TRACE("** Matrix **\n");
  TRACE("Address : %x\n",this);
  TRACE("Data    : (%g %g %g %g)\n",m_data[0][0],m_data[0][1],m_data[0][2],m_data[0][3]);
  TRACE("        : (%g %g %g %g)\n",m_data[1][0],m_data[1][1],m_data[1][2],m_data[1][3]);
  TRACE("        : (%g %g %g %g)\n",m_data[2][0],m_data[2][1],m_data[2][2],m_data[2][3]);
  TRACE("        : (%g %g %g %g)\n",m_data[3][0],m_data[3][1],m_data[3][2],m_data[3][3]);
}

//////////////////////////////////////////////
// OPERATORS
//////////////////////////////////////////////

//********************************************
// Operator +=
//********************************************
CMatrix44&
CMatrix44::operator+=(const CMatrix44& rMatrix)
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] += rMatrix[i][j];
  return *this;
}

//********************************************
// Operator +=
//********************************************
CMatrix44&
CMatrix44::operator+=(const CMatrix44* pMatrix)
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] += pMatrix->Get(i,j);
  return *this;
}

//********************************************
// Operator -=
//********************************************
CMatrix44&
CMatrix44::operator-=(const CMatrix44& rMatrix)
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] -= rMatrix[i][j];
  return *this;
}

//********************************************
// Operator -=
//********************************************
CMatrix44&
CMatrix44::operator-=(const CMatrix44 *pMatrix)
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] -= pMatrix->Get(i,j);
  return *this;
}

//********************************************
// Operator *=
//********************************************
CMatrix44&
CMatrix44::operator*=(const double d)
{
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      m_data[i][j] *= d;
  return *this;
}

//********************************************
// Operator + 
//********************************************
CMatrix44
operator+(const CMatrix44& u, const CMatrix44& v)
{
  CMatrix44 w;
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      w[i][j] = u[i][j] + v[i][j];
  return w;
}

//********************************************
// Operator -
//********************************************
CMatrix44
operator-(const CMatrix44& u, const CMatrix44& v)
{
  CMatrix44 w;
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      w[i][j] = u[i][j] - v[i][j];
  return w;
}

//********************************************
// Operator * 
//********************************************
CMatrix44
operator*(const double s, const CMatrix44& u)
{
  CMatrix44 w;
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      w[i][j] = s * u[i][j];
  return w;
}

//////////////////////////////////////////////
//////////////////////////////////////////////
// MATH
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// MultMatVec
// Multiply a matrix times a column vector
//********************************************
CVector3d
CMatrix44::MultMatVec(const CVector3d& v) const
{
  CVector3d res(0.f,0.f,0.f);
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
      res[i] += m_data[i][j] * v[j];
    res[i] += m_data[i][3];
  }
  double w = m_data[3][0]*v[0] + m_data[3][1]*v[1] +
            m_data[3][2]*v[2] + m_data[3][3];
  return (res/w);
}

//********************************************
// MultVecMat
// Multiply a row vector times a matrix
//********************************************
CVector3d
CMatrix44::MultVecMat(const CVector3d& v) const
{
  CVector3d res(0.f,0.f,0.f);
  for(int i = 0; i < 3; i++)
  {
    for(int j = 0; j < 3; j++)
      res[i] += v[j] * m_data[j][i];
    res[i] += m_data[3][i];
  }
  double w = v[0]*m_data[0][3] + v[1]*m_data[1][3] +
            v[2]*m_data[2][3] + m_data[3][3];
  return (res/w);
}

//********************************************
// MultMatDir
// Multiply a matrix times a column direction vector
//********************************************
CVector3d
CMatrix44::MultMatDir(const CVector3d& v) const
{
  CVector3d res(0.f,0.f,0.f);
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      res[i] += m_data[i][j] * v[j];
  return res;
}

//********************************************
// MultDirMat
// Multiply a row direction vector times a matrix
//********************************************
CVector3d
CMatrix44::MultDirMat(const CVector3d& v) const
{
  CVector3d res(0.f,0.f,0.f);
  for(int i = 0; i < 3; i++)
    for(int j = 0; j < 3; j++)
      res[i] += v[j] * m_data[j][i];
  return res;
}

//********************************************
// MultLeft
// Left multiply the matrix : mat * this
//********************************************
CMatrix44
CMatrix44::MultLeft(const CMatrix44& mat) const
{
  CMatrix44 res;
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      for(int k = 0; k < 4; k++)
	res[i][j] += mat[i][k] * m_data[k][j];
  return res;
}

//********************************************
// MultRight
// Right multiply the matrix : this * mat
//********************************************
CMatrix44
CMatrix44::MultRight(const CMatrix44& mat) const
{
  CMatrix44 res;
  for(int i = 0; i < 4; i++)
    for(int j = 0; j < 4; j++)
      for(int k = 0; k < 4; k++)
	res[i][j] += m_data[i][k] * mat[k][j];
  return res;
}

//////////////////////////////////////////////
//////////////////////////////////////////////
// PROCESSING
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// Det4
//********************************************
double
CMatrix44::Det4() const
{
  return (m_data[0][0] * Det3(1,2,3,1,2,3)
	  - m_data[0][1] * Det3(1,2,3,0,2,3)
	  + m_data[0][2] * Det3(1,2,3,0,1,3)
	  - m_data[0][3] * Det3(1,2,3,0,1,2));
}

//********************************************
// Det3
//********************************************
double
CMatrix44::Det3(int r1, int r2, int r3, int c1, int c2, int c3) const
{
  return (m_data[r1][c1] * Det2( r2, r3, c2, c3 )
          - m_data[r1][c2] * Det2( r2, r3, c1, c3 )
          + m_data[r1][c3] * Det2( r2, r3, c1, c2 ));
}

//********************************************
// Transpose
// Returns the transpose of the 4x4 matrix.
//********************************************
CMatrix44
CMatrix44::Transpose() const
{
  return CMatrix44(m_data[0][0], m_data[1][0], m_data[2][0], m_data[3][0],
		   m_data[0][1], m_data[1][1], m_data[2][1], m_data[3][1],
		   m_data[0][2], m_data[1][2], m_data[2][2], m_data[3][2],
		   m_data[0][3], m_data[1][3], m_data[2][3], m_data[3][3]);
}

//********************************************
// Adjoint
// Returns the adjoint of the 4x4 matrix.
// Adjoint_ij = (-1)^(i+j) * alpha_ji
//  where alpha_ij is the determinant of the 
//  submatrix of A without row i and column j
//********************************************
CMatrix44
CMatrix44::Adjoint() const
{
  CMatrix44 a;

  a[0][0] =  Det3(1,2,3,1,2,3);
  a[0][1] = -Det3(0,2,3,1,2,3);
  a[0][2] =  Det3(0,1,3,1,2,3);
  a[0][3] = -Det3(0,1,2,1,2,3);

  a[1][0] = -Det3(1,2,3,0,2,3);
  a[1][1] =  Det3(0,2,3,0,2,3);
  a[1][2] = -Det3(0,1,3,0,2,3);
  a[1][3] =  Det3(0,1,2,0,2,3);

  a[2][0] =  Det3(1,2,3,0,1,3);
  a[2][1] = -Det3(0,2,3,0,1,3);
  a[2][2] =  Det3(0,1,3,0,1,3);
  a[2][3] = -Det3(0,1,2,0,1,3);

  a[3][0] = -Det3(1,2,3,0,1,2);
  a[3][1] =  Det3(0,2,3,0,1,2);
  a[3][2] = -Det3(0,1,3,0,1,2);
  a[3][3] =  Det3(0,1,2,0,1,2);

  return a;
}

//********************************************
// Inverse
// Returns the inverse of the 4x4 matrix.
// A^-1 = Adjoint(A) / Determinant(A)
//********************************************
CMatrix44
CMatrix44::Inverse() const
{
  double det = Determinant();
  CMatrix44 mat = Adjoint();
  mat /= det;
  return mat;
}

//********************************************
// Inverse
// Returns the identity matrix.
//********************************************
CMatrix44
CMatrix44::Identity()
{
  return CMatrix44(1.f, 0.f, 0.f, 0.f,
		   0.f, 1.f, 0.f, 0.f,
		   0.f, 0.f, 1.f, 0.f,
		   0.f, 0.f, 0.f, 1.f);
}

// ** EOF **
