//********************************************
// Quaternion.cpp
//********************************************
// class CQuaternion
//********************************************
// mmeyer@gg.caltech.edu
// Created  : 09/07/00
// Modified : 09/07/00
//********************************************

#include "stdafx.h"
#include <math.h>

#include "Quaternion.h"

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
CQuaternion::CQuaternion(const double s, 
												 const double x, 
												 const double y, 
												 const double z)
{
  m_s = s;
  m_v.Set(x,y,z);
}

//********************************************
// Constructor
//********************************************
CQuaternion::CQuaternion(const CQuaternion &quat)
{
  Set(quat);
}

//********************************************
// Constructor
//********************************************
CQuaternion::CQuaternion(const CQuaternion *pQuat)
{
  Set(pQuat);
}

//********************************************
// Constructor
//********************************************
CQuaternion::CQuaternion(const double s, 
												 const CVector3d *pVector)
{
  m_s = s;
  m_v.Set(pVector);
}

//********************************************
// Constructor
//********************************************
CQuaternion::CQuaternion(const double s, 
												 const CVector3d& vector)
{
  m_s = s;
  m_v.Set(vector);
}

//********************************************
// Constructor
// This assumes unit vecFrom and vecTo
//********************************************
CQuaternion::CQuaternion(const CVector3d &vecFrom, 
												 const CVector3d &vecTo)
{
  CVector3d& vecHalf = vecTo + vecFrom;
  vecHalf.Normalize();
  m_s = vecHalf.Dot(vecTo);
  m_v = vecHalf.Cross(vecTo);
}


//////////////////////////////////////////////
// DATA
//////////////////////////////////////////////


//********************************************
// Clear
// Clears to a unit quaternion (no rotation)
//********************************************
void
CQuaternion::Clear()
{
  Set(1.0f,0.0f,0.0f,0.0f);
}

//********************************************
// Set
//********************************************
void
CQuaternion::Set(const double s, 
								 const double x, 
								 const double y, 
								 const double z)
{
  m_s = s;
  m_v.Set(x,y,z);
}

//********************************************
// Set
//********************************************
void
CQuaternion::Set(const double s, 
								 const CVector3d &v)
{
  m_s = s;
  m_v.Set(v);
}

//********************************************
// Set
//********************************************
void
CQuaternion::Set(const double s, 
								 const CVector3d *pV)
{
  m_s = s;
  m_v.Set(pV);
}

//********************************************
// Set
//********************************************
void
CQuaternion::Set(const CQuaternion *pQuat)
{
  Set(pQuat->s(),pQuat->v());
}

//********************************************
// Set
//********************************************
void
CQuaternion::Set(const CQuaternion &quat)
{
  Set(quat.s(),quat.v());
}

//********************************************
// SetRotation
//********************************************
void
CQuaternion::SetRotation(double ax, 
												 double ay, 
												 double az, 
												 double radAngle)
{
  double halfAngle = radAngle / 2.f;
  m_s = cos(halfAngle);
  m_v.Set(ax,ay,az);
  m_v *= sin(halfAngle);
}

//********************************************
// GetMatrix
// Construct rotation matrix from (possibly non-unit) quaternion.
// Assumes matrix is used to multiply column vector on the left:
// Vnew = Matrix * Vold.  Works correctly for right-handed 
// coordinate systems and right-handed rotations.
//********************************************
void
CQuaternion::GetMatrix(double *mat) const
{
  double Nq = LengthSquared();
  double c  = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
  double xc = x()*c,	      yc = y()*c,	  zc = z()*c;
  double sx = s()*xc,	      sy = s()*yc,	  sz = s()*zc;
  double xx = x()*xc,	      xy = x()*yc,	  xz = x()*zc;
  double yy = y()*yc,	      yz = y()*zc,	  zz = z()*zc;
  mat[0]  = double(1.0 - (yy + zz));
  mat[1]  = double(xy - sz);
  mat[2]  = double(xz + sy);
  mat[4]  = double(xy + sz);
  mat[5]  = double(1.0 - (xx + zz));
  mat[6]  = double(yz - sx);
  mat[8]  = double(xz - sy);
  mat[9]  = double(yz + sx);
  mat[10] = double(1.0 - (xx + yy));
  mat[3]  = mat[7] = mat[11] = mat[12] = mat[13] = mat[14] = 0.0f;
  mat[15] = 1.0f;
}

//********************************************
// GetMatrix
// Construct rotation matrix from (possibly non-unit) quaternion.
// Assumes matrix is used to multiply column vector on the left:
// Vnew = Matrix * Vold.  Works correctly for right-handed 
// coordinate systems and right-handed rotations.
//********************************************
CMatrix44
CQuaternion::GetMatrix() const
{
  CMatrix44 mat;
  double Nq = LengthSquared();
  double c  = (Nq > 0.0) ? (2.0 / Nq) : 0.0;
  double xc = x()*c,	      yc = y()*c,	  zc = z()*c;
  double sx = s()*xc,	      sy = s()*yc,	  sz = s()*zc;
  double xx = x()*xc,	      xy = x()*yc,	  xz = x()*zc;
  double yy = y()*yc,	      yz = y()*zc,	  zz = z()*zc;
  mat[0][0] = double(1.0 - (yy + zz));
  mat[0][1] = double(xy - sz);
  mat[0][2] = double(xz + sy);
  mat[1][0] = double(xy + sz);
  mat[1][1] = double(1.0 - (xx + zz));
  mat[1][2] = double(yz - sx);
  mat[2][0] = double(xz - sy);
  mat[2][1] = double(yz + sx);
  mat[2][2] = double(1.0 - (xx + yy));
  mat[0][3] = mat[1][3] = mat[2][3] = mat[3][0] = mat[3][1] = mat[3][2] = 0.0f;
  mat[3][3] = 1.0f;
  return mat;
}


//********************************************
// Trace
//********************************************
void
CQuaternion::Trace() const
{
  TRACE("\n");
  TRACE("** Quaternion **\n");
  TRACE("Address      : %x\n",this);
  TRACE("Coordinates : (%g %g %g %g)\n",m_s,m_v.x(),m_v.y(),m_v.z());
}

//////////////////////////////////////////////
// OPERATORS
//////////////////////////////////////////////

//********************************************
// Operator +=
//********************************************
CQuaternion&
CQuaternion::operator+=(const CQuaternion& rQuad)
{
  m_s += rQuad.s();
  m_v += rQuad.v();
  return *this;
}

//********************************************
// Operator +=
//********************************************
CQuaternion&
CQuaternion::operator+=(const CQuaternion* pQuad)
{
  m_s += pQuad->s();
  m_v += pQuad->m_v;
  return *this;
}

//********************************************
// Operator -=
//********************************************
CQuaternion&
CQuaternion::operator-=(const CQuaternion& rQuad)
{
  m_s -= rQuad.s();
  m_v -= rQuad.v();
  return *this;
}

//********************************************
// Operator -=
//********************************************
CQuaternion&
CQuaternion::operator-=(const CQuaternion* pQuad)
{
  m_s -= pQuad->s();
  m_v -= pQuad->v();
  return *this;
}

//********************************************
// Operator *=
//********************************************
CQuaternion&
CQuaternion::operator*=(const double d)
{
  m_s *= d;
  m_v *= d;
  return *this;
}

//********************************************
// Operator -
//  Nondestructive unary -
//  Returns a new vector.
//********************************************
CQuaternion
CQuaternion::operator -() const
{
  return CQuaternion(-m_s,-m_v);
}

//********************************************
// Operator + 
//********************************************
CQuaternion
operator+(const CQuaternion& u, const CQuaternion& v)
{
  return CQuaternion(u.m_s+v.m_s,u.m_v+v.m_v);
}

//********************************************
// Operator -
//********************************************
CQuaternion
operator-(const CQuaternion& u, const CQuaternion& v)
{
  return CQuaternion(u.m_s-v.m_s,u.m_v-v.m_v);
}

//********************************************
// Operator *
//********************************************
CQuaternion
operator*(const CQuaternion& u, const CQuaternion& v)
{
  CQuaternion w;
  CVector3d a  = u.m_v;
  CVector3d b  = v.m_v;
  double     ws = u.m_s * v.m_s - a.Dot(b);
  CVector3d wv = u.m_s*b + v.m_s*a + a.Cross(b);
  w.Set(ws,wv);
  return w;
}

//********************************************
// Operator * 
//********************************************
CQuaternion
operator*(const double s, 
					const CQuaternion& u)
{
  return CQuaternion(u.m_s * s, u.m_v * s);
}

//********************************************
// Operator ==
//********************************************
int
operator==(const CQuaternion& q1, 
					 const CQuaternion& q2)
{
  return (q1.m_s == q2.m_s &&
	  q1.m_v == q2.m_v);
}

//********************************************
// Equals
//  Determines if two quaternions are equal
//  within a tolerence (squared distance).
//********************************************
int
CQuaternion::Equals(const CQuaternion& q, 
										const double tolerence) const
{
  CQuaternion diff = *this - q;

  return diff.LengthSquared() <= tolerence;
}

//////////////////////////////////////////////
//////////////////////////////////////////////
// PROCESSING
//////////////////////////////////////////////
//////////////////////////////////////////////

//********************************************
// Normalize
//********************************************
double
CQuaternion::Normalize()
{
  double len = Length();
  if(len != 0.0f)
    (*this) *= (1.0/len);
  else
    Set(0.f,0.f,0.f,0.f);

  return len;
}

//********************************************
// Length
//********************************************
double
CQuaternion::Length()const
{
  return sqrt((double)m_s *(double)m_s + 
	      (double)m_v.x()*(double)m_v.x() + 
	      (double)m_v.y()*(double)m_v.y() + 
	      (double)m_v.z()*(double)m_v.z());
}
	
//********************************************
// LengthSquared
//********************************************
double
CQuaternion::LengthSquared()const
{
  return ((double)m_s    *(double)m_s + 
	  (double)m_v.x()*(double)m_v.x() + 
	  (double)m_v.y()*(double)m_v.y() + 
	  (double)m_v.z()*(double)m_v.z());
}

//********************************************
// Negate
//  Negate each component of the quaternion
//********************************************
void
CQuaternion::Negate()
{
  m_s = -m_s;
  m_v = -m_v;
}

//********************************************
// Conjugate
//  Return the conjugate of the quaternion
//********************************************
CQuaternion
CQuaternion::Conjugate()
{
  return CQuaternion(m_s,-m_v[0],-m_v[1],-m_v[2]);
}

// ** EOF **



