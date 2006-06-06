//********************************************
// Vector3d.cpp
//********************************************
// class CVector3d
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/07/00
// Modified : 09/07/00
//********************************************

#include "stdafx.h"
#include "Vector3d.h"

#include <math.h>
#include <stdio.h>


//////////////////////////////////////////////
// CONSTRUCTORS
//////////////////////////////////////////////


//////////////////////////////////////////////
// DATA
//////////////////////////////////////////////

//********************************************
// Get
//  Get the 3 components of the vector
//********************************************
void
CVector3d::Get(double& x, double& y, double& z) const
{
  x = vec[0];
  y = vec[1];
  z = vec[2];
}

//********************************************
// Trace
//********************************************
void
CVector3d::Trace() const
{
  TRACE("\n");
  TRACE("** Vector **\n");
  TRACE("Address      : %x\n",this);
  TRACE("Coordinates : (%g %g %g)\n",vec[0],vec[1],vec[2]);
}

//////////////////////////////////////////////
// OPERATORS
//////////////////////////////////////////////

//********************************************
// Operator +=
//********************************************
CVector3d&
CVector3d::operator+=(const CVector3d& rVector)
{
  vec[0] += rVector.x();
  vec[1] += rVector.y();
  vec[2] += rVector.z();
  return *this;
}

//********************************************
// Operator +=
//********************************************
CVector3d&
CVector3d::operator+=(const CVector3d* pVector)
{
  vec[0] += pVector->x();
  vec[1] += pVector->y();
  vec[2] += pVector->z();
  return *this;
}

//********************************************
// Operator -=
//********************************************
CVector3d&
CVector3d::operator-=(const CVector3d& rVector)
{
  vec[0] -= rVector.x();
  vec[1] -= rVector.y();
  vec[2] -= rVector.z();
  return *this;
}

//********************************************
// Operator -=
//********************************************
CVector3d&
CVector3d::operator-=(const CVector3d* pVector)
{
  vec[0] -= pVector->x();
  vec[1] -= pVector->y();
  vec[2] -= pVector->z();
  return *this;
}

//********************************************
// Operator *=
//********************************************
CVector3d&
CVector3d::operator*=(double d)
{
  vec[0] *= d;
  vec[1] *= d;
  vec[2] *= d;
  return *this;
}

//********************************************
// Operator -
//  Nondestructive unary -
//  Returns a new vector.
//********************************************
CVector3d
CVector3d::operator -() const
{
  return CVector3d(-vec[0],-vec[1],-vec[2]);
}

//********************************************
// Operator + 
//********************************************
CVector3d
operator+(const CVector3d& u, const CVector3d& v)
{
  return CVector3d(u.vec[0]+v.vec[0],u.vec[1]+v.vec[1],u.vec[2]+v.vec[2]);
}

//********************************************
// Operator -
//********************************************
CVector3d
operator-(const CVector3d& u, const CVector3d& v)
{
  return CVector3d(u.vec[0]-v.vec[0],u.vec[1]-v.vec[1],u.vec[2]-v.vec[2]);
}

//********************************************
// Operator * 
//********************************************
CVector3d
operator*(double s, const CVector3d& u)
{
  return CVector3d(u.vec[0] * s, u.vec[1] * s, u.vec[2] * s);
}

//********************************************
// Operator ^
//  Returns the cross product of u and v.
//********************************************
CVector3d
operator^(const CVector3d& u, const CVector3d& v)
{
  return CVector3d(u.vec[1] * v.vec[2] - u.vec[2] * v.vec[1],
		   u.vec[2] * v.vec[0] - u.vec[0] * v.vec[2],
		   u.vec[0] * v.vec[1] - u.vec[1] * v.vec[0]);
}

//********************************************
// Operator ==
//********************************************
int
operator==(const CVector3d& v1, const CVector3d& v2)
{
  return (v1.vec[0] == v2.vec[0] &&
	  v1.vec[1] == v2.vec[1] &&
	  v1.vec[2] == v2.vec[2]);
}

//********************************************
// Equals
//  Determines if two vectors are equal
//  within a tolerence (squared distance).
//********************************************
int
CVector3d::Equals(const CVector3d& v, double tolerence) const
{
  CVector3d diff = *this - v;

  return diff.LengthSquared() <= tolerence;
}


//////////////////////////////////////////////
//////////////////////////////////////////////
// PROCESSING
//////////////////////////////////////////////
//////////////////////////////////////////////


//********************************************
// Cross
//********************************************
CVector3d
CVector3d::Cross(const CVector3d& v) const
{
  return CVector3d(y() * v.z() - z() * v.y(),
		   z() * v.x() - x() * v.z(),
		   x() * v.y() - y() * v.x());
}

//********************************************
// Cross
//********************************************
CVector3d
CVector3d::Cross(const CVector3d* pV) const
{
  return CVector3d(y() * pV->z() - z() * pV->y(),
		   z() * pV->x() - x() * pV->z(),
		   x() * pV->y() - y() * pV->x());
}


//********************************************
// Normalize
//********************************************
double
CVector3d::Normalize()
{
  double len = Length();
  if(len != 0.0f)
    (*this) *= (1.0/len);
  else
    Set(0.f,0.f,0.f);

  return len;
}

//********************************************
// Normalize
//********************************************
double
CVector3d::Normalize(double value)
{
  double len = Length();
  if(len != 0.0f)
    (*this) *= (value/len);
  else
    Set(0.f,0.f,0.f);

  return len;
}

//********************************************
// LengthSquared
//********************************************
double
CVector3d::LengthSquared() const
{
  return ( (double)vec[0]*(double)vec[0]
	 + (double)vec[1]*(double)vec[1]
	 + (double)vec[2]*(double)vec[2]);
}
	
//********************************************
// Length
//********************************************
double
CVector3d::Length() const
{
  return sqrt( (double)vec[0]*(double)vec[0]
	     + (double)vec[1]*(double)vec[1]
	     + (double)vec[2]*(double)vec[2]);
}
	
//********************************************
// IsCollinear
//********************************************
int
CVector3d::IsCollinear(CVector3d *pVector) const
{
  double x = pVector->x() / vec[0];
  double y = pVector->y() / vec[1];
  double z = pVector->z() / vec[2];
  return ((x == y) && (y == z));
}

//********************************************
// IsCollinear
//********************************************
int
CVector3d::IsCollinear(CVector3d &vector) const
{
  double x = vector.x() / vec[0];
  double y = vector.y() / vec[1];
  double z = vector.z() / vec[2];
  return ((x == y) && (y == z));
}

//********************************************
// Negate
//********************************************
void
CVector3d::Negate()
{
  vec[0] = -vec[0];
  vec[1] = -vec[1];
  vec[2] = -vec[2];
}

//********************************************
// Rotate this vector around pAround by angle
// by Haeyoung Lee
//********************************************
CVector3d CVector3d::Rotate(double angle, 
														CVector3d Around) 
{
	double f1, f2, f3;
	CVector3d t1, t2;
	
	
	f1 = (double)cos((double)angle);
	f2 = (double)sin((double)angle);
	t1 = Projection(&Around);
	t2 = Around.Cross(this);
	f3 = Dot(Around);
	
	return CVector3d((double)(f1*t1.x()+f2*t2.x()+f3*Around.x()),
		(double)(f1*t1.y()+f2*t2.y()+f3*Around.y()),
		(double)(f1*t1.z()+f2*t2.z()+f3*Around.z()));
	
}

//********************************************
// Projection
// by Haeyoung Lee
//********************************************
CVector3d    
CVector3d::Projection(const CVector3d* pV) const
{
  double alpha = Dot(pV)/pV->Dot(pV);
	return CVector3d(x()-alpha* pV->x(), 
		               y()-alpha*pV->y(),
		               z()-alpha*pV->z());
}





// ** EOF **
