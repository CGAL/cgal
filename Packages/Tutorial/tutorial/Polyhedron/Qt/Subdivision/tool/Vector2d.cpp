//********************************************
// Vector2d.cpp pouet
//********************************************
// class CVector2d
//********************************************
// mmeyer@gg.caltech.edu
// Created :  09/07/00
// Modified : Pierre Alliez 21/02/02
//********************************************

#include "Vector2d.h"

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
//  Get the 2 components of the vector
//********************************************
void
CVector2d::Get(double& x, double& y) const
{
  x = vec[0];
  y = vec[1];
}

//********************************************
// Trace
//********************************************
void
CVector2d::Trace() const
{
  fprintf(stderr,"\n");
  fprintf(stderr,"** Vector **\n");
  fprintf(stderr,"Coordinates : (%g %g)\n",vec[0],vec[1]);
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
CVector2d::Normalize()
{
  double len = Length();
  if(len != 0.0f)
    (*this) *= (1.0/len);
  else
    Set(0.f,0.f);

  return len;
}

//********************************************
// Normalize
//********************************************
double
CVector2d::Normalize(double value)
{
  double len = Length();
  if(len != 0.0f)
    (*this) *= (value/len);
  else
    Set(0.f,0.f);

  return len;
}

//********************************************
// LengthSquared
//********************************************
double
CVector2d::LengthSquared() const
{
  return ( (double)vec[0]*(double)vec[0]
	 + (double)vec[1]*(double)vec[1]);
}
	
//********************************************
// Length
//********************************************
double
CVector2d::Length() const
{
  return sqrt( (double)vec[0]*(double)vec[0]
	     + (double)vec[1]*(double)vec[1]);
}
	
//********************************************
// IsCollinear
//********************************************
int
CVector2d::IsCollinear(CVector2d *pVector) const
{
  double x = pVector->x() / vec[0];
  double y = pVector->y() / vec[1];
  return (x == y);
}

//********************************************
// IsCollinear
//********************************************
int
CVector2d::IsCollinear(CVector2d &vector) const
{
  double x = vector.x() / vec[0];
  double y = vector.y() / vec[1];
  return (x == y);
}

//********************************************
// Negate
//********************************************
void
CVector2d::Negate()
{
  vec[0] = -vec[0];
  vec[1] = -vec[1];
}

//////////////////////////////////////////////
// OPERATORS
//////////////////////////////////////////////



//********************************************
// Operator +=
//********************************************
CVector2d&
CVector2d::operator+=(const CVector2d* pVector)
{
  vec[0] += pVector->x();
  vec[1] += pVector->y();
  return *this;
}


//********************************************
// Operator -=
//********************************************
CVector2d&
CVector2d::operator-=(const CVector2d* pVector)
{
  vec[0] -= pVector->x();
  vec[1] -= pVector->y();
  return *this;
}


//********************************************
// Operator -
//  Nondestructive unary -
//  Returns a new vector.
//********************************************
CVector2d
CVector2d::operator -() const
{
  return CVector2d(-vec[0],-vec[1]);
}


//********************************************
// Operator ==
//********************************************
int
operator==(const CVector2d& v1, const CVector2d& v2)
{
  return (v1.vec[0] == v2.vec[0] &&
	  v1.vec[1] == v2.vec[1]);
}

//********************************************
// Equals
//  Determines if two vectors are equal
//  within a tolerence (squared distance).
//********************************************
int
CVector2d::Equals(const CVector2d& v, double tolerence) const
{
  CVector2d diff = *this - v;

  return diff.LengthSquared() <= tolerence;
}


//********************************************
// Operator + 
//********************************************
CVector2d
operator+(const CVector2d& u, const CVector2d& v)
{
  return CVector2d(u.vec[0]+v.vec[0],u.vec[1]+v.vec[1]);
}

//********************************************
// Operator -
//********************************************
CVector2d
operator-(const CVector2d& u, const CVector2d& v)
{
  return CVector2d(u.vec[0]-v.vec[0],u.vec[1]-v.vec[1]);
}

//********************************************
// Operator * 
//********************************************
CVector2d
operator*(double s, const CVector2d& u)
{
  return CVector2d(u.vec[0] * s, u.vec[1] * s);
}


// ** EOF **
