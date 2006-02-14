///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CQuaternion                                                   //
//                                                                       //
//  Quaternion to represent rotations.  Each component of the            //
//  quaternion is stored as a floating point number.                     //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _QUATERNION_
#define _QUATERNION_

#include "Matrix44.h"
#include "Vector3d.h"

class CQuaternion
{

private :

  // Data
  float     m_s;
  CVector3d m_v;

public :

  // Constructors
  CQuaternion() { m_s = 1.0f; }
  CQuaternion(const float s, const float x, const float y, const float z);
  CQuaternion(const CQuaternion &quat);
  CQuaternion(const CQuaternion *pQuat);
  CQuaternion(const float s, const CVector3d *pVector);
  CQuaternion(const float s, const CVector3d &vector);
  CQuaternion(const CVector3d &vecFrom, const CVector3d &vecTo);

  virtual ~CQuaternion() { }

  void Copy(const CQuaternion &quat);
  void Copy(const CQuaternion *pQuat);

  // Data setting
  void Clear();
  void Set(const CQuaternion *pVector);
  void Set(const CQuaternion &vector);
  void Set(const float s, const CVector3d &v);
  void Set(const float s, const CVector3d *pV);
  void Set(const float s, const float x, const float y, const float z);
  void SetRotation(float ax, float ay, float az, float radAngle);

  void GetMatrix(float *mat) const;
  CMatrix44 GetMatrix() const;

  // Per coordinate (explicit inline functions)
  void s(const float s)     {	m_s = s; }
  void v(const CVector3d v) {	m_v = v; }
  void x(const float x)     {	m_v.x(x); }
  void y(const float y)     {	m_v.y(y); }
  void z(const float z)     {	m_v.z(z); } 

  // Data access (explicit inline functions)
  float     s() const { return m_s; }
  CVector3d v() const { return m_v; }
  float     x() const { return m_v.x(); }
  float     y() const { return m_v.y(); }
  float     z() const { return m_v.z(); }

  float& operator[](int i)
    { if(i==0) return (float&)m_s;
      else return (float&)m_v[i-1];}
  const float& operator[](int i) const
    { if(i==0) return (float&)m_s;
      else return (float&)m_v[i-1];}
	
  // Operators
  CQuaternion& operator+=(const CQuaternion& rQuad);
  CQuaternion& operator+=(const CQuaternion* pQuad);
  CQuaternion& operator-=(const CQuaternion& rQuad);
  CQuaternion& operator-=(const CQuaternion* pQuad);
  CQuaternion& operator*=(const float d);
  CQuaternion& operator/=(const float d)
    { return *this *= (1.f/d); }

  // Nondestructive unary negation - returns a new vector
  CQuaternion  operator -() const;

  // Binary operators
  friend CQuaternion operator+(const CQuaternion& u, const CQuaternion& v);
  friend CQuaternion operator-(const CQuaternion& u, const CQuaternion& v);
  friend CQuaternion operator*(const CQuaternion& u, const CQuaternion& v);
  friend CQuaternion operator*(const float s,        const CQuaternion& u);
  friend CQuaternion operator*(const CQuaternion& u, const float s)
    { return s * u; }
  friend CQuaternion operator/(const CQuaternion& u, const float s)
    { return (1.f/s) * u; }
  friend int         operator==(const CQuaternion& q1, const CQuaternion& q2);
  friend int         operator!=(const CQuaternion& q1, const CQuaternion& q2)
    { return !(q1 == q2); }

  int Equals(const CQuaternion& q, const float tolerence) const;


  // Misc
  double Normalize();
  double Length() const;
  double LengthSquared() const;
  void   Negate();
  CQuaternion Conjugate();
};

#endif // _QUATERNION_
