///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CQuaternion                                                   //
//                                                                       //
//  Quaternion to represent rotations.  Each component of the            //
//  quaternion is stored as a doubleing point number.                     //
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
  double    m_s;
  CVector3d m_v;

public :

  // Constructors
  CQuaternion() { m_s = 1.0f; }
  CQuaternion(const double s, const double x, const double y, const double z);
  CQuaternion(const CQuaternion &quat);
  CQuaternion(const CQuaternion *pQuat);
  CQuaternion(const double s, const CVector3d *pVector);
  CQuaternion(const double s, const CVector3d &vector);
  CQuaternion(const CVector3d &vecFrom, const CVector3d &vecTo);

  virtual ~CQuaternion() { }

  // Debug
  void Trace() const;

  // Data setting
  void Clear();
  void Set(const CQuaternion *pVector);
  void Set(const CQuaternion &vector);
  void Set(const double s, const CVector3d &v);
  void Set(const double s, const CVector3d *pV);
  void Set(const double s, const double x, const double y, const double z);
  void SetRotation(double ax, double ay, double az, double radAngle);

  void GetMatrix(double *mat) const;
  CMatrix44 GetMatrix() const;

  // Per coordinate (explicit inline functions)
  void s(const double s)     {	m_s = s; }
  void v(const CVector3d v) {	m_v = v; }
  void x(const double x)     {	m_v.x(x); }
  void y(const double y)     {	m_v.y(y); }
  void z(const double z)     {	m_v.z(z); } 

  // Data access (explicit inline functions)
  double     s() const { return m_s; }
  CVector3d v() const { return m_v; }
  double     x() const { return m_v.x(); }
  double     y() const { return m_v.y(); }
  double     z() const { return m_v.z(); }

  double& operator[](int i)
    { if(i==0) return (double&)m_s;
      else return (double&)m_v[i-1];}
  const double& operator[](int i) const
    { if(i==0) return (double&)m_s;
      else return (double&)m_v[i-1];}
	
  // Operators
  CQuaternion& operator+=(const CQuaternion& rQuad);
  CQuaternion& operator+=(const CQuaternion* pQuad);
  CQuaternion& operator-=(const CQuaternion& rQuad);
  CQuaternion& operator-=(const CQuaternion* pQuad);
  CQuaternion& operator*=(const double d);
  CQuaternion& operator/=(const double d)
    { return *this *= (1.f/d); }

  // Nondestructive unary negation - returns a new vector
  CQuaternion  operator -() const;

  // Binary operators
  friend CQuaternion operator+(const CQuaternion& u, const CQuaternion& v);
  friend CQuaternion operator-(const CQuaternion& u, const CQuaternion& v);
  friend CQuaternion operator*(const CQuaternion& u, const CQuaternion& v);
  friend CQuaternion operator*(const double s,        const CQuaternion& u);
  friend CQuaternion operator*(const CQuaternion& u, const double s)
    { return s * u; }
  friend CQuaternion operator/(const CQuaternion& u, const double s)
    { return (1.f/s) * u; }
  friend int         operator==(const CQuaternion& q1, const CQuaternion& q2);
  friend int         operator!=(const CQuaternion& q1, const CQuaternion& q2)
    { return !(q1 == q2); }

  int Equals(const CQuaternion& q, const double tolerence) const;


  // Misc
  double Normalize();
  double Length() const;
  double LengthSquared() const;
  void   Negate();
  CQuaternion Conjugate();
};

#endif // _QUATERNION_
