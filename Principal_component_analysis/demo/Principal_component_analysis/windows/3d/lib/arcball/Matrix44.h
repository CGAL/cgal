///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CMatrix44                                                     //
//                                                                       //
//  4x4 matrix to represent transformations.                             //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _MATRIX_44_
#define _MATRIX_44_

typedef double CMat44[4][4];

#include "Vector3d.h"

class CMatrix44
{

private :

  // Data
  CMat44 m_data;

public :

  // Constructors
  CMatrix44() { }
  CMatrix44(const CMat44& m);
  CMatrix44(const double *data);
  CMatrix44(double a11, double a12, double a13, double a14,
	    double a21, double a22, double a23, double a24,
	    double a31, double a32, double a33, double a34,
	    double a41, double a42, double a43, double a44);
  CMatrix44(const CMatrix44 &rMatrix);
  CMatrix44(const CMatrix44 *pMatrix);

  virtual ~CMatrix44() { }

  // Debug
  void Trace() const;

  // Data setting
  void Clear();
  void MakeIdentity();
  void Set(const CMat44& m);
  void Set(const double *data);
  void Set(const CMatrix44 &rMatrix);
  void Set(const CMatrix44 *pMatrix);

  void SetTranslate(double tx, double ty, double tz);
  void SetRotate(double ax, double ay, double az, double radAngle);
  void SetScale(double sx, double sy, double sz);

  // Per element (explicit inline functions)
  void Set(int i, int j, double data) { m_data[i][j] = data; }

  // Make it look like a usual matrix (so you can do m[3][2])
  double*       operator[](int i)       { return &m_data[i][0]; }
  const double* operator[](int i) const { return &m_data[i][0]; }

  // Data access (explicit inline functions)
  double  Get(int i, int j) const { return m_data[i][j]; }
  double* GetData()               { return &m_data[0][0]; }

  // Operators
  CMatrix44& operator+=(const CMatrix44& rMatrix);
  CMatrix44& operator+=(const CMatrix44* pMatrix);
  CMatrix44& operator-=(const CMatrix44& rMatrix);
  CMatrix44& operator-=(const CMatrix44* pMatrix);
  CMatrix44& operator*=(const double d);
  CMatrix44& operator/=(const double d)
    { return *this *= (1.f/d); }

  // Binary operators
  friend CMatrix44 operator+(const CMatrix44& u, const CMatrix44& v);
  friend CMatrix44 operator-(const CMatrix44& u, const CMatrix44& v);
  friend CMatrix44 operator*(const double s,      const CMatrix44& u);
  friend CMatrix44 operator*(const CMatrix44& u, const double s)
    { return s * u; }
  friend CMatrix44 operator/(const CMatrix44& u, const double s)
    { return (1.f/s) * u; }

  friend CVector3d operator*(const CMatrix44& m, const CVector3d& v)
    { return m.MultMatVec(v); }
  friend CVector3d operator*(const CVector3d& v, const CMatrix44& m)
    { return m.MultVecMat(v); }
  friend CMatrix44 operator*(const CMatrix44& u, const CMatrix44& v)
     { return u.MultLeft(v); }

  // Math
  CVector3d MultMatVec(const CVector3d& v) const;
  CVector3d MultVecMat(const CVector3d& v) const;
  CVector3d MultMatDir(const CVector3d& v) const;
  CVector3d MultDirMat(const CVector3d& v) const;
  CMatrix44 MultLeft(const CMatrix44& mat) const;
  CMatrix44 MultRight(const CMatrix44& mat) const;

  // Misc
  double Determinant() const { return Det4(); }
  double Det4() const;
  double Det3(int r1, int r2, int r3, int c1, int c2, int c3) const;
  double Det2(int r1, int r2, int c1, int c2) const
   { return (m_data[r1][c1]*m_data[r2][c2] - m_data[r2][c1]*m_data[r1][c2]); }
  CMatrix44 Transpose() const;
  CMatrix44 Adjoint() const;
  CMatrix44 Inverse() const;

  static CMatrix44  Identity();
};


#endif // _MATRIX_44_
