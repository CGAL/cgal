///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CVector2d                                                     //
//                                                                       //
//  2D Vector to represent points or directions.  Each component of the  //
//  vector is stored as a doubleing point number.                         //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _VECTOR2D_
#define _VECTOR2D_

class CVector2d
{

private :

protected:

 double vec[2];  // Storage for the vector components

public :

  // Constructors
  CVector2d() { vec[0] = vec[1] = 0; }
  CVector2d(const double x, const double y)
    { vec[0] = x; vec[1] = y; }
  CVector2d(const double v[2])
    { vec[0] = v[0]; vec[1] = v[1]; }
  CVector2d(const CVector2d &vector)  { Set(vector); }
  CVector2d(const CVector2d *pVector) { Set(pVector); }

  virtual ~CVector2d() { }

  // Debug
  void Trace() const;

  // Data setting
  void Clear()
    { vec[0] = 0.f; vec[1] = 0.f; }
  void Set(const CVector2d *pVector) { Set(pVector->GetArray()); }
  void Set(const CVector2d& vector)  { Set(vector.GetArray()); }
  void Set(const double x, const double y)
    { vec[0] = x; vec[1] = y; }
  void Set(const double v[2])
    { vec[0] = v[0]; vec[1] = v[1]; }

  // Data Access
  const double* GetArray() const { return vec; }
  void         Get(double& x, double& y) const;

  // Per coordinate (explicit inline functions)
  void x(double newX) { vec[0] = newX; }
  void y(double newY) { vec[1] = newY; }

  // Data access (explicit inline functions)
  double x() const { return (vec[0]); }
  double y() const { return (vec[1]); }

  // Data access using indices
  double&       operator[] (int i)       { return (vec[i]); }
  const double& operator[] (int i) const { return (vec[i]); }

  // Operators
  inline CVector2d& operator+=(const CVector2d& rVector);
  inline CVector2d& operator+=(const CVector2d* pVector);
  inline CVector2d& operator-=(const CVector2d& rVector);
  inline CVector2d& operator-=(const CVector2d* pVector);
  inline CVector2d& operator*=(double d);
  inline CVector2d& operator/=(double d) { return *this *= (1.f/d); }

  // Nondestructive unary negation - returns a new vector
  CVector2d  operator -() const;

  // Binary operators
  friend CVector2d operator+(const CVector2d& u, const CVector2d& v);
  friend CVector2d operator-(const CVector2d& u, const CVector2d& v);
  friend CVector2d operator*(double s,            const CVector2d& u);
  friend CVector2d operator*(const CVector2d& u, double s)
    { return s * u; }
  friend CVector2d operator/(const CVector2d& u, double s)
    { return (1.f/s) * u; }
  friend int operator==(const CVector2d& v1, const CVector2d& v2);
  friend int operator!=(const CVector2d& v1, const CVector2d& v2)
    { return !(v1 == v2); }

  int Equals(const CVector2d& v, double tolerence) const;

  inline double Dot(const CVector2d& v) const;
  inline double Dot(const CVector2d* pV) const;

  // Misc
  double Normalize();
  double Normalize(double value);
  double Length() const;
  double LengthSquared() const;
  int IsCollinear(CVector2d *pVector) const;
  int IsCollinear(CVector2d &vector) const;
  inline void Negate();
};

//********************************************
// Dot
//********************************************
double
CVector2d::Dot(const CVector2d& v) const
{
  return (x() * v.x() +
	  y() * v.y());
}

//********************************************
// Dot
//********************************************
double
CVector2d::Dot(const CVector2d* pV) const
{
  return (x() * pV->x() +
	  y() * pV->y());
}


//********************************************
// Operator +=
//********************************************
CVector2d&
CVector2d::operator+=(const CVector2d& rVector)
{
  vec[0] += rVector.x();
  vec[1] += rVector.y();
  return *this;
}

//********************************************
// Operator -=
//********************************************
CVector2d&
CVector2d::operator-=(const CVector2d& rVector)
{
  vec[0] -= rVector.x();
  vec[1] -= rVector.y();
  return *this;
}

//********************************************
// Operator *=
//********************************************
CVector2d&
CVector2d::operator*=(double d)
{
  vec[0] *= d;
  vec[1] *= d;
  return *this;
}


#endif // _VECTOR2D_


