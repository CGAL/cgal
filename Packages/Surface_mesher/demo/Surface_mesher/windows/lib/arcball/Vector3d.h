///////////////////////////////////////////////////////////////////////////
//                                                                       //
//  Class: CVector3d                                                     //
//                                                                       //
//  3D Vector to represent points or directions.  Each component of the  //
//  vector is stored as a doubleing point number.                        //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#ifndef _VECTOR_3D_
#define _VECTOR_3D_

class CVector3d
{

private :

protected:

 double vec[3];  // Storage for the vector components

public :

  // Constructors
  CVector3d() { }
  CVector3d(const double x, const double y, const double z)
    { vec[0] = x; vec[1] = y; vec[2] = z; }
   CVector3d(const double v[3])
    { vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; }
  CVector3d(const CVector3d &vector)  { Set(vector); }
  CVector3d(const CVector3d *pVector) { Set(pVector); }
  CVector3d(const CVector3d &a, const CVector3d& b)
    { Set(b - a); }
  CVector3d(const CVector3d *a, const CVector3d *b)
    { Set(*b - *a); }
  
  virtual ~CVector3d() { }

  // Debug
  void Trace() const;

  // Data setting
  void Clear()
    { vec[0] = 0.f; vec[1] = 0.f; vec[2] = 0.f; }
  void Set(const CVector3d *pVector) { Set(pVector->GetArray()); }
  void Set(const CVector3d& vector)  { Set(vector.GetArray()); }
  void Set(const double x, const double y, const double z)
    { vec[0] = x; vec[1] = y; vec[2] = z; }
  void Set(const double v[3])
    { vec[0] = v[0]; vec[1] = v[1]; vec[2] = v[2]; }
  void Set(const CVector3d& a, const CVector3d& b)
    { Set(b - a); }
  void Set(const CVector3d *a, const CVector3d *b)
    { Set(*b - *a); }

  // Data Access
  const double* GetArray() const { return vec; }
  void         Get(double& x, double& y, double& z) const;

  // Per coordinate (explicit inline functions)
  void x(double newX) { vec[0] = newX; }
  void y(double newY) { vec[1] = newY; }
  void z(double newZ) { vec[2] = newZ; }

  // Data access (explicit inline functions)
  double x() const { return (vec[0]); }
  double y() const { return (vec[1]); }
  double z() const { return (vec[2]); }

  // Data access using indices
  double&       operator[](int i)       { return (vec[i]); }
  const double& operator[](int i) const { return (vec[i]); }

  // Operators
  CVector3d& operator+=(const CVector3d& rVector);
  CVector3d& operator+=(const CVector3d* pVector);
  CVector3d& operator-=(const CVector3d& rVector);
  CVector3d& operator-=(const CVector3d* pVector);
  CVector3d& operator*=(double d);
  CVector3d& operator/=(double d)
    { return *this *= (1.f/d); }

  // Nondestructive unary negation - returns a new vector
  CVector3d  operator -() const;

  // Binary operators
  friend CVector3d operator+(const CVector3d& u, const CVector3d& v);
  friend CVector3d operator-(const CVector3d& u, const CVector3d& v);
  friend CVector3d operator*(double s,            const CVector3d& u);
  friend CVector3d operator*(const CVector3d& u, double s)
    { return s * u; }
  friend CVector3d operator/(const CVector3d& u, double s)
    { return (1.f/s) * u; }
  friend CVector3d operator^(const CVector3d& u, const CVector3d& v);
  friend int       operator==(const CVector3d& v1, const CVector3d& v2);
  friend int       operator!=(const CVector3d& v1, const CVector3d& v2)
    { return !(v1 == v2); }

  int Equals(const CVector3d& v, double tolerence) const;

  inline double     Dot(const CVector3d& v) const;
  inline double     Dot(const CVector3d* pV) const;
  CVector3d        Cross(const CVector3d& v) const;
  CVector3d        Cross(const CVector3d* pV) const;

  // Misc
  double Normalize();
  double Normalize(double value);
  double Length() const;
  double LengthSquared() const;
  int IsCollinear(CVector3d *pVector) const;
  int IsCollinear(CVector3d &vector) const;
	void Negate();
	CVector3d Rotate(double angle,CVector3d Around); 
	CVector3d Projection(const CVector3d* pV) const;
};


//********************************************
// Dot
//********************************************
double
CVector3d::Dot(const CVector3d& v) const
{
  return (x() * v.x() +
	  y() * v.y() +
	  z() * v.z());
}

//********************************************
// Dot
//********************************************
double
CVector3d::Dot(const CVector3d* pV) const
{
  return (x() * pV->x() +
	  y() * pV->y() +
	  z() * pV->z());
}

#endif // _VECTOR_3D_
