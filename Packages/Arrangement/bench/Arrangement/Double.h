#ifndef DOUBLE_H
#define DOUBLE_H

#include <math.h>
#include <iostream>

#define ZERO_EPSILON    0.000001
#define ERR_EPSILON     0.001

CGAL_BEGIN_NAMESPACE

class Double
{
 private:

  double           val;          // The value;

 public:

  // Constructors.
  Double () :
    val(0)
  {}

  Double (const double& value) :
    val(value)
  {}

  // Assignment from a double.
  const Double& operator= (const double& value)
  {
    val = value;
    return (*this);
  }

  // Arithmetic opertors.
  Double operator+ (const Double& x) const
  {
    return (Double(val + x.val));
  }

  Double operator- (const Double& x) const
  {
    return (Double(val - x.val));
  }

  Double operator* (const Double& x) const
  {
    return (Double(val * x.val));
  }

  Double operator/ (const Double& x) const
  {
    return (Double(val / x.val));
  }

  // Unary minus.
  Double operator- () const
  {
    return (Double(-val));
  }

  // Arithmetic opertors and assignment.
  void operator+= (const Double& x)
  {
    val += x.val;
  }

  void operator-= (const Double& x)
  {
    val -= x.val;
  }

  void operator*= (const Double& x)
  {
    val *= x.val;
  }

  void operator/= (const Double& x)
  {
    val /= x.val;
  }

  // Equality operators. Note that x equals y iff:
  //
  //     |x - y|
  //   ----------- < ERR_EPSILON 
  //    |x| + |y|
  //
  bool operator== (const Double& x) const
  {
    double numer = fabs(val - x.val);
    double denom = fabs(val) + fabs(x.val);

    if (denom < ZERO_EPSILON)
      return (true);           // The two numbers are very close to 0.
    else
      return (numer / denom < ERR_EPSILON);
  }

  bool operator!= (const Double& x) const
  {
    return (! (*this == x));
  }

  // Order operators.
  bool operator> (const Double& x) const
  {
    return (val > x.val && ! (*this == x));
  }

  bool operator>= (const Double& x) const
  {
    return (val > x.val || (*this == x));
  }

  bool operator< (const Double& x) const
  {
    return (val < x.val && ! (*this == x));
  }

  bool operator<= (const Double& x) const
  {
    return (val < x.val || (*this == x));
  }

  // Friend operators:
  friend Double operator+ (const double& x, const Double& y);
  friend Double operator- (const double& x, const Double& y);
  friend Double operator* (const double& x, const Double& y);
  friend Double operator/ (const double& x, const Double& y);

  // Friend functions:
  friend double to_double (const Double& x);
  friend bool is_finite (const Double& x);
  friend Double sqrt (const Double& x);
  friend Double pow (const Double& x, const Double& y);
  friend Double exp (const Double& x);
  friend Double log (const Double& x);
  friend Double sin (const Double& x);
  friend Double cos (const Double& x);
  friend Double tan (const Double& x);
  friend Double asin (const Double& x);
  friend Double acos (const Double& x);
  friend Double atan (const Double& x);
  friend Double atan2 (const Double& x, const Double& y);

  // I/O operations.
  friend std::istream& operator>> (std::istream& is, Double& x);
  friend std::ostream& operator<< (std::ostream& is, Double& x);

};

// Friend operators:
inline Double operator+ (const double& x, const Double& y)
{
  return (Double (x + y.val));
}

inline Double operator- (const double& x, const Double& y)
{
  return (Double (x - y.val));
}

inline Double operator* (const double& x, const Double& y)
{
  return (Double (x * y.val));
}

inline Double operator/ (const double& x, const Double& y)
{
  return (Double (x / y.val));
}

// Friend functions:
inline double to_double (const Double& x)
{
  return (x.val);
}

inline bool is_finite (const Double& x)
{
  return (::finite(x.val));
}

inline Double sqrt (const Double& x)
{
  return (Double (::sqrt(x.val)));
}

inline Double pow (const Double& x, const Double& y)
{
  return (Double (::pow(x.val, y.val)));
}

inline Double exp (const Double& x)
{
  return (Double (::exp(x.val)));
}

inline Double log (const Double& x)
{
  return (Double (::log(x.val)));
}

inline Double sin (const Double& x)
{
  return (Double (::sin(x.val)));
}

inline Double cos (const Double& x)
{
  return (Double (::cos(x.val)));
}

inline Double tan (const Double& x)
{
  return (Double (::tan(x.val)));
}

inline Double asin (const Double& x)
{
  return (Double (::asin(x.val)));
}

inline Double acos (const Double& x)
{
  return (Double (::acos(x.val)));
}

inline Double atan (const Double& x)
{
  return (Double (::atan(x.val)));
}

inline Double atan2 (const Double& x, const Double& y)
{
  return (Double (::atan2(x.val, y.val)));
}

// I/O operations.
inline std::istream& operator>> (std::istream& is, Double& x)
{
  is >> x.val;
  return (is);
}

inline std::ostream& operator<< (std::ostream& os, Double& x)
{
  os << x.val;
  return (os);
}

CGAL_END_NAMESPACE

#endif
