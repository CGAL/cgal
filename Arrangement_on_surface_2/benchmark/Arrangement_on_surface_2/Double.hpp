#ifndef DOUBLE_HPP
#define DOUBLE_HPP

#include <math.h>
#include <iostream>

#define ZERO_EPSILON    0.000001
#define ERR_EPSILON     0.001

namespace CGAL {

class Double {
public:
  double val;          // The value;

public:
  typedef Tag_false     Has_gcd;
  typedef Tag_true      Has_division;
  typedef Tag_true      Has_sqrt;

  typedef Tag_false     Has_exact_ring_operations;
  typedef Tag_false     Has_exact_division;
  typedef Tag_false     Has_exact_sqrt;

  typedef Tag_true      Is_real_embeddable;

  // Constructors.
  Double() : val(0) {}

  Double(const double & value) : val(value) {}

  // Assignment from a double.
  const Double & operator=(const double & value)
  {
    val = value;
    return(*this);
  }

  // Arithmetic operators.
  Double operator+(const Double & x) const { return Double(val + x.val); }

  Double operator-(const Double & x) const { return Double(val - x.val); }

  Double operator*(const Double & x) const { return Double(val * x.val); }

  Double operator/(const Double & x) const { return Double(val / x.val); }

  // Unary minus.
  Double operator-() const { return Double(-val); }

  // Arithmetic operators and assignment.
  void operator+=(const Double & x) { val += x.val; }

  void operator-=(const Double & x) { val -= x.val; }

  void operator*=(const Double & x) { val *= x.val; }

  void operator/=(const Double & x) { val /= x.val; }

  // Equality operators. Note that x equals y iff:
  //
  //     |x - y|
  //   ----------- <ERR_EPSILON
  //    |x| + |y|
  //
  bool operator==(const Double& x) const
  {
    double numer = fabs(val - x.val);
    double denom = fabs(val) + fabs(x.val);

    if (denom <ZERO_EPSILON)
      return true;           // The two numbers are very close to 0.
    else
      return (numer / denom <ERR_EPSILON);
  }

  bool operator!=(const Double& x) const
  {
    return !(*this == x);
  }

  // Order operators.
  bool operator>(const Double& x) const
  {
    return (val > x.val && !(*this == x));
  }

  bool operator>=(const Double& x) const
  {
    return (val > x.val ||(*this == x));
  }

  bool operator<(const Double & x) const
  {
    return (val <x.val && !(*this == x));
  }

  bool operator<=(const Double & x) const
  {
    return (val <x.val || (*this == x));
  }

  // Friend operators:
  friend Double operator+(const double & x, const Double & y);
  friend Double operator-(const double & x, const Double & y);
  friend Double operator*(const double & x, const Double & y);
  friend Double operator/(const double & x, const Double & y);

  // Friend functions:
  friend double to_double(const Double & x);
  friend std::pair<double,double> to_interval(const Double & x);
  friend bool is_finite(const Double & x);
  friend Double sqrt(const Double & x);
  friend Double pow(const Double & x, const Double & y);
  friend Double exp(const Double & x);
  friend Double log(const Double & x);
  friend Double sin(const Double& x);
  friend Double cos(const Double& x);
  friend Double tan(const Double& x);
  friend Double asin(const Double& x);
  friend Double acos(const Double& x);
  friend Double atan(const Double& x);
  friend Double atan2(const Double& x, const Double& y);
  friend Double fabs(const Double& x);

  // I/O operations.
  friend std::istream & operator>>(std::istream & is, Double & x);
  friend std::ostream & operator<<(std::ostream & os, Double & x);
};

inline io_Operator io_tag(const Double &) { return io_Operator(); }

// Friend operators:
inline Double operator+(const double & x, const Double & y)
{
  return (Double(x + y.val));
}

inline Double operator-(const double & x, const Double & y)
{
  return (Double(x - y.val));
}

inline Double operator*(const double & x, const Double & y)
{
  return (Double(x * y.val));
}

inline Double operator/(const double & x, const Double & y)
{
  return (Double(x / y.val));
}

// Order operators.
inline bool operator<(const double & a, const Double & b) { return a <b; }

inline bool operator>(const double & a, const Double & b) { return b <a; }

inline bool operator>=(const double & a, const Double & b) { return !(a <b); }

inline bool operator<=(const double & a, const Double & b) { return !(a > b); }

inline bool operator==(const double & a, const Double & b) { return (a == b); }

inline bool operator!=(const double & a, const Double & b) { return !(a == b); }

// Friend functions:
inline double to_double(const Double & x) { return x.val; }

inline std::pair<double,double> to_interval(const Double & x)
{
  return std::make_pair(x.val - ERR_EPSILON, x.val + ERR_EPSILON);
}

inline bool is_finite(const Double & x) { return ::finite(x.val); }

inline Double sqrt(const Double & x) { return Double(::sqrt(x.val)); }

inline Double pow(const Double & x, const Double & y)
{
  return Double(::pow(x.val, y.val));
}

inline Double exp(const Double & x) { return Double(::exp(x.val)); }

inline Double log(const Double & x) { return Double(::log(x.val)); }

inline Double sin(const Double & x) { return Double(::sin(x.val)); }

inline Double cos(const Double & x) { return Double(::cos(x.val)); }

inline Double tan(const Double & x) { return Double(::tan(x.val)); }

inline Double asin(const Double & x) { return Double(::asin(x.val)); }

inline Double acos(const Double & x) { return Double(::acos(x.val)); }

inline Double atan(const Double & x) { return Double(::atan(x.val)); }

inline Double atan2(const Double & x, const Double & y)
{ return Double(::atan2(x.val, y.val)); }

inline Double fabs(const Double & x) { return Double(::fabs(x.val)); }

// I/O operations.
inline std::istream & operator>>(std::istream & is, Double & x)
{
  is >> x.val;
  return is;
}

inline std::ostream & operator<<(std::ostream & os, const Double & x)
{
  os <<x.val;
  return os;
}

// Real embeddable traits
template <> class Real_embeddable_traits<Double>
  : public Real_embeddable_traits_base<Double> {
public:

// GCC is faster with std::fabs().
#ifdef __GNUG__
  class Abs : public Unary_function<Type, Type> {
  public:
    Type operator()( const Type& x ) const {
      return fabs( x );
    }
  };
#endif

  typedef INTERN_RET::To_double_by_conversion<Type>     To_double;
  typedef INTERN_RET::To_interval_by_conversion<Type>   To_interval;

// Is_finite depends on platform
#ifdef CGAL_CFG_IEEE_754_BUG
  class Is_finite : public Unary_function<Type, bool > {
  public:
    bool operator()( const Type& x ) const {
      Type d = x;
      IEEE_754_double* p = reinterpret_cast<IEEE_754_double*>(&d);
      return is_finite_by_mask_double( p->c.H );
    }
  };
#elif defined CGAL_CFG_NUMERIC_LIMITS_BUG
  class Is_finite : public Unary_function<Type, bool > {
  public:
    bool operator()( const Type& x ) const {
      return (x == x) && (is_valid(x-x));
    }
  };
#else
  class Is_finite : public Unary_function<Type, bool > {
  public:
    bool operator()( const Type& x ) const {
      return (x != std::numeric_limits<Type>::infinity())
        && (-x != std::numeric_limits<Type>::infinity())
        && is_valid(x);
    }
  };
#endif
};

} //namespace CGAL

#endif
