#ifndef CGAL_SQUARE_ROOT_1_H
#define CGAL_SQUARE_ROOT_1_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>

#include <iostream>

#define CHECK_CGAL_PRECONDITIONS 0

CGAL_BEGIN_NAMESPACE


template<class NT>
class Square_root_1
{
private:
  NT x, y, r;

private:
  typedef Square_root_1<NT>   Self;

public:
  typedef NT                  FT;
  typedef NT                  RT;

  Square_root_1() : x(0), y(0), r(0) {}
  Square_root_1(int i) : x(i), y(0), r(0) {}
  Square_root_1(const NT& a) : x(a), y(0), r(0) {}
  Square_root_1(const NT& a, const NT& b, const NT& c)
    : x(a), y(b), r(c)
  {
    //    CGAL_assertion( !(CGAL_NTS is_negative(r)) );
  }

  Square_root_1(const Square_root_1<NT>& other)
    : x(other.x), y(other.y), r(other.r) {}


  NT a() const { return x; }
  NT b() const { return y; }
  NT c() const { return r; }

  Self square() const
  {
    NT xy = x * y;

    return Self(CGAL_NTS square(x) + CGAL_NTS square(y) * r,
		xy + xy,
		r);
  }

  Sign sign() const
  {
    Sign sx = CGAL_NTS sign(x);

    if ( CGAL_NTS sign(r) == ZERO )  { return sx; }

    Sign sy = CGAL_NTS sign(y);

    if ( sx == sy )  { return sx; }
    if ( sx == ZERO )  { return sy; }

    return Sign( sx * CGAL_NTS compare( CGAL_NTS square(x),
					r * CGAL_NTS square(y) )
		 );
  }

  Self operator+() const
  {
    return (*this);
  }

  Self operator-() const
  {
    return Self(-x, -y, r);
  }

  double to_double() const
  {
    // THIS MUST BE CHECK WITH SYLVAIN FOR CORRECTNESS
    double xd = CGAL_NTS to_double(x);
    double yd = CGAL_NTS to_double(y);
    double rd = CGAL_NTS to_double(r);

    return (xd + yd * CGAL_NTS sqrt(rd));
  }

  std::pair<double,double> to_interval() const
  {
    // THIS MUST BE CHECK WITH SYLVAIN FOR CORRECTNESS
    std::pair<double,double> x_ivl = CGAL::to_interval(x);
    std::pair<double,double> y_ivl = CGAL::to_interval(y);
    std::pair<double,double> r_ivl = CGAL::to_interval(r);

    std::pair<double,double> sqrt_r_ivl(CGAL_NTS sqrt(r_ivl.first),
					CGAL_NTS sqrt(r_ivl.second));

    std::pair<double,double>
      ivl(x_ivl.first + y_ivl.first * sqrt_r_ivl.first,
	  x_ivl.second + y_ivl.second * sqrt_r_ivl.second);

    return ivl;
  }
};


// operator *
template<class NT>
inline
Square_root_1<NT>
operator*(const Square_root_1<NT>& x, const NT& n)
{
  return Square_root_1<NT>(x.a() * n, x.b() * n, x.c());
}


template<class NT>
inline
Square_root_1<NT>
operator*(const NT& n, const Square_root_1<NT>& x)
{
  return (x * n);
}

template<class NT>
inline
Square_root_1<NT>
operator*(const Square_root_1<NT>& x, const Square_root_1<NT>& y)
{
#if CHECK_CGAL_PRECONDITIONS
  CGAL_precondition( CGAL_NTS compare(x.c(), y.c()) == EQUAL );
#endif

  NT a = x.a() * y.a() + x.b() * y.b() * x.c();
  NT b = x.a() * y.b() + x.b() * y.a();

  return Square_root_1<NT>(a, b, x.c());
}


// operator +
template<class NT>
inline
Square_root_1<NT>
operator+(const Square_root_1<NT>& x, const NT& n)
{
  return Square_root_1<NT>(x.a() + n, x.b(), x.c());
}


template<class NT>
inline
Square_root_1<NT>
operator+(const NT& n, const Square_root_1<NT>& x)
{
  return (x + n);
}

template<class NT>
inline
Square_root_1<NT>
operator+(const Square_root_1<NT>& x, const Square_root_1<NT>& y)
{
#if CHECK_CGAL_PRECONDITIONS
  CGAL_precondition( CGAL_NTS compare(x.c(), y.c()) == EQUAL );
#endif

  return Square_root_1<NT>(x.a() + y.a(), x.b() + y.b(), x.c());
}



// operator -
template<class NT>
inline
Square_root_1<NT>
operator-(const Square_root_1<NT>& x, const NT& n)
{
  return x + (-n);
}


template<class NT>
inline
Square_root_1<NT>
operator-(const NT& n, const Square_root_1<NT>& x)
{
  return -(x - n);
}

template<class NT>
inline
Square_root_1<NT>
operator-(const Square_root_1<NT>& x, const Square_root_1<NT>& y)
{
  return (x + (-y));
}





template<class NT>
inline
std::pair<double,double>
to_interval(const Square_root_1<NT>& x)
{
  return x.to_interval();
}


namespace NTS {


  template<class NT>
  inline
  bool
  is_positive(const Square_root_1<NT>& x)
  {
    return sign(x) == POSITIVE;
  }

  template<class NT>
  inline
  bool
  is_negative(const Square_root_1<NT>& x)
  {
    return sign(x) == NEGATIVE;
  }

  template<class NT>
  inline
  bool
  is_zero(const Square_root_1<NT>& x)
  {
    return sign(x) == ZERO;
  }


  template<class NT>
  inline
  Sign
  sign(const Square_root_1<NT>& x)
  {
    return x.sign();
  }

  template<class NT>
  inline
  Square_root_1<NT>
  square(const Square_root_1<NT>& x)
  {
    return x.square();
  }

  template<class NT>
  inline
  Comparison_result
  compare(const Square_root_1<NT>& x, const Square_root_1<NT>& y)
  {
#if CHECK_CGAL_PRECONDITIONS
    CGAL_precondition( CGAL_NTS compare(x.c(), y.c()) == EQUAL );
#endif

    Sign s = CGAL_NTS sign(x - y);

    if ( s == ZERO ) { return EQUAL; }
    return (s == POSITIVE) ? LARGER : SMALLER;
  }

  template<class NT>
  inline
  double
  to_double(const Square_root_1<NT>& x)
  {
    return x.to_double();
  }

} // namespace NTS



// operator <<
template<class Stream, class NT>
inline
Stream&
operator<<(Stream& os, const Square_root_1<NT>& x)
{
  os << "(" << x.a()  << ")+(" << x.b() << ") sqrt{" << x.c() << "}";
  return os;

#if 0
  if ( CGAL_NTS is_zero(x) ) {
    os << "0";
    return os;
  }

  //  Square_root_1<NT> One(NT(1), NT(0), x.r());

  if ( CGAL_NTS sign(x.a()) != ZERO ) {
    os << x.a();
    if ( CGAL_NTS is_positive(x.b()) ) {
      os << "+";
    }
  }

  if ( CGAL_NTS sign(x.b()) != ZERO &&
       CGAL_NTS sign(x.c()) != ZERO ) {
    //    if ( CGAL_NTS sign(x.b() - One) != ZERO ) {
      os << x.b() << " ";
      //    }

    os << "sqrt{" << x.c() << "}";
  }

  return os;
#endif
}



CGAL_END_NAMESPACE



#endif // CGAL_SQUARE_ROOT_1_H
