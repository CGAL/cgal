#ifndef CGAL_SQUARE_ROOT_2_H
#define CGAL_SQUARE_ROOT_2_H

#include <CGAL/predicates/Square_root_1.h>




CGAL_BEGIN_NAMESPACE


template<class NT>
class Square_root_2
{
private:
  typedef Square_root_2<NT>  Self;
  typedef Square_root_1<NT>  Sqrt_1;

  NT a0_, a1_, a2_, a3_;
  NT A_, B_;

public:
  typedef NT                 FT;
  typedef NT                 RT;
  
public:
  Square_root_2()
    : a0_(0), a1_(0), a2_(0), a3_(0), A_(0), B_(0) {}
  Square_root_2(int i)
    : a0_(i), a1_(0), a2_(0), a3_(0), A_(0), B_(0) {}
  Square_root_2(const NT& a)
    : a0_(a), a1_(0), a2_(0), a3_(0), A_(0), B_(0) {}
  Square_root_2(const NT& a0, const NT& a1, const NT& a2,
		const NT& a3, const NT& A, const NT& B)
    : a0_(a0), a1_(a1), a2_(a2), a3_(a3), A_(A), B_(B)
  {
    //    CGAL_assertion( !(CGAL_NTS is_negative(A_)) );
    //    CGAL_assertion( !(CGAL_NTS is_negative(B_)) );
  }

  Square_root_2(const Square_root_2<NT>& other)
    : a0_(other.a0_), a1_(other.a1_), a2_(other.a2_),
      a3_(other.a3_), A_(other.A_), B_(other.B_) {}


  NT a() const { return a0_; }
  NT b() const { return a1_; }
  NT c() const { return a2_; }
  NT d() const { return a3_; }
  NT e() const { return A_; }
  NT f() const { return B_; }


  Self operator*(const Self& b) const
  {
#if CHECK_CGAL_PRECONDITIONS
    CGAL_precondition( CGAL_NTS compare(A_, b.A_) == EQUAL );
    CGAL_precondition( CGAL_NTS compare(B_, b.B_) == EQUAL );
#endif

    NT a0 = a0_ * b.a0_ + a1_ * b.a1_ * A_ + a2_ * b.a2_ * B_
      + a3_ * b.a3_ * A_ * B_;
    NT a1 = a0_ * b.a1_ + a1_ * b.a0_ + (a2_ * b.a3_ + a3_ * b.a2_) * B_;
    NT a2 = a0_ * b.a2_ + a2_ * b.a0_ + (a1_ * b.a3_ + a3_ * b.a1_) * A_;
    NT a3 = a0_ * b.a3_ + a3_ * b.a0_ + a1_ * b.a2_ + a2_ * b.a1_;

    return Self(a0, a1, a2, a3, A_, B_);
  }

  Self operator-() const
  {
    return Self(-a0_, -a1_, -a2_, -a3_, A_, B_);
  }

  Self operator+() const
  {
    return (*this);
  }

  Self square() const
  {
    NT a0 = CGAL_NTS square(a0_) + CGAL_NTS square(a1_) * A_
      + CGAL_NTS square(a2_) * B_ + CGAL_NTS square(a3_) * A_ * B_;
    NT a1_half = a0_ * a1_ + a2_ * a3_ * B_;
    NT a2_half = a0_ * a2_ + a1_ * a3_ * A_;
    NT a3_half = a0_ * a3_ + a1_ * a2_;

    NT a1 = a1_half + a1_half;
    NT a2 = a2_half + a2_half;
    NT a3 = a3_half + a3_half;

    return Self(a0, a1, a2, a3, A_, B_);
  }

  Sign sign() const
  {
    //    std::cout << "cp sgn 1" << std::flush;

    Sqrt_1 x(a0_, a1_, A_);
    Sqrt_1 y(a2_, a3_, A_);

#if 0
    std::cout << "a0: " << a0_ << std::endl;
    std::cout << "a1: " << a1_ << std::endl;
    std::cout << "A : " << A_ << std::endl;
    std::cout << "x: " << x << std::endl;

    std::cout << " 2" << std::flush;
#endif

    Sign s_x = CGAL_NTS sign(x);

    //    std::cout << " 3" << std::flush;

    Sign s_y = CGAL_NTS sign(y);

    //    std::cout << " 4" << std::flush;

    Sign s_B = CGAL_NTS sign(B_);

    //    std::cout << " 5" << std::flush;
    //    std::cout << std::endl;

    if ( s_B == ZERO ) {
      return s_x;
    } else if ( s_x == s_y ) {
      return s_x;
    } else if ( s_x == ZERO ) {
      return s_y;
    } else if ( s_y == ZERO ) {
      return s_x;
    } else {
      Sqrt_1 Q = CGAL_NTS square(x) - CGAL_NTS square(y) * B_;
      return Sign(s_x * CGAL_NTS sign(Q));
    }
  }

  double to_double() const
  {
    // THIS MUST BE CHECK WITH SYLVAIN FOR CORRECTNESS
    double a0d = CGAL_NTS to_double(a0_);
    double a1d = CGAL_NTS to_double(a1_);
    double a2d = CGAL_NTS to_double(a2_);
    double a3d = CGAL_NTS to_double(a3_);
    double Ad = CGAL_NTS to_double(A_);
    double Bd = CGAL_NTS to_double(B_);

    return (a0d + a1d * CGAL_NTS sqrt(Ad) + a2d * CGAL_NTS sqrt(Bd)
	    + a3d * CGAL_NTS sqrt(Ad * Bd));
  }
};





// operator *
template<class NT>
inline
Square_root_2<NT>
operator*(const Square_root_2<NT>& x, const NT& n)
{
  return Square_root_2<NT>(x.a() * n, x.b() * n, x.c() * n,
			   x.d() * n, x.e(), x.f());
}


template<class NT>
inline
Square_root_2<NT>
operator*(const NT& n, const Square_root_2<NT>& x)
{
  return (x * n);
}



// operator +
template<class NT>
inline
Square_root_2<NT>
operator+(const Square_root_2<NT>& x, const NT& n)
{
  return Square_root_2<NT>(x.a() + n, x.b(), x.c(), x.d(),
			   x.e(), x.f());
}


template<class NT>
inline
Square_root_2<NT>
operator+(const NT& n, const Square_root_2<NT>& x)
{
  return (x + n);
}

template<class NT>
inline
Square_root_2<NT>
operator+(const Square_root_2<NT>& x, const Square_root_2<NT>& y)
{
#if CHECK_CGAL_PRECONDITIONS
  CGAL_precondition( CGAL_NTS compare(x.e(), y.e()) == EQUAL );
  CGAL_precondition( CGAL_NTS compare(x.f(), y.f()) == EQUAL );
#endif

  return Square_root_2<NT>(x.a() + y.a(), x.b() + y.b(),
			   x.c() + y.c(), x.d() + y.d(),
			   x.e(), x.f());
}



// operator -
template<class NT>
inline
Square_root_2<NT>
operator-(const Square_root_2<NT>& x, const NT& n)
{
  return x + (-n);
}


template<class NT>
inline
Square_root_2<NT>
operator-(const NT& n, const Square_root_2<NT>& x)
{
  return -(x - n);
}

template<class NT>
inline
Square_root_2<NT>
operator-(const Square_root_2<NT>& x, const Square_root_2<NT>& y)
{
  return (x + (-y));
}





namespace NTS {


  template<class NT>
  inline
  bool
  is_positive(const Square_root_2<NT>& x)
  {
    return sign(x) == POSITIVE;
  }

  template<class NT>
  inline
  bool
  is_negative(const Square_root_2<NT>& x)
  {
    return sign(x) == NEGATIVE;
  }

  template<class NT>
  inline
  bool
  is_zero(const Square_root_2<NT>& x)
  {
    return sign(x) == ZERO;
  }


  template<class NT>
  inline
  Sign
  sign(const Square_root_2<NT>& x)
  {
    return x.sign();
  }

  template<class NT>
  inline
  Square_root_2<NT>
  square(const Square_root_2<NT>& x)
  {
    return x.square();
  }

  template<class NT>
  inline
  Comparison_result
  compare(const Square_root_2<NT>& x,
	  const Square_root_2<NT>& y)
  {
#if CHECK_CGAL_PRECONDITIONS
    CGAL_precondition( CGAL_NTS compare(x.e(), y.e()) == EQUAL );
    CGAL_precondition( CGAL_NTS compare(x.f(), y.f()) == EQUAL );
#endif

    Sign s = CGAL_NTS sign(x - y);

    if ( s == ZERO ) { return EQUAL; }
    return (s == POSITIVE) ? LARGER : SMALLER;
  }

} // namespace NTS



// operator <<
template<class Stream, class NT>
inline
Stream&
operator<<(Stream& os, const Square_root_2<NT>& x)
{
  os << "(" << x.a()  << ")+(" << x.b() << ") sqrt{" << x.e() << "}";
  os << "+(" << x.c() << ") sqrt{" << x.f() << "}";
  os << "+(" << x.d() << ") sqrt{(" << x.e() << ") (" << x.f() 
     << ")}";
  return os;

#if 0
  if ( CGAL_NTS is_zero(x) ) {
    os << "0";
    return os;
  }

  Square_root_2<NT> One(NT(1), NT(0), NT(0), NT(0), x.e(), x.f());

  if ( CGAL_NTS sign(x.a()) != ZERO ) {
    os << x.a();
    if ( CGAL_NTS is_positive(x.b()) ) {
      os << "+";
    }
  }

  if ( CGAL_NTS sign(x.b()) != ZERO &&
       CGAL_NTS sign(x.e()) != ZERO ) {
    if ( CGAL_NTS sign(x.b() - One) != ZERO ) {
      os << x.b() << " ";
    }
    os << "sqrt{" << x.e() << "}";
    if ( CGAL_NTS is_positive(x.c()) ) {
      os << "+";
    }
  }

  if ( CGAL_NTS sign(x.c()) != ZERO &&
       CGAL_NTS sign(x.f()) != ZERO ) {
    if ( CGAL_NTS sign(x.c() - One) != ZERO ) {
      os << x.c() << " ";
    }
    os << "sqrt{" << x.f() << "}";
    if ( CGAL_NTS is_positive(x.d()) ) {
      os << "+";
    }
  }

  if ( CGAL_NTS sign(x.d()) != ZERO &&
       CGAL_NTS sign(x.e()) != ZERO &&
       CGAL_NTS sign(x.f()) != ZERO ) {
    if ( CGAL_NTS sign(x.d() - One) != ZERO ) {
      os << x.d() << " ";
    }
    os << "sqrt{" << x.e() << " " << x.f() << "}";
  }

  return os;
#endif
}






CGAL_END_NAMESPACE


#endif // CGAL_SQUARE_ROOT_2_H
