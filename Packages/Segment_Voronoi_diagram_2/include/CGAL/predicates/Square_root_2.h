// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>




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
#if CHECK_CGAL_PRECONDITIONS
    CGAL_precondition( !(CGAL::is_negative(A_)) );
    CGAL_precondition( !(CGAL::is_negative(B_)) );
#endif
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
    CGAL_precondition( CGAL::compare(A_, b.A_) == EQUAL );
    CGAL_precondition( CGAL::compare(B_, b.B_) == EQUAL );
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
    NT a0 = CGAL::square(a0_) + CGAL::square(a1_) * A_
      + CGAL::square(a2_) * B_ + CGAL::square(a3_) * A_ * B_;
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
    Sqrt_1 x(a0_, a1_, A_);
    Sqrt_1 y(a2_, a3_, A_);

    Sign s_x = Number_type_traits<Sqrt_1>::sign(x);
    Sign s_y = Number_type_traits<Sqrt_1>::sign(y);
    Sign s_B = Number_type_traits<Sqrt_1>::sign(B_);

    if ( s_B == ZERO ) {
      return s_x;
    } else if ( s_x == s_y ) {
      return s_x;
    } else if ( s_x == ZERO ) {
      return s_y;
    } else if ( s_y == ZERO ) {
      return s_x;
    } else {
      Sqrt_1 Q = CGAL::square(x) - CGAL::square(y) * B_;
      return Sign(s_x * Number_type_traits<Sqrt_1>::sign(Q));
    }
  }

  double to_double() const
  {
    // THIS MUST BE CHECK WITH SYLVAIN FOR CORRECTNESS
    double a0d = CGAL::to_double(a0_);
    double a1d = CGAL::to_double(a1_);
    double a2d = CGAL::to_double(a2_);
    double a3d = CGAL::to_double(a3_);
    double Ad = CGAL::to_double(A_);
    double Bd = CGAL::to_double(B_);

    return (a0d + a1d * CGAL::sqrt(Ad) + a2d * CGAL::sqrt(Bd)
	    + a3d * CGAL::sqrt(Ad * Bd));
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
  CGAL_precondition( CGAL::compare(x.e(), y.e()) == EQUAL );
  CGAL_precondition( CGAL::compare(x.f(), y.f()) == EQUAL );
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



//===================================================================


template<class NT>
struct Number_type_traits< Square_root_2<NT> >
{
  static inline bool is_positive(const Square_root_2<NT>& x) {
    return x.sign() == POSITIVE;
  }

  static inline bool is_negative(const Square_root_2<NT>& x) {
    return x.sign() == NEGATIVE;
  }

  static inline bool is_zero(const Square_root_2<NT>& x) {
    return x.sign() == ZERO;
  }

  static inline Sign sign(const Square_root_2<NT>& x) {
    return x.sign();
  }

  static inline Square_root_2<NT> square(const Square_root_2<NT>& x) {
    return x.square();
  }

  static inline
  Comparison_result compare(const Square_root_2<NT>& x,
			    const Square_root_2<NT>& y)
  {
#if CHECK_CGAL_PRECONDITIONS
    CGAL_precondition( CGAL::compare(x.e(), y.e()) == EQUAL );
    CGAL_precondition( CGAL::compare(x.f(), y.f()) == EQUAL );
#endif

    //  Sign s = CGAL::sign(x - y);
    Sign s = (x - y).sign();

    if ( s == ZERO ) { return EQUAL; }
    return (s == POSITIVE) ? LARGER : SMALLER;
  }
};

template<class NT>
inline
bool
is_positive(const Square_root_2<NT>& x)
{
  return Number_type_traits< Square_root_2<NT> >::is_positive(x);
}

template<class NT>
inline
bool
is_negative(const Square_root_2<NT>& x)
{
  return Number_type_traits< Square_root_2<NT> >::is_negative(x);
}

template<class NT>
inline
bool
is_zero(const Square_root_2<NT>& x)
{
  return Number_type_traits< Square_root_2<NT> >::is_zero(x);
}


template<class NT>
inline
Sign
sign(const Square_root_2<NT>& x)
{
  return Number_type_traits< Square_root_2<NT> >::sign(x);
}

template<class NT>
inline
Square_root_2<NT>
square(const Square_root_2<NT>& x)
{
  return Number_type_traits< Square_root_2<NT> >::square(x);
}

template<class NT>
inline
Comparison_result
compare(const Square_root_2<NT>& x,
	const Square_root_2<NT>& y)
{
  return Number_type_traits< Square_root_2<NT> >::compare(x, y);
}

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
}






CGAL_END_NAMESPACE


#endif // CGAL_SQUARE_ROOT_2_H
