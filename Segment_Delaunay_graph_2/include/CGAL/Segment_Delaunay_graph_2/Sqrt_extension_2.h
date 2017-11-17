// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>




#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_SQRT_EXTENSION_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_SQRT_EXTENSION_2_H

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Sqrt_extension.h>




namespace CGAL {


template<class NT>
class Sqrt_extension_2
{
private:
  typedef Sqrt_extension_2<NT>  Self;
  typedef Sqrt_extension<NT,NT,Tag_true>  Sqrt_1;

  NT a0_, a1_, a2_, a3_;
  NT A_, B_;

public:
  typedef NT                 FT;
  typedef NT                 RT;
  
public:
  Sqrt_extension_2()
    : a0_(0), a1_(0), a2_(0), a3_(0), A_(0), B_(0) {}
  Sqrt_extension_2(int i)
    : a0_(i), a1_(0), a2_(0), a3_(0), A_(0), B_(0) {}
  Sqrt_extension_2(const NT& a)
    : a0_(a), a1_(0), a2_(0), a3_(0), A_(0), B_(0) {}
  Sqrt_extension_2(const NT& a0, const NT& a1, const NT& a2,
		   const NT& a3, const NT& A, const NT& B)
    : a0_(a0), a1_(a1), a2_(a2), a3_(a3), A_(A), B_(B)
  {
    CGAL_exactness_precondition( !(CGAL::is_negative(A_)) );
    CGAL_exactness_precondition( !(CGAL::is_negative(B_)) );
  }

  Sqrt_extension_2(const Sqrt_extension_2<NT>& other)
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
    CGAL_exactness_precondition( CGAL::compare(A_, b.A_) == EQUAL );
    CGAL_exactness_precondition( CGAL::compare(B_, b.B_) == EQUAL );

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

    Sign s_x = CGAL_NTS sign(x);
    Sign s_y = CGAL_NTS sign(y);
    Sign s_B = CGAL_NTS sign(B_);

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
      return s_x * CGAL_NTS sign(Q);
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
Sqrt_extension_2<NT>
operator*(const Sqrt_extension_2<NT>& x, const NT& n)
{
  return Sqrt_extension_2<NT>(x.a() * n, x.b() * n, x.c() * n,
			      x.d() * n, x.e(), x.f());
}


template<class NT>
inline
Sqrt_extension_2<NT>
operator*(const NT& n, const Sqrt_extension_2<NT>& x)
{
  return (x * n);
}



// operator +
template<class NT>
inline
Sqrt_extension_2<NT>
operator+(const Sqrt_extension_2<NT>& x, const NT& n)
{
  return Sqrt_extension_2<NT>(x.a() + n, x.b(), x.c(), x.d(), x.e(), x.f());
}


template<class NT>
inline
Sqrt_extension_2<NT>
operator+(const NT& n, const Sqrt_extension_2<NT>& x)
{
  return (x + n);
}

template<class NT>
inline
Sqrt_extension_2<NT>
operator+(const Sqrt_extension_2<NT>& x, const Sqrt_extension_2<NT>& y)
{
  CGAL_exactness_precondition( CGAL::compare(x.e(), y.e()) == EQUAL );
  CGAL_exactness_precondition( CGAL::compare(x.f(), y.f()) == EQUAL );

  return Sqrt_extension_2<NT>(x.a() + y.a(), x.b() + y.b(),
			      x.c() + y.c(), x.d() + y.d(),
			      x.e(), x.f());
}



// operator -
template<class NT>
inline
Sqrt_extension_2<NT>
operator-(const Sqrt_extension_2<NT>& x, const NT& n)
{
  return x + (-n);
}


template<class NT>
inline
Sqrt_extension_2<NT>
operator-(const NT& n, const Sqrt_extension_2<NT>& x)
{
  return -(x - n);
}

template<class NT>
inline
Sqrt_extension_2<NT>
operator-(const Sqrt_extension_2<NT>& x, const Sqrt_extension_2<NT>& y)
{
  return (x + (-y));
}



//===================================================================


template <class NT> 
class Algebraic_structure_traits<Sqrt_extension_2<NT> >
    :public Algebraic_structure_traits_base<Sqrt_extension_2<NT>,CGAL::Integral_domain_without_division_tag>{
private:
    typedef Algebraic_structure_traits<NT> AST_NT;
public:
    typedef Sqrt_extension_2<NT> Algebraic_structure;
    typedef typename AST_NT::Is_exact Is_exact;
};

template<class NT>
class Real_embeddable_traits<Sqrt_extension_2<NT> >{
private:
    typedef Real_embeddable_traits<NT> RET_NT;
public:
    
    typedef Sqrt_extension_2<NT> Real_embeddable;
    
    class Abs 
        : public CGAL::unary_function< Real_embeddable, Real_embeddable >{
    public:
        Real_embeddable operator()(const Real_embeddable& x) const {
            return (x>=0)?x:-x;
        }
    };    

    class Sgn 
        : public CGAL::unary_function< Real_embeddable, CGAL::Sign >{
    public:
        CGAL::Sign operator()(const Real_embeddable& x) const {
            return x.sign();
        }
    };
    
    class Compare 
        : public CGAL::binary_function< Real_embeddable,
                                  Real_embeddable, 
                                  CGAL::Comparison_result >{
    public:
        CGAL::Comparison_result operator()(
                const Real_embeddable& x, 
                const Real_embeddable& y) const {
            CGAL_exactness_precondition( CGAL::compare(x.e(), y.e()) == EQUAL );
            CGAL_exactness_precondition( CGAL::compare(x.f(), y.f()) == EQUAL );
            return (x - y).sign();
        }
    };
    
    class To_double 
        : public CGAL::unary_function< Real_embeddable, double >{
    public:
        double operator()(const Real_embeddable& x) const {
            return x.to_double();
        }
    };
    
    class To_interval 
        : public CGAL::unary_function< Real_embeddable, std::pair< double, double > >{
    public:
        std::pair<double,double> operator()(const Real_embeddable& x) const {
            return x.to_interval();
        }
    };   
};

// operator <<
template<class Stream, class NT>
inline
Stream&
operator<<(Stream& os, const Sqrt_extension_2<NT>& x)
{
  os << "(" << x.a()  << ")+(" << x.b() << ") sqrt{" << x.e() << "}";
  os << "+(" << x.c() << ") sqrt{" << x.f() << "}";
  os << "+(" << x.d() << ") sqrt{(" << x.e() << ") (" << x.f() 
     << ")}";
  return os;
}






} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_SQRT_EXTENSION_2_H
