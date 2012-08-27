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
// 
//
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>


#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_SQRT_EXTENSION_1_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_SQRT_EXTENSION_1_H

#include <CGAL/basic.h>
#include <CGAL/enum.h>

#include <iostream>

namespace CGAL {

template<class NT>
class Sqrt_extension_1;

template<class NT>
class Sqrt_extension_1
{
private:
  NT x, y, r;

private:
  typedef Sqrt_extension_1<NT>   Self;

public:
  typedef NT                  FT;
  typedef NT                  RT;

  Sqrt_extension_1() : x(0), y(0), r(0) {}
  Sqrt_extension_1(int i) : x(i), y(0), r(0) {}
  Sqrt_extension_1(const NT& a) : x(a), y(0), r(0) {}
  Sqrt_extension_1(const NT& a, const NT& b, const NT& c)
    : x(a), y(b), r(c)
  {
    CGAL_exactness_assertion( !(CGAL::is_negative(r)) );
  }

  Sqrt_extension_1(const Sqrt_extension_1<NT>& other)
    : x(other.x), y(other.y), r(other.r) {}


  NT a() const { return x; }
  NT b() const { return y; }
  NT c() const { return r; }

  Self square() const
  {
    NT xy = x * y;

    return Self(CGAL::square(x) + CGAL::square(y) * r,
		xy + xy,
		r);
  }

  Sign sign() const
  {
    Sign sx = CGAL::sign(x);

    if ( CGAL::sign(r) == ZERO )  { return sx; }

    Sign sy = CGAL::sign(y);

    if ( sx == sy )  { return sx; }
    if ( sx == ZERO )  { return sy; }

    return sx * CGAL::compare( CGAL::square(x),
			       r * CGAL::square(y) );
  }

  const Self& operator+() const
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
    double xd = CGAL::to_double(x);
    double yd = CGAL::to_double(y);
    double rd = CGAL::to_double(r);

    return (xd + yd * CGAL::sqrt(rd));
  }

  std::pair<double,double> to_interval() const
  {
    // THIS MUST BE CHECK WITH SYLVAIN FOR CORRECTNESS
    std::pair<double,double> x_ivl = CGAL::to_interval(x);
    std::pair<double,double> y_ivl = CGAL::to_interval(y);
    std::pair<double,double> r_ivl = CGAL::to_interval(r);

    std::pair<double,double> sqrt_r_ivl(CGAL::sqrt(r_ivl.first),
					CGAL::sqrt(r_ivl.second));

    std::pair<double,double>
      ivl(x_ivl.first + y_ivl.first * sqrt_r_ivl.first,
	  x_ivl.second + y_ivl.second * sqrt_r_ivl.second);

    return ivl;
  }
};


// operator *
template<class NT>
inline
Sqrt_extension_1<NT>
operator*(const Sqrt_extension_1<NT>& x, const NT& n)
{
  return Sqrt_extension_1<NT>(x.a() * n, x.b() * n, x.c());
}


template<class NT>
inline
Sqrt_extension_1<NT>
operator*(const NT& n, const Sqrt_extension_1<NT>& x)
{
  return (x * n);
}

template<class NT>
inline
Sqrt_extension_1<NT>
operator*(const Sqrt_extension_1<NT>& x, const Sqrt_extension_1<NT>& y)
{
  CGAL_exactness_precondition( CGAL::compare(x.c(), y.c()) == EQUAL );

  NT a = x.a() * y.a() + x.b() * y.b() * x.c();
  NT b = x.a() * y.b() + x.b() * y.a();

  return Sqrt_extension_1<NT>(a, b, x.c());
}


// operator +
template<class NT>
inline
Sqrt_extension_1<NT>
operator+(const Sqrt_extension_1<NT>& x, const NT& n)
{
  return Sqrt_extension_1<NT>(x.a() + n, x.b(), x.c());
}


template<class NT>
inline
Sqrt_extension_1<NT>
operator+(const NT& n, const Sqrt_extension_1<NT>& x)
{
  return (x + n);
}

template<class NT>
inline
Sqrt_extension_1<NT>
operator+(const Sqrt_extension_1<NT>& x, const Sqrt_extension_1<NT>& y)
{
  CGAL_exactness_precondition( CGAL::compare(x.c(), y.c()) == EQUAL );

  return Sqrt_extension_1<NT>(x.a() + y.a(), x.b() + y.b(), x.c());
}



// operator -
template<class NT>
inline
Sqrt_extension_1<NT>
operator-(const Sqrt_extension_1<NT>& x, const NT& n)
{
  return x + (-n);
}


template<class NT>
inline
Sqrt_extension_1<NT>
operator-(const NT& n, const Sqrt_extension_1<NT>& x)
{
  return -(x - n);
}

template<class NT>
inline
Sqrt_extension_1<NT>
operator-(const Sqrt_extension_1<NT>& x, const Sqrt_extension_1<NT>& y)
{
  return (x + (-y));
}


//=============================================================

template <class NT> 
class Algebraic_structure_traits<Sqrt_extension_1<NT> >
    :public Algebraic_structure_traits_base<Sqrt_extension_1<NT>,CGAL::Integral_domain_without_division_tag>{
    // I haven't found division 
private:
    typedef Algebraic_structure_traits<NT> AST_NT;
public:
    typedef Sqrt_extension_1<NT> Algebraic_structure;
    typedef typename AST_NT::Is_exact Is_exact;
};

template<class NT>
class Real_embeddable_traits<Sqrt_extension_1<NT> >{
private:
    typedef Real_embeddable_traits<NT> RET_NT;
public:
    
    typedef Sqrt_extension_1<NT> Real_embeddable;
    
    class Abs 
        : public std::unary_function< Real_embeddable, Real_embeddable >{
    public:
        Real_embeddable operator()(const Real_embeddable& x) const {
            return (x>=0)?x:-x;
        }
    };    

    class Sgn 
        : public std::unary_function< Real_embeddable, CGAL::Sign >{
    public:
        CGAL::Sign operator()(const Real_embeddable& x) const {
            return x.sign();
        }
    };
    
    class Compare 
        : public std::binary_function< Real_embeddable, 
                                  Real_embeddable, 
                                  CGAL::Comparison_result >{
    public:
        CGAL::Comparison_result operator()(
                const Real_embeddable& x, 
                const Real_embeddable& y) const {
            CGAL_exactness_precondition( CGAL::compare(x.c(), y.c()) == EQUAL );
            return (x - y).sign();
            
// This is not needed due to equality of CGAL::Sign CGAL::Comparison_result
//             CGAL::Sign s = (x - y).sign();
//             if ( s == ZERO ) { return EQUAL; }
//             return (s == POSITIVE) ? LARGER : SMALLER;
        }
    };
    
    class To_double 
        : public std::unary_function< Real_embeddable, double >{
    public:
        double operator()(const Real_embeddable& x) const {
            return x.to_double();
        }
    };
    
    class To_interval 
        : public std::unary_function< Real_embeddable, std::pair< double, double > >{
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
operator<<(Stream& os, const Sqrt_extension_1<NT>& x)
{
  os << "(" << x.a()  << ")+(" << x.b() << ") sqrt{" << x.c() << "}";
  return os;
}



} //namespace CGAL



#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_SQUARE_ROOT_1_H
