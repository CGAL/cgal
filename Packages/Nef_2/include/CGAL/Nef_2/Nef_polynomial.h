// Copyright (c) 1997-2000  Max-Planck-Institute Saarbruecken (Germany).
// All rights reserved.
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
// Author(s)     : Michael Seel <seel@mpi-sb.mpg.de>

#ifndef CGAL_NEF_POLYNOMIAL_H
#define CGAL_NEF_POLYNOMIAL_H

#include <CGAL/Nef_2/Polynomial.h>
//#include <CGAL/basic.h>
//#include <CGAL/kernel_assertions.h>
//#include <CGAL/Handle_for.h>
//#include <CGAL/number_type_basic.h>
//#include <CGAL/number_utils.h>
//#include <CGAL/Number_type_traits.h>
//#include <CGAL/IO/io.h>
#include <cstddef>
#undef CGAL_NEF_DEBUG
#define CGAL_NEF_DEBUG 3
#include <CGAL/Nef_2/debug.h>
#include <vector>

#include <CGAL/Kernel/mpl.h>

#include <boost/operators.hpp>

CGAL_BEGIN_NAMESPACE

#define CGAL_int(T)    typename First_if_different<int,    T>::Type
#define CGAL_double(T) typename First_if_different<double, T>::Type

template <class NT> 
class Nef_polynomial
  : boost::ordered_field_operators1< Nef_polynomial<NT>
  , boost::ordered_field_operators2< Nef_polynomial<NT>, int
  > >
  , public Polynomial<NT>
{
  typedef typename CGAL::Polynomial<NT>  Base;
  typedef typename Base::size_type       size_type;

 protected:
  Nef_polynomial(size_type s) : Base(s) {}

 public:
  Nef_polynomial() : Base() {}
  Nef_polynomial(const NT& a0) : Base(a0) {}
  Nef_polynomial(const NT& a0, const NT& a1) : Base(a0,a1) {}
  Nef_polynomial(const NT& a0, const NT& a1, const NT& a2) : Base(a0,a1,a2) {}

  template <class Fwd_iterator>
  Nef_polynomial(std::pair<Fwd_iterator, Fwd_iterator> poly) : Base(poly) {}

  Nef_polynomial(CGAL_double(NT) n) : Base(n) {}
  Nef_polynomial(CGAL_double(NT) n1, CGAL_double(NT) n2) : Base(n1, n2) {}
  Nef_polynomial(CGAL_int(NT) n) : Base(NT(n)) {}
  Nef_polynomial(CGAL_int(NT) n1, CGAL_int(NT) n2) : Base(n1,n2) {}

  Nef_polynomial(const Base& p) : Base(p) {}

  Base & polynomial() { return static_cast<Base&>(*this); }
  const Base & polynomial() const  { return static_cast<const Base&>(*this); }

    //  static NT R_; // for visualization only
    //  static void set_R(const NT& R) { R_ = R; }
    static NT& infi_maximal_value() {
      static NT R_ = 1;
      return R_;
    }
};

template <class NT> 
inline
Nef_polynomial<NT> operator-(const Nef_polynomial<NT> &a)
{
  return - a.polynomial();
}

template <class NT> 
inline
bool operator<(const Nef_polynomial<NT> &a, const Nef_polynomial<NT> &b)
{
  return a.polynomial() < b.polynomial();
}

template <class NT> 
inline
bool operator==(const Nef_polynomial<NT> &a, const Nef_polynomial<NT> &b)
{
  return a.polynomial() == b.polynomial();
}

template <class NT> 
inline
bool operator==(const Nef_polynomial<NT> &a, int b)
{
  return a.polynomial() == b;
}

template <class NT> 
inline
bool operator<(const Nef_polynomial<NT> &a, int b)
{
  return a.polynomial() < b;
}

template <class NT> 
inline
bool operator>(const Nef_polynomial<NT> &a, int b)
{
  return a.polynomial() > b;
}


// template <class NT> NT Nef_polynomial<NT>::R_ = 1;
// int                    Nef_polynomial<int>::R_ = 1;
// double                 Nef_polynomial<double>::R_ = 1.0;

template <class NT>
double to_double(const Nef_polynomial<NT>& p)
{
  return CGAL::to_double(p.eval_at(Nef_polynomial<NT>::infi_maximal_value()));
}

template <class NT>
Nef_polynomial<NT>
gcd(const Nef_polynomial<NT>& p1, const Nef_polynomial<NT>& p2)
{
  return gcd(p1.polynomial(), p2.polynomial());
}

#undef CGAL_double
#undef CGAL_int

CGAL_END_NAMESPACE

#endif  // CGAL_NEF_POLYNOMIAL_H
