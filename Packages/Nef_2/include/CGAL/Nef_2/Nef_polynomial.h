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
#undef _DEBUG
#define _DEBUG 3
#include <CGAL/Nef_2/debug.h>
#include <vector>


CGAL_BEGIN_NAMESPACE

template <class NT> 
class Nef_polynomial : public Polynomial<NT> {

  typedef typename CGAL::Polynomial<NT>  Base;
  typedef typename Base::size_type       size_type;

 protected:
  Nef_polynomial(size_type s) : Base(s) {}

 public:
  Nef_polynomial() : Base() {}
  Nef_polynomial(const NT& a0) : Base(a0) {}
  Nef_polynomial(NT a0, NT a1) : Base(a0,a1) {}
  Nef_polynomial(const NT& a0, const NT& a1,const NT& a2) : Base(a0,a1,a2) {}

  template <class Fwd_iterator>
  Nef_polynomial(std::pair<Fwd_iterator, Fwd_iterator> poly) : Base(poly) {}

  Nef_polynomial(double n) : Base(n) {}
  Nef_polynomial(double n1, double n2) : Base(n1, n2) {}
  Nef_polynomial(int n) : Base(NT(n)) {}
  Nef_polynomial(int n1, int n2) : Base(n1,n2) {}

  Nef_polynomial(const Base& p) : Base(p) {}
  Nef_polynomial(const Nef_polynomial<NT>& p) : Base(p) {}

  static NT R_; // for visualization only
  static void set_R(const NT& R) { R_ = R; }
};

template <> 
class Nef_polynomial<int> : public Polynomial<int> {

  typedef CGAL::Polynomial<int>           Base;
  typedef Base::size_type                 size_type;

 protected:
  Nef_polynomial(size_type s) : Base(s) {}

 public:
  Nef_polynomial() : Base() {}
  Nef_polynomial(const int& a0) : Base(a0) {}
  Nef_polynomial(int a0, int a1) : Base(a0,a1) {}
  Nef_polynomial(const int& a0, const int& a1,const int& a2) : Base(a0,a1,a2) {}

  template <class Fwd_iterator>
  Nef_polynomial(std::pair<Fwd_iterator, Fwd_iterator> poly) : Base(poly) {}

  Nef_polynomial(double n) : Base(n) {}
  Nef_polynomial(double n1, double n2) : Base(n1, n2) {}

  Nef_polynomial(const Base& p) : Base(p) {}
  Nef_polynomial(const Nef_polynomial<int>& p) : Base(p) {}

  static int R_; // for visualization only
  static void set_R(const int& R) { R_ = R; }
};

template <> 
class Nef_polynomial<double> : public Polynomial<double> {

  typedef CGAL::Polynomial<double>  Base;
  typedef Base::size_type           size_type;

 protected:
  Nef_polynomial(size_type s) : Base(s) {}

 public:
  Nef_polynomial() : Base() {}
  Nef_polynomial(const double& a0) : Base(a0) {}
  Nef_polynomial(double a0, double a1) : Base(a0,a1) {}
  Nef_polynomial(const double& a0, const double& a1,const double& a2) 
    : Base(a0,a1,a2) {}

  template <class Fwd_iterator>
  Nef_polynomial(std::pair<Fwd_iterator, Fwd_iterator> poly) : Base(poly) {}

  Nef_polynomial(int n) : Base(NT(n)) {}
  Nef_polynomial(int n1, int n2) : Base(n1,n2) {}

  Nef_polynomial(const Base& p) : Base(p) {}
  Nef_polynomial(const Nef_polynomial<double>& p) : Base(p) {}

  static double R_; // for visualization only
  static void set_R(const double& R) { R_ = R; }
};

template <class NT> NT Nef_polynomial<NT>::R_ = 1;
int                    Nef_polynomial<int>::R_ = 1;
double                 Nef_polynomial<double>::R_ = 1.0;

template <class NT> double to_double(const ::CGAL::Nef_polynomial<NT>& p) { 
  return (CGAL::to_double(p.eval_at(CGAL::Nef_polynomial<NT>::R_))); 
}

inline
double to_double(const Nef_polynomial<int>& p) { 
  return (CGAL::to_double(p.eval_at(Nef_polynomial<int>::R_))); 
}

inline
double to_double(const Nef_polynomial<double>& p) { 
  return (CGAL::to_double(p.eval_at(Nef_polynomial<double>::R_))); 
}

template <class NT> 
Nef_polynomial<NT> 
gcd(const Nef_polynomial<NT>& p1, const Nef_polynomial<NT>& p2) {
  return gcd(Polynomial<NT>(p1), Polynomial<NT>(p2));
}

CGAL_END_NAMESPACE

#endif  // CGAL_NEF_POLYNOMIAL_H


