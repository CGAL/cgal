// ======================================================================
//
// Copyright (c) 1997-2000 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : $CGAL_Revision: CGAL-2.4-I-78 $
// release_date  : $CGAL_Date: 2002/04/17 $
//
// file          : include/CGAL/Nef_2/Nef_polynomial.h
// package       : Nef_2
// maintainer    : Peter Hachenberger <hachenberger@mpi-sb.mpg.de>
// chapter       : Nef Polyhedra
//
// revision      : $Revision$
// revision_date : $Date$
//
// author(s)     : Michael Seel <seel@mpi-sb.mpg.de>
// coordinator   : Michael Seel <seel@mpi-sb.mpg.de>
//
// implementation: Polynomials in one variable
// ======================================================================

#ifndef CGAL_NEF_POLYNOMIAL_H
#define CGAL_NEF_POLYNOMIAL_H

#include <CGAL/Nef_2/Polynomial.h>
#include <CGAL/basic.h>
#include <CGAL/kernel_assertions.h>
#include <CGAL/Handle_for.h>
#include <CGAL/number_type_basic.h>
#include <CGAL/number_utils.h>
#include <CGAL/Number_type_traits.h>
#include <CGAL/IO/io.h>
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

template <class NT> NT Nef_polynomial<NT>::R_;
int                    Nef_polynomial<int>::R_;
double                 Nef_polynomial<double>::R_;

template <class NT> double to_double(const Nef_polynomial<NT>& p) { 
  return (CGAL::to_double(p.eval_at(Nef_polynomial<NT>::R_))); 
}

double to_double(const Nef_polynomial<int>& p) { 
  return (CGAL::to_double(p.eval_at(Nef_polynomial<int>::R_))); 
}

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


