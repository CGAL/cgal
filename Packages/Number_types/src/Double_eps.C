// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
// 
// release       : 
// release_date  : 2000, December 13
// 
// source        : Double_eps.fw
// file          : Double_eps.C
// package       : Number_types (4.2)
// maintainer    : Stefan Schirra <stschirr@mpi-sb.mpg.de>
// revision      : 4.2
// revision_date : 13 Dec 2000 
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ======================================================================
 

#ifndef CGAL_BASIC_H
#include <CGAL/basic.h>
#endif // CGAL_BASIC_H
#ifndef CGAL_DOUBLE_EPS_H
#include <CGAL/Double_eps.h>
#endif // CGAL_DOUBLE_EPS_H

CGAL_BEGIN_NAMESPACE

double Double_eps::_eps = 0.0;


double set_eps(double eps)
  {
    double e = Double_eps::_eps;
    Double_eps::_eps = eps;
    return e;
  }

Double_eps sqrt(const Double_eps &de)
  {
    return Double_eps(sqrt(de.d()));
  }

std::ostream&
operator<<(std::ostream& os, const Double_eps &de)
{
  os << de.d();
  return os;
}

std::istream&
operator>>(std::istream& is, Double_eps &de)
{
  is >> de._d;
  return is;
}
CGAL_END_NAMESPACE

