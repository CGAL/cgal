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
// release       : $CGAL_Revision: CGAL-2.4-I-65 $
// release_date  : $CGAL_Date: 2002/03/19 $
// 
// file          : include/CGAL/Rational_traits.h
// package       : Number_types (4.46)
// maintainer    : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Andreas Fabri
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_RATIONAL_TRAITS_H
#define CGAL_RATIONAL_TRAITS_H

CGAL_BEGIN_NAMESPACE


template <class Rational>
struct Rational_traits {
  typedef typename Rational::NT RT;

 RT numerator   (const Rational & r) const { return r.numerator(); }
 RT denominator (const Rational & r) const { return r.denominator(); }

 Rational make_rational(const RT & n, const RT & d) const
 { return Rational(n, d); } 
};


CGAL_END_NAMESPACE

#endif // CGAL_RATIONAL_TRAITS_H
