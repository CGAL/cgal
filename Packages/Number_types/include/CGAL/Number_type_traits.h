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
// file          : include/CGAL/Number_type_traits.h
// package       : Number_types (4.46)
// maintainer    : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Susan Hert, Michael Hoffmann
//
// coordinator   : MPI, Saarbruecken
// ======================================================================
 

#ifndef CGAL_NUMBER_TYPE_TRAITS_H
#define CGAL_NUMBER_TYPE_TRAITS_H

CGAL_BEGIN_NAMESPACE

template < class NT >
struct Number_type_traits {
  typedef typename NT::Has_gcd       Has_gcd;
  typedef typename NT::Has_division  Has_division;
  typedef typename NT::Has_sqrt      Has_sqrt;
};

template < class Rational >
struct Rational_traits {
  typedef typename Rational::NT RT;

  RT numerator   (const Rational & r) const { return r.numerator(); }
  RT denominator (const Rational & r) const { return r.denominator(); }
  
  Rational make_rational(const RT & n, const RT & d) const
  { return Rational(n, d); } 
};

// number type tags
struct Ring_tag {};
struct Euclidean_ring_tag {};
struct Field_tag {};
struct Sqrt_field_tag {};

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_TYPE_TRAITS_H
