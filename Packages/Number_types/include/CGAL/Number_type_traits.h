// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
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
// Author(s)     : Susan Hert, Michael Hoffmann
 

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
