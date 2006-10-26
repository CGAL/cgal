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
// $URL$
// $Id$
// 
//
// Author(s)     : Susan Hert, Michael Hoffmann
 

#ifndef CGAL_NUMBER_TYPE_TRAITS_H
#define CGAL_NUMBER_TYPE_TRAITS_H

#include <CGAL/basic.h>

//#include <CGAL/Algebraic_structure_traits.h>

CGAL_BEGIN_NAMESPACE

/*
// This is a proposal for backward compatibility. 
// there is a problem with the division. 
// The old definition defines it as Tag_true even for integers!, 
// while the use seems to be meant for MP_float only 
// The new implementation thus defines has_division as true for 
// Fields only. 
// (NOTE: MP_float is not concidered as a Field in the new layout, AST)
template < class NT >
struct Number_type_traits {
private:
    typedef Algebraic_structure_traits<NT> AST; 
    static const bool is_exact =   ::boost::is_same< CGAL::Tag_true, typename AST::Is_exact>::value;
    static const bool has_gcd  = ! ::boost::is_same< CGAL::Null_functor, typename AST::Gcd>::value;
    static const bool has_sqrt = ! ::boost::is_same< CGAL::Null_functor, typename AST::Sqrt >::value;

    static const bool has_division  = 
        ::boost::is_base_and_derived<CGAL::Field_tag , typename AST::Algebraic_structure_tag>::value 
    ||  ::boost::is_same<CGAL::Field_tag, typename AST::Algebraic_structure_tag>::value ;
    
    static const bool has_ring_operations = 
        ::boost::is_base_and_derived<CGAL::Integral_domain_without_division_tag , typename AST::Algebraic_structure_tag>::value 
    ||  ::boost::is_same<CGAL::Integral_domain_without_division_tag, typename AST::Algebraic_structure_tag>::value ;
public:
    typedef typename ::boost::mpl::if_c<has_gcd     ,CGAL::Tag_true, CGAL::Tag_false>::type 
    Has_gcd;
    typedef typename ::boost::mpl::if_c<has_division,CGAL::Tag_true, CGAL::Tag_false>::type 
    Has_division;
    typedef typename ::boost::mpl::if_c<has_sqrt    ,CGAL::Tag_true, CGAL::Tag_false>::type 
    Has_sqrt;

    typedef typename ::boost::mpl::if_c<is_exact && has_ring_operations ,CGAL::Tag_true, CGAL::Tag_false>::type 
    Has_exact_ring_operations;
    typedef typename ::boost::mpl::if_c<is_exact && has_division        ,CGAL::Tag_true, CGAL::Tag_false>::type 
    Has_exact_division;
    typedef typename ::boost::mpl::if_c<is_exact && has_sqrt            ,CGAL::Tag_true, CGAL::Tag_false>::type 
    Has_exact_sqrt;
};
*/

template < class NT >
struct Number_type_traits {
  typedef typename NT::Has_exact_ring_operations  Has_exact_ring_operations;
  typedef typename NT::Has_exact_division         Has_exact_division;
  typedef typename NT::Has_exact_sqrt             Has_exact_sqrt;

  typedef typename NT::Has_gcd       Has_gcd;
  typedef typename NT::Has_division  Has_division;
  typedef typename NT::Has_sqrt      Has_sqrt;
};

template < class Rational >
struct Rational_traits {
  typedef Rational RT;

  const RT& numerator   (const Rational& r) const { return r; }
  RT denominator (const Rational&) const { return RT(1); }
  
  Rational make_rational(const RT & n, const RT & d) const
  { return n / d; }
};

// number type tags
struct Ring_tag {};
//struct Euclidean_ring_tag {};
//struct Field_tag {};
struct Sqrt_field_tag {};

CGAL_END_NAMESPACE

#endif // CGAL_NUMBER_TYPE_TRAITS_H
