// Copyright (c) 1999  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
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
// Author(s)     : Stefan Schirra, Sylvain Pion
 
#ifndef CGAL_SIMPLE_HOMOGENEOUS_H
#define CGAL_SIMPLE_HOMOGENEOUS_H

#include <CGAL/Homogeneous/Homogeneous_base.h>
#include <CGAL/Simple_Handle_for.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/Quotient.h>

CGAL_BEGIN_NAMESPACE

template < typename RT_, typename FT_, typename Kernel >
struct Homogeneous_base_no_ref_count
  : public Homogeneous_base< Kernel >
{
    typedef RT_                                           RT;
    typedef FT_                                           FT;

    // The mecanism that allows to specify reference-counting or not.
    template < typename T >
    struct Handle { typedef Simple_Handle_for<T>    type; };

    template < typename Kernel2 >
    struct Base {
        typedef Homogeneous_base_no_ref_count<RT_,FT_,Kernel2> Type;
    };

    // TODO: cleanup (use Rational_traits<> instead)
    static FT make_FT(const RT & num, const RT& denom)
    { return FT(num, denom); }

    static FT make_FT(const RT & num)
    { return FT(num); }

    static RT FT_numerator(const FT &r)
    { return r.numerator(); }

    static RT FT_denominator(const FT &r)
    { return r.denominator(); }
};

template < typename RT_, typename FT_ = Quotient<RT_> >
struct Simple_homogeneous
  : public Type_equality_wrapper<
                Homogeneous_base_no_ref_count<RT_, FT_,
                                              Simple_homogeneous<RT_, FT_> >,
                Simple_homogeneous<RT_, FT_> >
{};

CGAL_END_NAMESPACE

CGAL_ITERATOR_TRAITS_POINTER_SPEC_TEMPLATE(CGAL::Simple_homogeneous)

#endif // CGAL_SIMPLE_HOMOGENEOUS_H
