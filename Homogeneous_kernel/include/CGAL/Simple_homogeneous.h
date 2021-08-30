// Copyright (c) 1999
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra, Sylvain Pion

#ifndef CGAL_SIMPLE_HOMOGENEOUS_H
#define CGAL_SIMPLE_HOMOGENEOUS_H

#include <CGAL/Homogeneous/Homogeneous_base.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>
#include <CGAL/Quotient.h>

namespace CGAL {

template < typename RT_, typename FT_, typename Kernel >
struct Homogeneous_base_no_ref_count
  : public Homogeneous_base<RT_, FT_, Kernel >
{
    typedef RT_                                           RT;
    typedef FT_                                           FT;

    // The mechanism that allows to specify reference-counting or not.
    template < typename T >
    struct Handle { typedef T    type; };

    template < typename Kernel2 >
    struct Base {
        typedef Homogeneous_base_no_ref_count<RT_,FT_,Kernel2> Type;
    };
};

template < typename RT_, typename FT_ = Quotient<RT_> >
struct Simple_homogeneous
  : public Type_equality_wrapper<
                Homogeneous_base_no_ref_count<RT_, FT_,
                                              Simple_homogeneous<RT_, FT_> >,
                Simple_homogeneous<RT_, FT_> >
{};

} //namespace CGAL

#endif // CGAL_SIMPLE_HOMOGENEOUS_H
