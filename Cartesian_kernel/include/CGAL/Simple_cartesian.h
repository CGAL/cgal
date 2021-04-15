// Copyright (c) 2000-2004
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
// Author(s)     : Herve Bronnimann, Sylvain Pion

#ifndef CGAL_SIMPLE_CARTESIAN_H
#define CGAL_SIMPLE_CARTESIAN_H

#include <CGAL/Cartesian/Cartesian_base.h>
#include <CGAL/Handle_for.h>
#include <CGAL/Kernel/Type_equality_wrapper.h>

namespace CGAL {

template < typename FT_, typename Kernel_ >
struct Cartesian_base_no_ref_count
  : public Cartesian_base< Kernel_, FT_ >
{
    typedef FT_                                           RT;
    typedef FT_                                           FT;

    // The mechanism that allows to specify reference-counting or not.
    template < typename T >
    struct Handle { typedef T   type; };

    template < typename Kernel2 >
    struct Base { typedef Cartesian_base_no_ref_count<FT_, Kernel2>  Type; };

    typedef Kernel_ K;
#define CGAL_Kernel_pred(Y,Z) typedef CartesianKernelFunctors::Y<K> Y; \
                              Y Z() const { return Y(); }
#define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

#include <CGAL/Kernel/interface_macros.h>
};

template < typename FT_ >
struct Simple_cartesian
  : public Type_equality_wrapper<
                Cartesian_base_no_ref_count<FT_, Simple_cartesian<FT_> >,
                Simple_cartesian<FT_> >
{};

} //namespace CGAL

#endif // CGAL_SIMPLE_CARTESIAN_H
