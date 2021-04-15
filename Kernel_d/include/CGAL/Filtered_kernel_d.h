// Copyright (c) 2009-2014 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)    : Samuel Hornus, Olivier Devillers

#ifndef CGAL_FILTERED_KERNEL_D_H
#define CGAL_FILTERED_KERNEL_D_H

#include <CGAL/Filtered_predicate.h>
#include <CGAL/internal/Exact_type_selector.h>
#include <CGAL/Kernel_d/Cartesian_converter_d.h>
#include <CGAL/Cartesian_d.h>

namespace CGAL {

template<typename Kernel> // a dD kernel we want to filter
struct Filtered_kernel_d : public Cartesian_d<typename Kernel::FT>
{
    typedef typename Kernel::LA                        LA;
    typedef typename Kernel::RT                        RT; // Ring type
    typedef typename Kernel::FT                        FT; // Field type

    typedef Cartesian_d<FT>                        Base;
    typedef Filtered_kernel_d<Kernel>                Self;

    // an exact number type
    typedef typename internal::Exact_type_selector<RT>::Type        Exact_nt;

    // the corresponding exact kernel
    //typedef Linear_algebraCd< Exact_nt, boost::pool_allocator<Exact_nt> > Exact_linalg;
    typedef Linear_algebraCd< Exact_nt > Exact_linalg;
    typedef Cartesian_d<Exact_nt, Exact_linalg>                Exact_kernel;

    // the kernel used for filtered predicates
    typedef Interval_nt<false> IA;
    //typedef Linear_algebraCd<IA, boost::pool_allocator<IA> > Interval_linalg;
    typedef Linear_algebraCd<IA> Interval_linalg;
    typedef Cartesian_d<IA, Interval_linalg >        Approximate_kernel;

    // the converter
    typedef Cartesian_converter_d<Base, Exact_kernel>        C2E;
    typedef Cartesian_converter_d<Base, Approximate_kernel>        C2F;

    // we change the predicates.
#define CGAL_Kernel_pred(P, Pf) \
    typedef Filtered_predicate<typename Exact_kernel::P, typename Approximate_kernel::P, C2E, C2F> P; \
    P Pf() const { return P(); }

    // we don't touch the constructions.
#define CGAL_Kernel_cons(Y,Z)

#include <CGAL/Kernel_d/interface_macros_d.h>
};

} //namespace CGAL

#endif // CGAL_FILTERED_KERNEL_D_H
