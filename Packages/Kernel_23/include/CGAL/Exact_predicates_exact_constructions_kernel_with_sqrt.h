// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h
// package       : Kernel_23
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas, Sylvain Pion
//
// coordinator   :
//
// ======================================================================

#ifndef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_WITH_SQRT_H
#define CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_WITH_SQRT_H

#include <CGAL/Simple_cartesian.h>

#ifdef CGAL_USE_LEDA
#  include <CGAL/leda_real.h>
#elif defined CGAL_USE_CORE
#  include <CGAL/CORE_Expr.h>
#else
#  error "You need LEDA or CORE installed."
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_LEDA
typedef Simple_cartesian< leda_real >
        Exact_predicates_exact_constructions_kernel_with_sqrt;
#elif defined CGAL_USE_CORE
typedef Simple_cartesian< CORE::Expr >
        Exact_predicates_exact_constructions_kernel_with_sqrt;
#endif

CGAL_END_NAMESPACE

#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_WITH_SQRT_H
