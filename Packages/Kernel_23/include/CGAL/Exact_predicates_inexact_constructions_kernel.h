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
// file          : include/CGAL/Exact_predicates_inexact_constructions_kernel.h
// package       : Kernel_23
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas, Sylvain Pion
//
// coordinator   :
//
// ======================================================================

#ifndef CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
#define CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H

#include <CGAL/Simple_cartesian.h>
// We don't use Filtered_kernel at the moment because it's slower
// that Filtered_exact.
// #include <CGAL/Filtered_kernel.h>

#ifdef CGAL_CFG_MATCHING_BUG_2
#  define CGAL_IA_CT        double
#  define CGAL_IA_PROTECTED true
#  define CGAL_IA_CACHE     No_Filter_Cache
#  define CGAL_IA_ET        CGAL::MP_Float
#endif

#include <CGAL/MP_Float.h>
#include <CGAL/Filtered_exact.h>

CGAL_BEGIN_NAMESPACE

#if 0
typedef Filtered_kernel< Simple_cartesian<double> >
        Exact_predicates_inexact_constructions_kernel;
#else
typedef Simple_cartesian<Filtered_exact<double, MP_Float> >
        Exact_predicates_inexact_constructions_kernel;
#endif

CGAL_END_NAMESPACE

#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
