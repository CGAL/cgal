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
// file          : include/CGAL/Exact_predicates_exact_constructions_kernel.h
// package       : Kernel_23
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas, Sylvain Pion
//
// coordinator   :
//
// ======================================================================

#ifndef CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
#define CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Lazy_nt_exact.h>

#ifdef CGAL_USE_GMP
#  include <CGAL/Gmpq.h>
#else
#  include <CGAL/Quotient.h>
#  include <CGAL/MP_Float.h>
#endif

CGAL_BEGIN_NAMESPACE

#ifdef CGAL_USE_GMP
typedef Simple_cartesian< Lazy_exact_nt< Gmpq > >               Exact_kernel;
#else
typedef Simple_cartesian< Lazy_exact_nt< Quotient<MP_Float> > > Exact_kernel;
#endif

CGAL_END_NAMESPACE

#endif // CGAL_EXACT_PREDICATES_EXACT_CONSTRUCTIONS_KERNEL_H
