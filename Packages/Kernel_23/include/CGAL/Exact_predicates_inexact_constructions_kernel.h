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
#include <CGAL/Filtered_kernel.h>

CGAL_BEGIN_NAMESPACE

typedef Filtered_kernel< Simple_cartesian<double> >
        Exact_predicates_inexact_constructions_kernel;

CGAL_END_NAMESPACE

#endif // CGAL_EXACT_PREDICATES_INEXACT_CONSTRUCTIONS_KERNEL_H
