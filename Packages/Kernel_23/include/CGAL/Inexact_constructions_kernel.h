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
// file          : include/CGAL/Inexact_constructions_kernel.h
// package       : Kernel_23
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================

#ifndef CGAL_INEXACT_CONSTRUCTIONS_KERNEL_H
#define CGAL_INEXACT_CONSTRUCTIONS_KERNEL_H


#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Filtered_kernel.h>

CGAL_BEGIN_NAMESPACE


typedef Filtered_kernel< Simple_cartesian<double> >
/*                                       */ Inexact_constructions_kernel;


CGAL_END_NAMESPACE


#endif // CGAL_INEXACT_CONSTRUCTIONS_KERNEL_H


