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
// file          : include/CGAL/Exact_kernel.h
// package       : Kernel_23
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   :
//
// ======================================================================

#ifndef CGAL_EXACT_KERNEL_H
#define CGAL_EXACT_KERNEL_H


#include <CGAL/basic.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Lazy_nt_exact.h>
#include <CGAL/Quotient.h>
#include <CGAL/MP_Float.h>

CGAL_BEGIN_NAMESPACE


typedef Simple_cartesian<  Lazy_exact_nt< Quotient<MP_Float> >  >
/*                                                       */ Exact_kernel;


CGAL_END_NAMESPACE


#endif // CGAL_EXACT_KERNEL_H


