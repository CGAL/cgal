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
// file          : include/CGAL/Exact_kernel_with_sqrt.h
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



#if defined CGAL_USE_LEDA

#include <CGAL/leda_real.h>

CGAL_BEGIN_NAMESPACE

typedef Simple_cartesian<  Lazy_exact_nt< leda_real >  >
/*                                             */ Exact_kernel_with_sqrt;

CGAL_END_NAMESPACE

#elif

#include <CORE/CGAL_Expr.h>

CGAL_BEGIN_NAMESPACE

typedef Simple_cartesian<  Lazy_exact_nt< CORE::Expr >  >
/*                                             */ Exact_kernel_with_sqrt;

CGAL_END_NAMESPACE

#endif



#endif // CGAL_EXACT_KERNEL_H


