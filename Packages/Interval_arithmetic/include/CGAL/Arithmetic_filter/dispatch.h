// ======================================================================
//
// Copyright (c) 1999 The CGAL Consortium
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
// file          : include/CGAL/Arithmetic_filter/dispatch.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Interval_arithmetic
// author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
//
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
// ======================================================================

#ifndef CGAL_ARITHMETIC_FILTER_DISPATCH_H
#define CGAL_ARITHMETIC_FILTER_DISPATCH_H

#if !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_BUILTIN_H )
#include <CGAL/Arithmetic_filter/predicates/builtin.h>
#endif

#if defined( CGAL_PREDICATES_SIGN_OF_DETERMINANT_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_SIGN_OF_DETERMINANT_H )
#include <CGAL/Arithmetic_filter/predicates/sign_of_determinant.h>
#endif

#if defined( CGAL_PREDICATES_KERNEL_FTC2_H ) && \
   !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_KERNEL_FTC2_H )
#include <CGAL/Arithmetic_filter/predicates/kernel_ftC2.h>
#endif

#if defined( CGAL_PREDICATES_KERNEL_FTC3_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_PREDICATES_KERNEL_FTC3_H )
#include <CGAL/Arithmetic_filter/predicates/kernel_ftC3.h>
#endif

#if defined( CGAL_REGULAR_TRIANGULATION_FTC2_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC2_H )
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_ftC2.h>
#endif

    /*
#if defined( CGAL_REGULAR_TRIANGULATION_RTH2_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_RTH2_H )
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_rtH2.h>
#endif
*/

#if defined( CGAL_REGULAR_TRIANGULATION_FTC3_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_FTC3_H )
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_ftC3.h>
#endif

    /*
#if defined( CGAL_REGULAR_TRIANGULATION_RTH3_H ) && \
       !defined( CGAL_ARITHMETIC_FILTER_REGULAR_TRIANGULATION_RTH3_H )
#include <CGAL/Arithmetic_filter/predicates/Regular_triangulation_rtH3.h>
#endif
*/

#endif // CGAL_ARITHMETIC_FILTER_DISPATCH_H
