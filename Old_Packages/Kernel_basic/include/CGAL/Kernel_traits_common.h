// ======================================================================
//
// Copyright (c) 2000 The CGAL Consortium
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
// file          : include/CGAL/Kernel_traits_common.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Bronnimann, Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (Mariette.Yvinec@sophia.inria.fr)
//
// ======================================================================

// This file is intentionally not protected for re-inclusion.
// It's aimed at being included from within a kernel traits class, this
// way we share more code.

#define CGAL_Kernel_pred(X,Y,Z) typedef X Y; Y Z() const { return Y(); }
#define CGAL_Kernel_cons(X,Y,Z) CGAL_Kernel_pred(X,Y,Z)
#define CGAL_Kernel_pred2(W,X,Y,Z) typedef W,X Y; Y Z() const { return Y(); }
#define CGAL_Kernel_cons2(W,X,Y,Z) CGAL_Kernel_pred2(W,X,Y,Z)

#include <CGAL/Kernel/interface_macros.h>

// The following are like constructive predicate, and are deprecated !
#ifndef CGAL_NO_DEPRECATED_CODE

typedef CGAL ::p_Less_dist_to_line_2p<Point_2>  Less_signed_distance_to_line_2;
Less_signed_distance_to_line_2
less_signed_distance_to_line_2_object(const Point_2& p,
				      const Point_2& q)  const   // XXX XXX
{ return Less_signed_distance_to_line_2(p,q); }

typedef CGAL ::p_Less_rotate_ccw<Point_2>          Less_rotate_ccw_2;
Less_rotate_ccw_2
less_rotate_ccw_2_object(const Point_2& p) const   // XXX XXX
{ return Less_rotate_ccw_2(p); }

typedef CGAL ::p_Left_of_line_2p<Point_2>          Left_of_line_2;
Left_of_line_2
left_of_line_2_object(const Point_2& p, const Point_2& q) const  // XXX XXX
{ return Left_of_line_2(p,q); }

#endif // CGAL_NO_DEPRECATED_CODE

