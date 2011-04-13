// ============================================================================
//
// Copyright (c) 1999 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------------
// 
// release       :
// release_date  :
// 
// file          : user_classes.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Herve Broennimann
//                 Stefan Schirra
//
// coordinator   : MPI, Saarbruecken  (<Stefan.Schirra@mpi-sb.mpg.de>)
// ============================================================================
 

#ifndef CGAL_USER_CLASSES_H
#define CGAL_USER_CLASSES_H

namespace CGAL {
/* 2D */
template < class R, class T = typename R::Rep_tag > class Point_2;
template < class R, class T = typename R::Rep_tag > class Vector_2;
template < class R, class T = typename R::Rep_tag > class Direction_2;
template < class R, class T = typename R::Rep_tag > class Line_2;
template < class R, class T = typename R::Rep_tag > class Ray_2;
template < class R, class T = typename R::Rep_tag > class Segment_2;
template < class R, class T = typename R::Rep_tag > class Triangle_2;
template < class R, class T = typename R::Rep_tag > class Circle_2;
template < class R, class T = typename R::Rep_tag > class Data_accessor_2;
template < class R, class T = typename R::Rep_tag > class Iso_rectangle_2;
template < class R, class T = typename R::Rep_tag > class Aff_transformation_2;

/* 3D */
template < class R, class T = typename R::Rep_tag > class Plane_3;
template < class R, class T = typename R::Rep_tag > class Point_3;
template < class R, class T = typename R::Rep_tag > class Vector_3;
template < class R, class T = typename R::Rep_tag > class Direction_3;
template < class R, class T = typename R::Rep_tag > class Line_3;
template < class R, class T = typename R::Rep_tag > class Ray_3;
template < class R, class T = typename R::Rep_tag > class Segment_3;
template < class R, class T = typename R::Rep_tag > class Triangle_3;
template < class R, class T = typename R::Rep_tag > class Tetrahedron_3;
template < class R, class T = typename R::Rep_tag > class Iso_cuboid_3;
template < class R, class T = typename R::Rep_tag > class Aff_transformation_3;

/* dD */
template < class R> class Point_d;

inline
bool
advanced_kernel_enabled()
{
#ifndef CGAL_CFG_NO_ADVANCED_KERNEL
  return true;
#else
  return false;
#endif // CGAL_CFG_NO_ADVANCED_KERNEL
}

} // namespace CGAL
#endif // CGAL_USER_CLASSES_H
