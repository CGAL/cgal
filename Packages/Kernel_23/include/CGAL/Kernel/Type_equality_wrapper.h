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
// file          : include/CGAL/Kernel/Type_equality_wrapper.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis
//
// ======================================================================

#ifndef CGAL_KERNEL_TYPE_EQUALITY_WRAPPER_H
#define CGAL_KERNEL_TYPE_EQUALITY_WRAPPER_H

#include <CGAL/user_classes.h>

CGAL_BEGIN_NAMESPACE

// This is a kernel wrapper which provides the type equality between
// Kernel::Point_2 and CGAL::Point_2<Kernel>
// (and similar for the other types).

template < typename K_base >
class Type_equality_wrapper
  : public K_base
{
    typedef typename K_base::Kernel                 K_;    

public:

    typedef K_base                                  Kernel_base;    

    typedef CGAL::Point_2<K_>                       Point_2;
    typedef CGAL::Vector_2<K_>                      Vector_2;
    typedef CGAL::Direction_2<K_>                   Direction_2;
    typedef CGAL::Line_2<K_>                        Line_2;
    typedef CGAL::Ray_2<K_>                         Ray_2;
    typedef CGAL::Segment_2<K_>                     Segment_2;
    typedef CGAL::Triangle_2<K_>                    Triangle_2;
    typedef CGAL::Circle_2<K_>                      Circle_2;
    typedef CGAL::Iso_rectangle_2<K_>               Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<K_>          Aff_transformation_2;

    typedef CGAL::Point_3<K_>                       Point_3;
    typedef CGAL::Vector_3<K_>                      Vector_3;
    typedef CGAL::Direction_3<K_>                   Direction_3;
    typedef CGAL::Line_3<K_>                        Line_3;
    typedef CGAL::Plane_3<K_>                       Plane_3;
    typedef CGAL::Ray_3<K_>                         Ray_3;
    typedef CGAL::Segment_3<K_>                     Segment_3;
    typedef CGAL::Triangle_3<K_>                    Triangle_3;
    typedef CGAL::Tetrahedron_3<K_>                 Tetrahedron_3;
    typedef CGAL::Sphere_3<K_>                      Sphere_3;
    typedef CGAL::Iso_cuboid_3<K_>                  Iso_cuboid_3;
    typedef CGAL::Aff_transformation_3<K_>          Aff_transformation_3;
};

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_TYPE_EQUALITY_WRAPPER_H
