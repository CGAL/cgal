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
// Kernel::Point_2 and CGAL::Point_2<Kernel>, by deriving from
// K_base::Point_2 (and similar for the other types).

template < typename K_base, typename Kernel >
struct Type_equality_wrapper
  : public K_base
{
    typedef K_base                                  Kernel_base;

    typedef CGAL::Point_2<Kernel>                   Point_2;
    typedef CGAL::Vector_2<Kernel>                  Vector_2;
    typedef CGAL::Direction_2<Kernel>               Direction_2;
    typedef CGAL::Line_2<Kernel>                    Line_2;
    typedef CGAL::Ray_2<Kernel>                     Ray_2;
    typedef CGAL::Segment_2<Kernel>                 Segment_2;
    typedef CGAL::Triangle_2<Kernel>                Triangle_2;
    typedef CGAL::Circle_2<Kernel>                  Circle_2;
    typedef CGAL::Iso_rectangle_2<Kernel>           Iso_rectangle_2;
    typedef CGAL::Aff_transformation_2<Kernel>      Aff_transformation_2;

    typedef CGAL::Point_3<Kernel>                   Point_3;
    typedef CGAL::Vector_3<Kernel>                  Vector_3;
    typedef CGAL::Direction_3<Kernel>               Direction_3;
    typedef CGAL::Line_3<Kernel>                    Line_3;
    typedef CGAL::Plane_3<Kernel>                   Plane_3;
    typedef CGAL::Ray_3<Kernel>                     Ray_3;
    typedef CGAL::Segment_3<Kernel>                 Segment_3;
    typedef CGAL::Triangle_3<Kernel>                Triangle_3;
    typedef CGAL::Tetrahedron_3<Kernel>             Tetrahedron_3;
    typedef CGAL::Sphere_3<Kernel>                  Sphere_3;
    typedef CGAL::Iso_cuboid_3<Kernel>              Iso_cuboid_3;
    typedef CGAL::Aff_transformation_3<Kernel>      Aff_transformation_3;
};

CGAL_END_NAMESPACE

#endif // CGAL_KERNEL_TYPE_EQUALITY_WRAPPER_H
