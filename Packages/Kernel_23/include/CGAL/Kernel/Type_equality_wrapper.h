// Copyright (c) 2003  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbrucken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion

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

    // Undocumented stuff.
    typedef CGAL::Conic_2<Kernel>                   Conic_2;

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
