// Copyright (c) 1999-2004
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stefan Schirra, Sylvain Pion

#ifndef CGAL_HOMOGENEOUS_BASE_H
#define CGAL_HOMOGENEOUS_BASE_H

#include <CGAL/config.h>
#include <CGAL/basic_classes.h>

#include <CGAL/Kernel/global_functions.h>

#include <CGAL/Homogeneous/Aff_transformationH2.h>
#include <CGAL/Cartesian/Circle_2.h>

//#include <CGAL/Cartesian/Direction_2.h>
#include <CGAL/Homogeneous/DirectionH2.h>
#include <CGAL/Homogeneous/Iso_rectangleH2.h>
#include <CGAL/Homogeneous/LineH2.h>
#include <CGAL/Homogeneous/PointH2.h>
#include <CGAL/Homogeneous/Weighted_point_2.h>
#include <CGAL/Cartesian/Ray_2.h>
#include <CGAL/Cartesian/Segment_2.h>
#include <CGAL/Cartesian/Triangle_2.h>
#include <CGAL/Homogeneous/VectorH2.h>
#include <CGAL/Homogeneous/Data_accessorH2.h>
#include <CGAL/Homogeneous/ConicHPA2.h>

#include <CGAL/Homogeneous/Aff_transformationH3.h>
#include <CGAL/Homogeneous/DirectionH3.h>
#include <CGAL/Homogeneous/Iso_cuboidH3.h>
#include <CGAL/Cartesian/Line_3.h>
#include <CGAL/Homogeneous/PlaneH3.h>
#include <CGAL/Homogeneous/PointH3.h>
#include <CGAL/Homogeneous/Weighted_point_3.h>
#include <CGAL/Homogeneous/RayH3.h>
#include <CGAL/Cartesian/Segment_3.h>
#include <CGAL/Homogeneous/SphereH3.h>
#include <CGAL/Cartesian/Tetrahedron_3.h>
#include <CGAL/Cartesian/Triangle_3.h>
#include <CGAL/Cartesian/Circle_3.h>
#include <CGAL/Homogeneous/VectorH3.h>

#include <CGAL/Homogeneous/basic_constructionsH2.h>
#include <CGAL/Homogeneous/distance_predicatesH2.h>
#include <CGAL/Homogeneous/predicates_on_directionsH2.h>
#include <CGAL/Homogeneous/predicates_on_pointsH2.h>

#include <CGAL/Homogeneous/basic_constructionsH3.h>
#include <CGAL/Homogeneous/distance_predicatesH3.h>
#include <CGAL/Homogeneous/predicates_on_pointsH3.h>
#include <CGAL/Homogeneous/predicates_on_pointsH2.h>

#include <CGAL/representation_tags.h>
#include <CGAL/Homogeneous/function_objects.h>

#include <CGAL/Kernel_d/Cartesian_const_iterator_d.h>

namespace CGAL {

template <typename RT_, typename FT_, typename K_ >
struct Homogeneous_base
{
    typedef K_                                      Kernel;
    typedef RT_                                     RT;
    typedef FT_                                     FT;

    typedef Homogeneous_tag                         Rep_tag;
    typedef Homogeneous_tag                         Kernel_tag;

    enum { Has_filtered_predicates = false };
    typedef Boolean_tag<Has_filtered_predicates> Has_filtered_predicates_tag;

    typedef CGAL::Object                            Object_2;
    typedef CGAL::Object                            Object_3;

    typedef typename Same_uncertainty_nt<bool, FT>::type
                                                    Boolean;
    typedef typename Same_uncertainty_nt<CGAL::Sign, FT>::type
                                                    Sign;
    typedef typename Same_uncertainty_nt<CGAL::Comparison_result, FT>::type
                                                    Comparison_result;
    typedef typename Same_uncertainty_nt<CGAL::Orientation, FT>::type
                                                    Orientation;
    typedef typename Same_uncertainty_nt<CGAL::Oriented_side, FT>::type
                                                    Oriented_side;
    typedef typename Same_uncertainty_nt<CGAL::Bounded_side, FT>::type
                                                    Bounded_side;
    typedef typename Same_uncertainty_nt<CGAL::Angle, FT>::type
                                                    Angle;

    template <typename T>
    struct Ambient_dimension {
      typedef typename T::Ambient_dimension type;
    };

    template <typename T>
    struct Feature_dimension {
      typedef typename T::Feature_dimension type;
    };

    typedef PointH2<Kernel>                         Point_2;
    typedef VectorH2<Kernel>                        Vector_2;
    typedef DirectionH2<Kernel>                     Direction_2;
    typedef SegmentC2<Kernel>                       Segment_2;
    typedef LineH2<Kernel>                          Line_2;
    typedef RayC2<Kernel>                           Ray_2;
    typedef CircleC2<Kernel>                        Circle_2;
    typedef TriangleC2<Kernel>                      Triangle_2;
    typedef Iso_rectangleH2<Kernel>                 Iso_rectangle_2;
    typedef Aff_transformationH2<Kernel>            Aff_transformation_2;
    typedef Weighted_pointH2<Kernel>                Weighted_point_2;

    typedef PointH3<Kernel>                         Point_3;
    typedef VectorH3<Kernel>                        Vector_3;
    typedef DirectionH3<Kernel>                     Direction_3;
    typedef SegmentC3<Kernel>                       Segment_3;
    typedef PlaneH3<Kernel>                         Plane_3;
    typedef LineC3<Kernel>                          Line_3;
    typedef RayH3<Kernel>                           Ray_3;
    typedef TriangleC3<Kernel>                      Triangle_3;
    typedef TetrahedronC3<Kernel>                   Tetrahedron_3;
    typedef Iso_cuboidH3<Kernel>                    Iso_cuboid_3;
    typedef SphereH3<Kernel>                        Sphere_3;
    typedef CircleC3<Kernel>                        Circle_3;
    typedef Aff_transformationH3<Kernel>            Aff_transformation_3;
    typedef Weighted_pointH3<Kernel>                Weighted_point_3;

    typedef Cartesian_const_iterator_d<typename std::array<RT, 3>::const_iterator> Cartesian_const_iterator_2;
    typedef Cartesian_const_iterator_d<typename std::array<RT, 4>::const_iterator> Cartesian_const_iterator_3;

    typedef FT_                                     Cartesian_coordinate_type;
    typedef const RT_&                              Homogeneous_coordinate_type;
    // Undocumented stuff.
    typedef Data_accessorH2<Kernel>                 Data_accessor_2;
    typedef ConicHPA2<Point_2, Data_accessor_2>     Conic_2;
    // Functors types and access functions.
#define CGAL_Kernel_pred(Y,Z) typedef HomogeneousKernelFunctors::Y<Kernel> Y; \
                              Y Z() const { return Y(); }
#define CGAL_Kernel_cons(Y,Z) CGAL_Kernel_pred(Y,Z)

#include <CGAL/Kernel/interface_macros.h>

};

} //namespace CGAL

#endif // CGAL_HOMOGENEOUS_BASE_H
