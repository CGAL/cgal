// ======================================================================
//
// Copyright (c) 2000,2001,2002,2003 The CGAL Consortium
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
// file          : include/CGAL/Cartesian/Cartesian_base.h
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis
//
// ======================================================================

#ifndef CGAL_CARTESIAN_BASE_H
#define CGAL_CARTESIAN_BASE_H

#define CGAL_REP_CLASS_DEFINED

#include <CGAL/basic.h>
#include <CGAL/basic_classes.h>

#include <CGAL/Cartesian/Point_2.h>
#include <CGAL/Cartesian/Vector_2.h>
#include <CGAL/Cartesian/Direction_2.h>
#include <CGAL/Cartesian/Line_2.h>
#include <CGAL/Cartesian/Ray_2.h>
#include <CGAL/Cartesian/Segment_2.h>
#include <CGAL/Cartesian/Triangle_2.h>
#include <CGAL/Cartesian/Circle_2.h>
#include <CGAL/Cartesian/Iso_rectangle_2.h>
#include <CGAL/Cartesian/Aff_transformation_2.h>
#include <CGAL/Cartesian/Data_accessor_2.h>
#include <CGAL/ConicCPA2.h>

#include <CGAL/Cartesian/global_operators_2.h>
#include <CGAL/Cartesian/predicates_on_points_2.h>
#include <CGAL/Cartesian/predicates_on_directions_2.h>
#include <CGAL/Cartesian/predicates_on_lines_2.h>
#include <CGAL/Cartesian/predicates_on_segments_2.h>
#include <CGAL/Cartesian/distance_predicates_2.h>
#include <CGAL/Cartesian/basic_constructions_2.h>

#include <CGAL/Cartesian/Point_3.h>
#include <CGAL/Cartesian/Vector_3.h>
#include <CGAL/Cartesian/Direction_3.h>
#include <CGAL/Cartesian/Line_3.h>
#include <CGAL/Cartesian/Plane_3.h>
#include <CGAL/Cartesian/Ray_3.h>
#include <CGAL/Cartesian/Segment_3.h>
#include <CGAL/Cartesian/Triangle_3.h>
#include <CGAL/Cartesian/Tetrahedron_3.h>
#include <CGAL/Cartesian/Iso_cuboid_3.h>
#include <CGAL/Cartesian/Sphere_3.h>
#include <CGAL/Cartesian/Aff_transformation_3.h>

#include <CGAL/Cartesian/global_operators_3.h>
#include <CGAL/Cartesian/predicates_on_points_3.h>
#include <CGAL/Cartesian/predicates_on_planes_3.h>
#include <CGAL/Cartesian/distance_predicates_3.h>
#include <CGAL/Cartesian/basic_constructions_3.h>

#include <CGAL/representation_tags.h>

CGAL_BEGIN_NAMESPACE

template < typename K_ >
struct Cartesian_base
{
    typedef K_                                          Kernel;

    typedef Cartesian_tag                               Rep_tag;
    typedef Cartesian_tag                               Kernel_tag;

    typedef CGAL::Object                                Object_2;
    typedef CGAL::Object                                Object_3;

    typedef PointC2<Kernel>                             Point_2;
    typedef VectorC2<Kernel>                            Vector_2;
    typedef DirectionC2<Kernel>                         Direction_2;
    typedef SegmentC2<Kernel>                           Segment_2;
    typedef LineC2<Kernel>                              Line_2;
    typedef RayC2<Kernel>                               Ray_2;
    typedef TriangleC2<Kernel>                          Triangle_2;
    typedef CircleC2<Kernel>                            Circle_2;
    typedef Iso_rectangleC2<Kernel>                     Iso_rectangle_2;
    typedef Aff_transformationC2<Kernel>                Aff_transformation_2;
    typedef Data_accessorC2<Kernel>                     Data_accessor_2;
    typedef ConicCPA2<Point_2,Data_accessor_2>          Conic_2;

    typedef PointC3<Kernel>                             Point_3;
    typedef VectorC3<Kernel>                            Vector_3;
    typedef DirectionC3<Kernel>                         Direction_3;
    typedef LineC3<Kernel>                              Line_3;
    typedef PlaneC3<Kernel>                             Plane_3;
    typedef RayC3<Kernel>                               Ray_3;
    typedef SegmentC3<Kernel>                           Segment_3;
    typedef TriangleC3<Kernel>                          Triangle_3;
    typedef TetrahedronC3<Kernel>                       Tetrahedron_3;
    typedef Iso_cuboidC3<Kernel>                        Iso_cuboid_3;
    typedef SphereC3<Kernel>                            Sphere_3;
    typedef Aff_transformationC3<Kernel>                Aff_transformation_3;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_BASE_H
