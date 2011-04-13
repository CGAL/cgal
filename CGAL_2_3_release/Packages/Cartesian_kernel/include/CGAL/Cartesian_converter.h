// ============================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// file          : include/CGAL/Cartesian_converter.h
// revision      : $Revision$
// revision_date : $Date$
// package       : Cartesian_kernel
// author(s)     : Sylvain Pion
// coordinator   : INRIA Sophia-Antipolis (<Mariette.Yvinec@sophia.inria.fr>)
//
// ============================================================================

#ifndef CGAL_CARTESIAN_CONVERTER_H
#define CGAL_CARTESIAN_CONVERTER_H

// This file contains the definition of a kernel converter, based on Cartesian
// representation.  It should work between *Cartesian<A> and *Cartesian<B>,
// provided you give a NT converter from A to B.
// There's a Homogeneous counterpart.

#include <CGAL/basic.h>
#include <CGAL/NT_converter.h>

CGAL_BEGIN_NAMESPACE

template < class K1, class K2,
           class Converter = NT_converter<CGAL_TYPENAME_MSVC_NULL K1::FT,
	                                  CGAL_TYPENAME_MSVC_NULL K2::FT> >
class Cartesian_converter
{
public:

    typename K2::Point_2
    operator()(const typename K1::Point_2 &a) const
    {
	return k.construct_point_2_object()(c(a.x()), c(a.y()));
    }

    typename K2::Vector_2
    operator()(const typename K1::Vector_2 &a) const
    {
	return k.construct_vector_2_object()(c(a.x()), c(a.y()));
    }

    typename K2::Direction_2
    operator()(const typename K1::Direction_2 &a) const
    {
	return k.construct_direction_2_object()(c(a.dx()), c(a.dy()));
    }

    typename K2::Segment_2
    operator()(const typename K1::Segment_2 &a) const
    {
	return k.construct_segment_2_object()(operator()(a.source()),
		                              operator()(a.target()));
    }

    typename K2::Line_2
    operator()(const typename K1::Line_2 &a) const
    {
	return k.construct_line_2_object()(c(a.a()), c(a.b()), c(a.c()));
    }

    typename K2::Ray_2
    operator()(const typename K1::Ray_2 &a) const
    {
	return k.construct_ray_2_object()(operator()(a.source()),
		                          operator()(a.second_point()));
    }

    typename K2::Circle_2
    operator()(const typename K1::Circle_2 &a) const
    {
	return k.construct_circle_2_object()(operator()(a.center()),
		                             c(a.squared_radius()),
					     a.orientation());
    }

    typename K2::Triangle_2
    operator()(const typename K1::Triangle_2 &a) const
    {
	return k.construct_triangle_2_object()(operator()(a.vertex(0)),
		                               operator()(a.vertex(1)),
		                               operator()(a.vertex(2)));
    }

    typename K2::Iso_rectangle_2
    operator()(const typename K1::Iso_rectangle_2 &a) const
    {
	return k.construct_iso_rectangle_2_object()(operator()(a.min()),
		                                    operator()(a.max()));
    }


    typename K2::Point_3
    operator()(const typename K1::Point_3 &a) const
    {
	return k.construct_point_3_object()(c(a.x()), c(a.y()), c(a.z()));
    }

    typename K2::Vector_3
    operator()(const typename K1::Vector_3 &a) const
    {
	return k.construct_vector_3_object()(c(a.x()), c(a.y()), c(a.z()));
    }

    typename K2::Direction_3
    operator()(const typename K1::Direction_3 &a) const
    {
	return k.construct_direction_3_object()(c(a.dx()), c(a.dy()),
		                                c(a.dz()));
    }

    typename K2::Segment_3
    operator()(const typename K1::Segment_3 &a) const
    {
	return k.construct_segment_3_object()(operator()(a.source()),
		                              operator()(a.target()));
    }

    typename K2::Line_3
    operator()(const typename K1::Line_3 &a) const
    {
	return k.construct_line_3_object()(operator()(a.point()),
		                           operator()(a.direction()));
    }

    typename K2::Ray_3
    operator()(const typename K1::Ray_3 &a) const
    {
	return k.construct_ray_3_object()(operator()(a.source()),
		                          operator()(a.second_point()));
    }

    typename K2::Sphere_3
    operator()(const typename K1::Sphere_3 &a) const
    {
	return k.construct_sphere_3_object()(operator()(a.center()),
		                             c(a.squared_radius()),
					     a.orientation());
    }

    typename K2::Triangle_3
    operator()(const typename K1::Triangle_3 &a) const
    {
	return k.construct_triangle_3_object()(operator()(a.vertex(0)),
		                               operator()(a.vertex(1)),
		                               operator()(a.vertex(2)));
    }

    typename K2::Tetrahedron_3
    operator()(const typename K1::Tetrahedron_3 &a) const
    {
	return k.construct_tetrahedron_3_object()(operator()(a.vertex(0)),
		                                  operator()(a.vertex(1)),
		                                  operator()(a.vertex(2)),
		                                  operator()(a.vertex(3)));
    }

    typename K2::Plane_3
    operator()(const typename K1::Plane_3 &a) const
    {
	return k.construct_plane_3_object()(c(a.a()), c(a.b()), c(a.c()),
		                            c(a.d()));
    }

    typename K2::Iso_cuboid_3
    operator()(const typename K1::Iso_cuboid_3 &a) const
    {
	return k.construct_iso_cuboid_3_object()(operator()(a.min()),
		                                 operator()(a.max()));
    }

private:
    Converter c;
    K2 k;
};

CGAL_END_NAMESPACE

#endif // CGAL_CARTESIAN_CONVERTER_H
