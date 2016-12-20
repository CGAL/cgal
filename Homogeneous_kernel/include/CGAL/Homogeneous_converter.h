// Copyright (c) 2001-2004  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved. 
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Sylvain Pion
//                 Menelaos Karavelas <mkaravel@cse.nd.edu>

#ifndef CGAL_HOMOGENEOUS_CONVERTER_H
#define CGAL_HOMOGENEOUS_CONVERTER_H

// This file contains the definition of a kernel converter, based on
// Homogeneous representation.  It should work between *Homogeneous<A,B>
// and *Homogeneous<C,D>, provided you give an RT converter from A to C,
// and an FT converter from B to D.

#include <CGAL/basic.h>
#include <CGAL/NT_converter.h>
#include <CGAL/Enum_converter.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>

namespace CGAL {

template <class K1, class K2,
          class RT_Converter = NT_converter<typename K1::RT, typename K2::RT>,
          class FT_Converter = NT_converter<typename K1::FT, typename K2::FT> >
class Homogeneous_converter : public Enum_converter
{
private:
    typedef Enum_converter   Base;

public:
    typedef K1            Source_kernel;
    typedef K2            Target_kernel;
    typedef RT_Converter  Ring_number_type_converter;
    typedef FT_Converter  Field_number_type_converter;

    using Base::operator();

    Bbox_2
    operator()(const Bbox_2& b)
    {
        return b;
    }

    Bbox_3
    operator()(const Bbox_3& b)
    {
        return b;
    }

    typename K2::RT
    operator()(const typename K1::RT &a) const
    {
        return c(a);
    }

    typename K2::FT
    operator()(const typename K1::FT &a) const
    {
        return c(a);
    }

    typename K2::Point_2
    operator()(const typename K1::Point_2 &a) const
    {
	return k.construct_point_2_object()(rc(a.hx()), rc(a.hy()),
		rc(a.hw()));
    }

    typename K2::Weighted_point_2
    operator()(const typename K1::Weighted_point_2 &a) const
    {
      return k.construct_weighted_point_2_object()(operator()(a.point()), operator()(a.weight()),
		rc(a.hw()));
    }

    typename K2::Vector_2
    operator()(const typename K1::Vector_2 &a) const
    {
	return k.construct_vector_2_object()(rc(a.hx()), rc(a.hy()),
		rc(a.hw()));
    }

    typename K2::Direction_2
    operator()(const typename K1::Direction_2 &a) const
    {
	return k.construct_direction_2_object()(rc(a.dx()), rc(a.dy()));
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
	return k.construct_line_2_object()(rc(a.a()), rc(a.b()), rc(a.c()));
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
		                             fc(a.squared_radius()),
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
	return k.construct_iso_rectangle_2_object()(operator()((a.min)()),
		                                    operator()((a.max)()), 0);
    }


    typename K2::Point_3
    operator()(const typename K1::Point_3 &a) const
    {
	return k.construct_point_3_object()(rc(a.hx()), rc(a.hy()),
		rc(a.hz()), rc(a.hw()));
    }

    typename K2::Vector_3
    operator()(const typename K1::Vector_3 &a) const
    {
	return k.construct_vector_3_object()(rc(a.hx()), rc(a.hy()),
		rc(a.hz()), rc(a.hw()));
    }

    typename K2::Direction_3
    operator()(const typename K1::Direction_3 &a) const
    {
	return k.construct_direction_3_object()(rc(a.dx()), rc(a.dy()),
		                                rc(a.dz()));
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
		                             fc(a.squared_radius()),
					     a.orientation());
    }

    typename K2::Circle_3
    operator()(const typename K1::Circle_3 &a) const
    {
        return k.construct_circle_3_object()(operator()(a.center()),
                                             fc(a.squared_radius()),
                                             operator()(a.supporting_plane()));
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
	return k.construct_plane_3_object()(rc(a.a()), rc(a.b()), rc(a.c()),
		                            rc(a.d()));
    }

    typename K2::Iso_cuboid_3
    operator()(const typename K1::Iso_cuboid_3 &a) const
    {
	return k.construct_iso_cuboid_3_object()(operator()((a.min)()),
		                                 operator()((a.max)()), 0);
    }

private:
    RT_Converter rc;
    FT_Converter fc;
    K2 k;
};

// Specialization when converting to the same kernel,
// to avoid making copies.
template < class K, class C1, class C2 >
class Homogeneous_converter <K, K, C1, C2>
{
public:
  template < typename T >
  const T& operator()(const T&t) const { return t; }
};

} //namespace CGAL

#endif // CGAL_HOMOGENEOUS_CONVERTER_H
