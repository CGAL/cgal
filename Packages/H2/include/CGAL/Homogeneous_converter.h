// Copyright (c) 2001  Utrecht University (The Netherlands),
// ETH Zurich (Switzerland), Freie Universitaet Berlin (Germany),
// INRIA Sophia-Antipolis (France), Martin-Luther-University Halle-Wittenberg
// (Germany), Max-Planck-Institute Saarbruecken (Germany), RISC Linz (Austria),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Sylvain Pion <Sylvain.Pion@sophia.inria.fr>
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

CGAL_BEGIN_NAMESPACE

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

#ifdef CGAL_CFG_USING_BASE_MEMBER_BUG
    bool operator()(bool b) const { return Base::operator()(b); }
    Sign operator()(Sign s) const { return Base::operator()(s); }

    Oriented_side operator()(Oriented_side os) const {
      return Base::operator()(os);
    }

    Bounded_side operator()(Bounded_side bs) const {
      return Base::operator()(bs);
    }

    Comparison_result operator()(Comparison_result cr) const {
      return Base::operator()(cr);
    }

    Angle operator()(Angle a) const { return Base::operator()(a); }
#else
    using Base::operator();
#endif

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
	return k.construct_iso_rectangle_2_object()(operator()(a.min()),
		                                    operator()(a.max()));
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
	return k.construct_iso_cuboid_3_object()(operator()(a.min()),
		                                 operator()(a.max()));
    }

private:
    RT_Converter rc;
    FT_Converter fc;
    K2 k;
};

CGAL_END_NAMESPACE

#endif // CGAL_HOMOGENEOUS_CONVERTER_H
