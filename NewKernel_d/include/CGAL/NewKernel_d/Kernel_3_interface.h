// Copyright (c) 2014
// INRIA Saclay-Ile de France (France)
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
// Author(s)     : Marc Glisse

#ifndef CGAL_KD_KERNEL_3_INTERFACE_H
#define CGAL_KD_KERNEL_3_INTERFACE_H

#include <CGAL/NewKernel_d/functor_tags.h>
#include <CGAL/transforming_iterator.h>
#include <CGAL/NewKernel_d/utils.h>
#include <CGAL/tuple.h>


namespace CGAL {
template <class Base_> struct Kernel_3_interface : public Base_ {
	typedef Base_ Base;
	typedef Kernel_3_interface<Base> Kernel;
        typedef typename Get_type<Base, RT_tag>::type RT;
        typedef typename Get_type<Base, FT_tag>::type FT;
        typedef typename Get_type<Base, Bool_tag>::type Boolean;
        typedef typename Get_type<Base, Sign_tag>::type Sign;
        typedef typename Get_type<Base, Comparison_result_tag>::type Comparison_result;
        typedef typename Get_type<Base, Orientation_tag>::type Orientation;
        typedef typename Get_type<Base, Oriented_side_tag>::type Oriented_side;
        typedef typename Get_type<Base, Bounded_side_tag>::type Bounded_side;
        typedef typename Get_type<Base, Angle_tag>::type Angle;
	typedef typename Get_type<Base, Point_tag>::type	Point_3;
	typedef typename Get_type<Base, Vector_tag>::type	Vector_3;
	typedef typename Get_type<Base, Segment_tag>::type	Segment_3;
	typedef cpp0x::tuple<Point_3,Point_3,Point_3>		Triangle_3; // placeholder
	typedef cpp0x::tuple<Point_3,Point_3,Point_3,Point_3>	Tetrahedron_3; // placeholder
	struct Compare_xyz_3 {
		typedef typename Get_functor<Base, Compare_lexicographically_tag>::type CL;
		typedef typename CL::result_type result_type;
		CL cl;
		Compare_xyz_3(Kernel const&k):cl(k){}
		result_type operator()(Point_3 const&a, Point_3 const&b) {
			return cl(a,b);
		}
	};
	struct Compare_distance_3 {
		typedef typename Get_functor<Base, Compare_distance_tag>::type CD;
		typedef typename CD::result_type result_type;
		CD cd;
		Compare_distance_3(Kernel const&k):cd(k){}
		result_type operator()(Point_3 const&a, Point_3 const&b, Point_3 const&c) {
			return cd(a,b,c);
		}
		result_type operator()(Point_3 const&a, Point_3 const&b, Point_3 const&c, Point_3 const&d) {
			return cd(a,b,c,d);
		}
	};
	struct Orientation_3 {
		typedef typename Get_functor<Base, Orientation_of_points_tag>::type O;
		typedef typename O::result_type result_type;
		O o;
		Orientation_3(Kernel const&k):o(k){}
		result_type operator()(Point_3 const&a, Point_3 const&b, Point_3 const&c, Point_3 const&d) {
			//return o(a,b,c,d);
			Point_3 const* t[4]={&a,&b,&c,&d};
			return o(make_transforming_iterator<Dereference_functor>(t+0),make_transforming_iterator<Dereference_functor>(t+4));

		}
	};
	struct Side_of_oriented_sphere_3 {
		typedef typename Get_functor<Base, Side_of_oriented_sphere_tag>::type SOS;
		typedef typename SOS::result_type result_type;
		SOS sos;
		Side_of_oriented_sphere_3(Kernel const&k):sos(k){}
		result_type operator()(Point_3 const&a, Point_3 const&b, Point_3 const&c, Point_3 const&d, Point_3 const&e) {
			//return sos(a,b,c,d);
			Point_3 const* t[5]={&a,&b,&c,&d,&e};
			return sos(make_transforming_iterator<Dereference_functor>(t+0),make_transforming_iterator<Dereference_functor>(t+5));
		}
	};

	// I don't have the Coplanar predicates (yet)


	Compare_xyz_3 compare_xyz_3_object()const{ return Compare_xyz_3(*this); }
	Compare_distance_3 compare_distance_3_object()const{ return Compare_distance_3(*this); }
	Orientation_3 orientation_3_object()const{ return Orientation_3(*this); }
	Side_of_oriented_sphere_3 side_of_oriented_sphere_3_object()const{ return Side_of_oriented_sphere_3(*this); }
};
}

#endif
