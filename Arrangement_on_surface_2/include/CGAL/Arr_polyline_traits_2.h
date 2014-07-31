// Copyright (c) 2000-2014  
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
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
// Author(s)     : Waqar Khan <wkhan@mpi-inf.mpg.de>
#ifndef CGAL_ARR_POLYLINE_TRAITS_2_H
#define CGAL_ARR_POLYLINE_TRAITS_2_H

/*! \file
 * The traits-class for the linear piece-wiese(polyline) type of curves of the
 * arrangement package.
 */

#include <iterator>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

#include <CGAL/basic.h>
#include <CGAL/tags.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_geometry_traits/Polyline_2.h>
#include <CGAL/Arr_polycurve_traits_2.h>
#include <CGAL/Arr_tags.h>
#include <CGAL/Arr_enums.h>

namespace CGAL{

// If no template instantiation is provided, it will be instantiated 
// with Arr_segment_traits.	
template < typename GeometryTraits_2 = Arr_segment_traits_2<> >

class Arr_polyline_traits_2 : public Arr_polycurve_traits_2<GeometryTraits_2>
{
public:	
	typedef GeometryTraits_2 								Geometry_traits_2;

	// Tag definitions:
	typedef Tag_true                                 		Has_left_category;
	typedef Tag_true                                 		Has_merge_category;
	typedef Tag_false                                		Has_do_intersect_category;

	typedef typename Geometry_traits_2::Left_side_category   Left_side_category;
	typedef typename Geometry_traits_2::Bottom_side_category Bottom_side_category;
	typedef typename Geometry_traits_2::Top_side_category    Top_side_category;
	typedef typename Geometry_traits_2::Right_side_category  Right_side_category;

	typedef typename Geometry_traits_2::Has_construct_x_monotone_curve_from_two_points_category  Has_construct_x_monotone_curve_from_two_points_category;

	typedef typename Arr_are_all_sides_oblivious_tag
	<Left_side_category, Bottom_side_category,
	 Top_side_category, Right_side_category>::result
	Are_all_sides_oblivious_tag;

	private:
	typedef Arr_polyline_traits_2<Geometry_traits_2>  Self;
	typedef Arr_polycurve_traits_2<Geometry_traits_2> Base;

	// Data members:
	const Geometry_traits_2* m_geom_traits;    // The base segment-traits class.
	bool m_own_traits;

	private:
	enum { INVALID_INDEX = 0xffffffff };

	public:
	/*! Default constructor */
	Arr_polyline_traits_2() :
	  m_geom_traits(new Geometry_traits_2()), m_own_traits(true) {}

	/*! Constructor with given segment traits
	 * \param geom_traits an already existing segment tarits which is passed will
	 *        be used by the class.
	 */
	Arr_polyline_traits_2(const Geometry_traits_2* geom_traits) :
	  m_geom_traits(geom_traits), m_own_traits(false){ }

	/* Destructor
	 * Deletes the segment tarits class in case it was constructed during the
	 * construction of this.
	 */
	~Arr_polyline_traits_2(){
	  if (m_own_traits)
	    delete m_geom_traits;
	}

	/*! Obtain the segment traits.
	 * \return the segment traits.
	 */
	const Geometry_traits_2* geometry_traits_2() const { return m_geom_traits; }

	/// \name Types and functors inherited from the base segment traits.
	//@{

	// Traits types:
	typedef typename Geometry_traits_2::Point_2            Point_2;
	typedef typename Geometry_traits_2::X_monotone_curve_2 X_monotone_segment_2;
	typedef typename Geometry_traits_2::Curve_2            Segment_2;

	typedef typename Base::Curve_2 						   Curve_2;
	typedef typename Base::X_monotone_curve_2 			   X_monotone_curve_2;

	typedef typename Geometry_traits_2::Multiplicity       Multiplicity;

// class Compare_x_2 {};
// class Compare_xy_2 {};
// class Number_of_points_2 {};
// class Construct_max_vertex_2 {};
// class Is_vertical_2 {};
// class Compare_y_at_x_2 {};
// class Compare_y_at_x_left_2 {};
// class Compare_y_at_x_right_2 {};
// class Equal_2 {};
// class Compare_endpoints_xy_2 {};
// class Construct_opposite_2 {};
// class Make_x_monotone_2 {};
// class Push_back_2 {}:
// class Push_front_2 {};
// class Split_2 {};
// class Intersect_2 {};
// class Are_mergeable_2 {};
// class Merge_2 {};
// class Construct_curve_2 {};
// class Construct_x_monotone_curve_2 {};
// class Trim_2{};	

};

} // namespace CGAL
#endif