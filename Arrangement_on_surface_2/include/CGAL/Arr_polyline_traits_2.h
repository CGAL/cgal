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
template < typename SegmentTraits_2 >

class Arr_polyline_traits_2 : public Arr_polycurve_traits_2<SegmentTraits_2>
{
public:	
	typedef SegmentTraits_2 								Segment_traits_2;

	// Tag definitions:
	typedef Tag_true                                 		Has_left_category;
	typedef Tag_true                                 		Has_merge_category;
	typedef Tag_false                                		Has_do_intersect_category;

	typedef typename Segment_traits_2::Left_side_category   Left_side_category;
	typedef typename Segment_traits_2::Bottom_side_category Bottom_side_category;
	typedef typename Segment_traits_2::Top_side_category    Top_side_category;
	typedef typename Segment_traits_2::Right_side_category  Right_side_category;

	typedef typename Segment_traits_2::Has_construct_x_monotone_curve_from_two_points_category  Has_construct_x_monotone_curve_from_two_points_category;

	typedef typename Arr_are_all_sides_oblivious_tag
	<Left_side_category, Bottom_side_category,
	 Top_side_category, Right_side_category>::result
	Are_all_sides_oblivious_tag;

	private:
	typedef Arr_polyline_traits_2<Segment_traits_2>  Self;
	typedef Arr_polycurve_traits_2<Segment_traits_2> Base;

	// Data members:
	const Segment_traits_2* m_geom_traits;    // The base segment-traits class.
	bool m_own_traits;

	private:
	enum { INVALID_INDEX = 0xffffffff };

	public:
	/*! Default constructor */
	Arr_polyline_traits_2() :
	  m_geom_traits(new Segment_traits_2()), m_own_traits(true) {}

	/*! Constructor with given segment traits
	 * \param geom_traits an already existing segment tarits which is passed will
	 *        be used by the class.
	 */
	Arr_polyline_traits_2(const Segment_traits_2* geom_traits) :
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
	const Segment_traits_2* segment_traits_2() const { return m_geom_traits; }

	/// \name Types and functors inherited from the base segment traits.
	//@{

	// Traits types:
	typedef typename Segment_traits_2::Point_2            Point_2;
	typedef typename Segment_traits_2::X_monotone_curve_2 X_monotone_segment_2;
	typedef typename Segment_traits_2::Curve_2            Segment_2;

	typedef typename Base::Curve_2 						  Curve_2;
	typedef typename Base::X_monotone_curve_2 			  X_monotone_curve_2;

	typedef typename Segment_traits_2::Multiplicity       Multiplicity;



	class Push_back_2 : public Base::Push_back_2{
	protected:
	  typedef Arr_polyline_traits_2<SegmentTraits_2>     Polyline_traits_2;
	  /*! The traits (in case it has state) */
	  const Polyline_traits_2& m_poly_traits;

	public:
	  /*! Constructor. */
	  Push_back_2(const Polyline_traits_2& traits) : Base::Push_back_2(traits), m_poly_traits(traits) {}

	  //http://stackoverflow.com/questions/21168635/inheritance-and-overloading-of-the-function-call-operator
	  using Base::Push_back_2::operator();

	  /* Append a point `p` to an existing polyline `cv` at the back. */
	  void operator()(Curve_2& cv, const Point_2& p) const
	  {
	    typedef typename Curve_2::Segments_size_type size_type;
	    size_type num_seg = cv.number_of_segments();
	    CGAL_precondition(num_seg > 0);
	    int last_seg = num_seg-1;

	    const SegmentTraits_2* seg_traits = m_poly_traits.segment_traits_2();
	    typename SegmentTraits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
	      seg_traits->compare_endpoints_xy_2_object();

	    /*
	     * Since it is desired to maintain `cv` well-oriented, we have
	     * to append the segment [cv[last_seg].target(),p]. The
	     * following test determines which end of the last segment is
	     * the target, i.e. the actual end of `cv`.
	     */
	    if (cmp_seg_endpts(cv[last_seg]) == SMALLER) {
	      typename SegmentTraits_2::Construct_max_vertex_2 get_max_v =
	        seg_traits->construct_max_vertex_2_object();
	      cv.push_back(Segment_2(get_max_v(cv[last_seg]), p));
	    }
	    else {
	      typename SegmentTraits_2::Construct_min_vertex_2 get_min_v =
	        seg_traits->construct_min_vertex_2_object();
	      cv.push_back(Segment_2(get_min_v(cv[last_seg]), p));
	    }
	  }


	  /* Append a segment `seg` to an existing polyline `cv`. If `cv` is
	     empty, `seg` will be its first segment. */
	  // void operator()(Curve_2& cv, const Segment_2& seg) const
	  // { 
	  // 	Base::Push_back_2::operator() (cv, seg);
	  // }

	  /* Append a point `p` to an existing polyline `xcv` at the back. */
	  void operator()(X_monotone_curve_2& xcv, const Point_2& p) const
	  {
	    typedef typename X_monotone_curve_2::Segments_size_type size_type;
	    size_type num_seg = xcv.number_of_segments();
	    CGAL_precondition(num_seg > 0);

	    const SegmentTraits_2* seg_traits = m_poly_traits.segment_traits_2();
	    CGAL_precondition_code
	      (
	       typename SegmentTraits_2::Compare_x_2 comp_x =
	         seg_traits->compare_x_2_object();
	       typename SegmentTraits_2::Compare_xy_2 comp_xy =
	         seg_traits->compare_xy_2_object();
	       typename Base::Is_vertical_2 is_vertical = m_poly_traits.is_vertical_2_object();
	       );

	    if (seg_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
	      // xcv is oriented left-to-right
	      typename SegmentTraits_2::Construct_max_vertex_2 get_max_v =
	        seg_traits->construct_max_vertex_2_object();
	      CGAL_precondition
	        (
	         (!is_vertical(xcv) &&
	          (comp_x(get_max_v(xcv[num_seg-1]), p) == SMALLER)) ||
	         (is_vertical(xcv) &&
	          (comp_x(get_max_v(xcv[num_seg-1]), p) == EQUAL) &&
	          (comp_xy(get_max_v(xcv[num_seg-1]), p) == SMALLER))
	         );
	      xcv.push_back(X_monotone_segment_2(get_max_v(xcv[num_seg-1]), p));
	    }
	    else {
	      // xcv is oriented right-to-left
	      typename SegmentTraits_2::Construct_min_vertex_2 get_min_v =
	        seg_traits->construct_min_vertex_2_object();
	      CGAL_precondition
	        (
	         (!is_vertical(xcv) &&
	          (comp_x(get_min_v(xcv[num_seg-1]), p) == LARGER)) ||
	         (is_vertical(xcv) &&
	          (comp_x(get_min_v(xcv[num_seg-1]), p) == EQUAL) &&
	          (comp_xy(get_min_v(xcv[num_seg-1]), p) == LARGER))
	         );
	      xcv.push_back(X_monotone_segment_2(get_min_v(xcv[num_seg-1]), p));
	    }
	  }

	  // /* Append a segment `seg` to an existing polyline `xcv` at the back. */
	  // void operator()(X_monotone_curve_2& xcv,
	  //                 const X_monotone_segment_2& seg) const
	  // {
	  //   Base::Push_back_2::operator()(xcv,seg);
	  // }
	};

	/*! Get a Push_back_2 functor object. */
	Push_back_2 push_back_2_object() const
	{ return Push_back_2(*this); }


	 /* Functor to augment a polyline by either adding a vertex or a segment
	 * at the front.
	 * TODO: Test all the operator()'s. (Don't forget vertical cases!)
	 */
	class Push_front_2 : public Base::Push_front_2 {
	protected:
	  typedef Arr_polyline_traits_2<SegmentTraits_2>     Polyline_traits_2;
	  /*! The traits (in case it has state) */
	  const Polyline_traits_2& m_poly_traits;

	public:
	  /*! Constructor. */
	  Push_front_2(const Polyline_traits_2& traits) : m_poly_traits(traits), Base::Push_front_2(traits) {}

	  using Base::Push_front_2::operator();

	  /* Append a point `p` to an existing polyline `cv` at the front. */
	  void operator()(Curve_2& cv, const Point_2& p) const
	  {
	    //waqar
	    //CGAL_static_assertion_msg((Has_construct_x_monotone_curve_from_two_points_category()), "X_monotone_curve_2 does not support construction from two points!");

	    CGAL_precondition_code
	      (
	       typedef typename Curve_2::Segments_size_type size_type;
	       size_type num_seg = cv.number_of_segments();
	       );
	    CGAL_precondition(num_seg > 0);

	    const SegmentTraits_2* geom_traits = m_poly_traits.geometry_traits_2();
	    typename SegmentTraits_2::Compare_endpoints_xy_2 cmp_seg_endpts =
	      geom_traits->compare_endpoints_xy_2_object();

	    if (cmp_seg_endpts(cv[0]) == SMALLER) {
	      typename SegmentTraits_2::Construct_min_vertex_2 get_min_v =
	        geom_traits->construct_min_vertex_2_object();
	      cv.push_front(Segment_2(p, get_min_v(cv[0])));
	    }
	    else {
	      typename SegmentTraits_2::Construct_max_vertex_2 get_max_v =
	        geom_traits->construct_max_vertex_2_object();
	      cv.push_front(Segment_2(p, get_max_v(cv[0])));
	    }
	  }


	  /* Append a segment `seg` to an existing polyline `cv` at the front. */
	  void operator()(Curve_2& cv, const Segment_2& seg) const
	  { 
	  	//cv.push_front(seg); 
	  	Base::Push_front_2::operator()(cv, seg);
	  }

	  /* Append a point `p` to an existing polyline `xcv` at the front. */
	  void operator()(const X_monotone_curve_2& xcv, Point_2& p) const
	  {
	    //waqar
	    //CGAL_static_assertion_msg((Has_construct_x_monotone_curve_from_two_points_category()), "X_monotone_curve_2 does not support construction from two points!");

	    const SegmentTraits_2* geom_traits = m_poly_traits.geometry_traits_2();
	    CGAL_precondition_code
	      (
	       typedef typename X_monotone_curve_2::Segments_size_type size_type;
	       size_type num_seg = xcv.number_of_segments();
	       typename SegmentTraits_2::Compare_x_2 comp_x =
	         geom_traits->compare_x_2_object();
	       typename SegmentTraits_2::Compare_xy_2 comp_xy =
	         geom_traits->compare_xy_2_object();
	       typename Base::Is_vertical_2 is_vertical = m_poly_traits.is_vertical_2_object();
	       );
	    CGAL_precondition(num_seg > 0);

	    if (geom_traits->compare_endpoints_xy_2_object()(xcv[0]) == SMALLER) {
	      // xcv is oriented left-to-right
	      typename SegmentTraits_2::Construct_max_vertex_2 get_max_v =
	        geom_traits->construct_max_vertex_2_object();
	      CGAL_precondition
	        (
	         (!is_vertical(xcv) &&
	          (comp_x(get_max_v(xcv[0]), p) == LARGER)) ||
	         (is_vertical(xcv) &&
	          (comp_x(get_max_v(xcv[0]), p) == EQUAL) &&
	          (comp_xy(get_max_v(xcv[0]), p) == LARGER))
	         );
	      xcv.push_front(X_monotone_segment_2(p, get_max_v(xcv[0])));
	    }
	    else {
	      // xcv is oriented right-to-left
	      typename SegmentTraits_2::Construct_min_vertex_2 get_min_v =
	        geom_traits->construct_min_vertex_2_object();
	      CGAL_precondition
	        (
	         (!is_vertical(xcv) &&
	          (comp_x(get_min_v(xcv[0]), p) == SMALLER)) ||
	         (is_vertical(xcv) &&
	          (comp_x(get_min_v(xcv[0]), p) == EQUAL) &&
	          (comp_xy(get_min_v(xcv[0]), p) == SMALLER))
	         );
	      xcv.push_front(X_monotone_segment_2(p, get_min_v(xcv[0])));
	    }
	  }

	  /* Append a segment `seg` to an existing polyline `xcv` at the front. */
	  void operator()(X_monotone_curve_2& xcv,
	                  const X_monotone_segment_2& seg) const
	  {
	   	Base::Push_front_2::operator()(xcv, seg);
	  }
	};

	/*! Get a Push_front_2 functor object. */
	Push_front_2 push_front_2_object() const
	{ return Push_front_2(*this); }

// These functors are defined in the base class.
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