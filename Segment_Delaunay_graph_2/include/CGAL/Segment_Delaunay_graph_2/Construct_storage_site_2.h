// Copyright (c) 2006 Foundation for Research and Technology-Hellas (Greece).
// All rights reserved.
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
// Author(s)     : Menelaos Karavelas <mkaravel@iacm.forth.gr>

#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCT_STORAGE_SITE_2_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCT_STORAGE_SITE_2_H 1

#include <CGAL/license/Segment_Delaunay_graph_2.h>


#include <CGAL/Segment_Delaunay_graph_2/basic.h>


namespace CGAL {

namespace SegmentDelaunayGraph_2 {

template<class STraits>
class Construct_storage_site_2
{
public:
  typedef STraits                                    Storage_traits;
  typedef typename Storage_traits::Storage_site_2    Storage_site_2;
  typedef typename Storage_traits::Point_handle      Point_handle;
  typedef typename Storage_traits::Geom_traits       Geom_traits;

  typedef Storage_site_2                             result_type;

protected:
  typedef typename Geom_traits::Intersections_tag    ITag;

  result_type construct(const Point_handle& h1,
			const Point_handle& h2,
			const Point_handle& h3,
			const Point_handle& h4, const Tag_true&) const {
    return Storage_site_2::construct_storage_site_2(h1, h2, h3, h4);
  }

  inline
  result_type construct(const Point_handle& h1,
			const Point_handle& h2,
			const Point_handle& h3,
			const Point_handle& h4,
			const Point_handle& h5,
			const Point_handle& h6, const Tag_true&) const {
    return Storage_site_2::construct_storage_site_2(h1, h2, h3, h4, h5, h6);
  }

  inline
  result_type construct(const Point_handle& h1,
			const Point_handle& h2,
			const Point_handle& h3,
			const Point_handle& h4,
			bool is_first_exact, const Tag_true&) const {
    return Storage_site_2::construct_storage_site_2(h1, h2, h3, h4,
						    is_first_exact);
  }

  
  result_type construct(const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&, const Tag_false&) const {
    CGAL_error();
    return Storage_site_2();
  }

  inline
  result_type construct(const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&, const Tag_false&) const {
    CGAL_error();
    return Storage_site_2();
  }

  inline
  result_type construct(const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			const Point_handle&,
			bool /* is_first_exact */, const Tag_false&) const {
    CGAL_error();
    return Storage_site_2();
  }

public:
  inline
  result_type operator()(const Point_handle& h) const {
    return Storage_site_2::construct_storage_site_2(h);
  }

  inline
  result_type operator()(const Point_handle& h1,
			 const Point_handle& h2) const {
    return Storage_site_2::construct_storage_site_2(h1, h2);
  }

  inline
  result_type operator()(const Point_handle& h1,
			 const Point_handle& h2,
			 const Point_handle& h3,
			 const Point_handle& h4) const {
    return construct(h1, h2, h3, h4, ITag());
  }

  inline
  result_type operator()(const Point_handle& h1,
			 const Point_handle& h2,
			 const Point_handle& h3,
			 const Point_handle& h4,
			 const Point_handle& h5,
			 const Point_handle& h6) const {
    return construct(h1, h2, h3, h4, h5, h6, ITag());
  }

  inline
  result_type operator()(const Point_handle& h1,
			 const Point_handle& h2,
			 const Point_handle& h3,
			 const Point_handle& h4,
			 bool is_first_exact) const {
    return construct(h1, h2, h3, h4, is_first_exact, ITag());
  }

  // constructs the point of intersection
  inline
  result_type operator()(const Storage_site_2& ss0,
			 const Storage_site_2& ss1) const {
    CGAL_precondition( ss0.is_segment() && ss1.is_segment() );
    return Storage_site_2::construct_storage_site_2
      ( ss0.source_of_supporting_site(),
	ss0.target_of_supporting_site(),
	ss1.source_of_supporting_site(),
	ss1.target_of_supporting_site() );
  }

  Storage_site_2
  split_on_point_first_subsegment(const Storage_site_2& s,
				  const Storage_site_2& p) const
  {
    // Splits the first storage site which is a segment using the
    // second storage site which is an exact point
    // Two new storage sites are created, corresponding to the two
    // subsegments
    CGAL_precondition( s.is_segment() && p.is_point() );

    // computing the first sub-segment
    if ( s.is_input(0) ) {
      if ( p.is_input() ) {
	// both segment and point are input
	return operator()(s.source_of_supporting_site(), p.point());
      } else {
	// segment is input but point is intersection of two segments
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  Geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}
	return operator()( supp0.source_of_supporting_site(),
			   supp0.target_of_supporting_site(),
			   supp1.source_of_supporting_site(),
			   supp1.target_of_supporting_site(),
			   true );
      }
    } else {
      if ( p.is_input() ) {
	// point is input but source of segment is not
	return operator()( s.source_of_supporting_site(),
			   p.point(),
			   s.source_of_crossing_site(0),
			   s.target_of_crossing_site(0),
			   false );
      } else {
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  Geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}

	return operator()( supp0.source_of_supporting_site(),
			   supp0.target_of_supporting_site(),
			   s.source_of_crossing_site(0),
			   s.target_of_crossing_site(0),
			   supp1.source_of_supporting_site(),
			   supp1.target_of_supporting_site() );
      }
    }
  }

  Storage_site_2
  split_on_point_second_subsegment(const Storage_site_2& s,
				   const Storage_site_2& p) const
  {
    CGAL_precondition( s.is_segment() && p.is_point() );

    // computing the second sub-segment
    if ( s.is_input(1) ) {
      if ( p.is_input() ) {
	// both segment and point are input
	return operator()(p.point(), s.target_of_supporting_site());
      } else {
	// segment is input but point is intersection of two segments
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  Geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}
	return operator()( supp0.source_of_supporting_site(),
			   supp0.target_of_supporting_site(),
			   supp1.source_of_supporting_site(),
			   supp1.target_of_supporting_site(),
			   false );
      }
    } else {
      if ( p.is_input() ) {
	// point is input but source of segment is not
	return operator()( p.point(),
			   s.target_of_supporting_site(),
			   s.source_of_crossing_site(1),
			   s.target_of_crossing_site(1),
			   true );
      } else {
	Storage_site_2 supp0 = s.supporting_site();
	Storage_site_2 supp1 = p.supporting_site(0);

	typename Geom_traits::Are_parallel_2 are_parallel =
	  Geom_traits().are_parallel_2_object();

	if ( are_parallel(supp0.site(), supp1.site()) ) {
	  supp1 = p.supporting_site(1);
	}

	return operator()( supp0.source_of_supporting_site(),
			   supp0.target_of_supporting_site(),
			   supp1.source_of_supporting_site(),
			   supp1.target_of_supporting_site(),
			   s.source_of_crossing_site(1),
			   s.target_of_crossing_site(1) );
      }
    }
  }

  Storage_site_2
  construct_first_subsegment(const Storage_site_2& ss0,
			     const Storage_site_2& ss1,
			     const Tag_true&) const
  {
    if ( ss0.is_input(0) ) {
      return Storage_site_2::construct_storage_site_2
	( ss0.source_of_supporting_site(),
	  ss0.target_of_supporting_site(),
	  ss1.source_of_supporting_site(),
	  ss1.target_of_supporting_site(), true );
    } else {
      return Storage_site_2::construct_storage_site_2
	( ss0.source_of_supporting_site(),
	  ss0.target_of_supporting_site(),
	  ss0.source_of_crossing_site(0),
	  ss0.target_of_crossing_site(0),
	  ss1.source_of_supporting_site(),
	  ss1.target_of_supporting_site() );
    }
  }

  Storage_site_2
  construct_second_subsegment(const Storage_site_2& ss0,
			      const Storage_site_2& ss1,
			      const Tag_true&) const
  {
    if ( ss0.is_input(1) ) {
      return Storage_site_2::construct_storage_site_2
	( ss0.source_of_supporting_site(),
	  ss0.target_of_supporting_site(),
	  ss1.source_of_supporting_site(),
	  ss1.target_of_supporting_site(), false );
    } else {
      return Storage_site_2::construct_storage_site_2
	( ss0.source_of_supporting_site(),
	  ss0.target_of_supporting_site(),
	  ss1.source_of_supporting_site(),
	  ss1.target_of_supporting_site(),
	  ss0.source_of_crossing_site(1),
	  ss0.target_of_crossing_site(1) );
    }
  }

  Storage_site_2
  construct_first_subsegment(const Storage_site_2&,
			     const Storage_site_2&,
			     const Tag_false&) const
  {
    CGAL_error();
    return Storage_site_2();
  }

  Storage_site_2
  construct_second_subsegment(const Storage_site_2&,
			      const Storage_site_2&,
			      const Tag_false&) const
  {
    CGAL_error();
    return Storage_site_2();
  }

  // constructs the subsegment with supporting segment ss0 and
  // endpoints the point of intersection of ss1 and ss0; the boolean
  // determines if the first or segment subsegment is constructed
  inline
  result_type operator()(const Storage_site_2& ss0,
			 const Storage_site_2& ss1,
			 bool first) const {
    //    CGAL_precondition( ss0.is_segment() && ss1.is_segment() );
    CGAL_precondition( ss0.is_segment() );
    if ( ss1.is_point() ) {
      if ( first ) {
	return split_on_point_first_subsegment(ss0, ss1);
      } else {
	return split_on_point_second_subsegment(ss0, ss1);
      }
    }

    if ( first ) {
      return construct_first_subsegment(ss0, ss1, ITag());
    } else {
      return construct_second_subsegment(ss0, ss1, ITag());
    }
  }

};

//----------------------------------------------------------------------
//----------------------------------------------------------------------


} //namespace SegmentDelaunayGraph_2

} //namespace CGAL

#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_CONSTRUCT_STORAGE_SITE_2_H
