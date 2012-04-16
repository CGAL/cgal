// Copyright (c) 2003,2004,2005,2006  INRIA Sophia-Antipolis (France).
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



#ifndef CGAL_SEGMENT_DELAUNAY_GRAPH_2_CARTESIAN_CONVERTER_H
#define CGAL_SEGMENT_DELAUNAY_GRAPH_2_CARTESIAN_CONVERTER_H

#include <CGAL/Segment_Delaunay_graph_2/basic.h>


namespace CGAL {

namespace SegmentDelaunayGraph_2 {


template<class K1, class K2, class Converter>
class Cartesian_converter : public Converter
{
private:
  typedef typename K1::Site_2     K1_Site_2;
  typedef typename K1::Point_2    K1_Point_2;

  typedef typename K2::Site_2     K2_Site_2;
  typedef typename K2::Point_2    K2_Point_2;

  typedef Converter               Base;

  typedef typename K1::Intersections_tag  Intersections_tag;

private:
  static const Intersections_tag&  intersections_tag()
  {
    static Intersections_tag itag;
    return itag;
  }

private:
  // with intersections
  K2_Site_2 convert_site(const K1_Site_2& t, const Tag_true&) const
  {
    if ( t.is_point() ) {
      if ( t.is_input() ) {
	return K2_Site_2::construct_site_2( Base::operator()(t.point()) );
      } else {
	return K2_Site_2::construct_site_2
	  ( Base::operator()(t.source_of_supporting_site(0)),
	    Base::operator()(t.target_of_supporting_site(0)),
	    Base::operator()(t.source_of_supporting_site(1)),
	    Base::operator()(t.target_of_supporting_site(1)) );
      }
    }

    if ( t.is_input() ) {
      return K2_Site_2::construct_site_2
	( Base::operator()(t.source_of_supporting_site()),
	  Base::operator()(t.target_of_supporting_site()) );
    } else {
      if ( t.is_input(0) ) {
	return K2_Site_2::construct_site_2
	  ( Base::operator()(t.source_of_supporting_site()),
	    Base::operator()(t.target_of_supporting_site()),
	    Base::operator()(t.source_of_crossing_site(1)),
	    Base::operator()(t.target_of_crossing_site(1)),
	    true );
      } else if ( t.is_input(1) ) {
	return K2_Site_2::construct_site_2
	  ( Base::operator()(t.source_of_supporting_site()),
	    Base::operator()(t.target_of_supporting_site()),
	    Base::operator()(t.source_of_crossing_site(0)),
	    Base::operator()(t.target_of_crossing_site(0)),
	    false );
      } else {
	return K2_Site_2::construct_site_2
	  ( Base::operator()(t.source_of_supporting_site()),
	    Base::operator()(t.target_of_supporting_site()),
	    Base::operator()(t.source_of_crossing_site(0)),
	    Base::operator()(t.target_of_crossing_site(0)),
	    Base::operator()(t.source_of_crossing_site(1)),
	    Base::operator()(t.target_of_crossing_site(1)) );
      }
    }
  }

  // without intersections
  K2_Site_2 convert_site(const K1_Site_2& t, const Tag_false&) const
  {
    if ( t.is_point() ) {
      return K2_Site_2::construct_site_2( Base::operator()(t.point()) );
    }

    // t is a segment
    return K2_Site_2::construct_site_2
      ( Base::operator()(t.source_of_supporting_site()),
	Base::operator()(t.target_of_supporting_site()) );    
  }

public:
  K2_Site_2
  operator()(const K1_Site_2& t) const
  {
    return convert_site(t, intersections_tag());
  }

#if defined(_MSC_VER)
  K2_Point_2
  operator()(const K1_Point_2& p) const
  {
    return Base::operator()(p);
  }

  Sign
  operator()(const Sign& s) const
  {
    return s;
  }
#else
  using Base::operator();
#endif
};


} //namespace SegmentDelaunayGraph_2

} //namespace CGAL


#endif // CGAL_SEGMENT_DELAUNAY_GRAPH_2_CARTESIAN_CONVERTER_H
