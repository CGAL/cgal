// Copyright (c) 2003,2004  INRIA Sophia-Antipolis (France) and
// Notre Dame University (U.S.A.).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
// Author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>



#ifndef CGAL_SEGMENT_VORONOI_DIAGRAM_KERNEL_WRAPPER_2_H
#define CGAL_SEGMENT_VORONOI_DIAGRAM_KERNEL_WRAPPER_2_H

#include <CGAL/Segment_Voronoi_diagram_short_names_2.h>

#include <CGAL/Segment_Voronoi_diagram_site_2.h>
#include <CGAL/Segment_Voronoi_diagram_simple_site_2.h>



CGAL_BEGIN_NAMESPACE

namespace CGALi {

  template<class K, class ITag> struct SVD_Which_site;

  // If the ITag is Tag_true we fully support intersections and
  // therefore we need the full-fletched site.
  template<class K>
  struct SVD_Which_site<K,Tag_true>
  {
    typedef K          Kernel;
    typedef Tag_true   Intersections_tag;

    typedef CGAL::Segment_Voronoi_diagram_site_2<K> Site_2;
  };

  // If the ITag is Tag_false we are happy with the simple site.
  template<class K>
  struct SVD_Which_site<K,Tag_false>
  {
    typedef K          Kernel;
    typedef Tag_false  Intersections_tag;

    typedef CGAL::Segment_Voronoi_diagram_simple_site_2<K> Site_2;
  };

} // namespace CGALi



template<class Kernel_base_2, class ITag>
class Segment_Voronoi_diagram_kernel_wrapper_2
  : public Kernel_base_2
{
public:
  typedef Kernel_base_2    Kernel_base;
  typedef ITag             Intersections_tag;

  typedef typename
  CGALi::SVD_Which_site<Kernel_base,Intersections_tag>::Site_2  Site_2;
};




template<class K1, class K2, class Converter>
class Svd_cartesian_converter : public Converter
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
      if ( t.is_exact() ) {
	return K2_Site_2( Base::operator()(t.point()) );
      } else {
	K1_Site_2 s1 = t.supporting_site(0);
	K1_Site_2 s2 = t.supporting_site(1);
	return K2_Site_2( Base::operator()(s1.point(0)),
			  Base::operator()(s1.point(1)),
			  Base::operator()(s2.point(0)),
			  Base::operator()(s2.point(1)) );
      }
    }

    if ( t.is_exact() ) {
      return K2_Site_2( Base::operator()(t.point(0)),
			Base::operator()(t.point(1)) );
    } else {
      K1_Site_2 supp = t.supporting_site();
      if ( t.is_exact(0) ) {
	K1_Site_2 cs = t.crossing_site(1);
	return K2_Site_2(Base::operator()(supp.point(0)),
			 Base::operator()(supp.point(1)),
			 Base::operator()(cs.point(0)),
			 Base::operator()(cs.point(1)), true);
      } else if ( t.is_exact(1) ) {
	K1_Site_2 cs = t.crossing_site(0);
	return K2_Site_2(Base::operator()(supp.point(0)),
			 Base::operator()(supp.point(1)),
			 Base::operator()(cs.point(0)),
			 Base::operator()(cs.point(1)), false);
      } else {
	K1_Site_2 cs1 = t.crossing_site(0);
	K1_Site_2 cs2 = t.crossing_site(1);
	return K2_Site_2(Base::operator()(supp.point(0)),
			 Base::operator()(supp.point(1)),
			 Base::operator()(cs1.point(0)),
			 Base::operator()(cs1.point(1)),
			 Base::operator()(cs2.point(0)),
			 Base::operator()(cs2.point(1)));
      }
    }
  }

  // without intersections
  K2_Site_2 convert_site(const K1_Site_2& t, const Tag_false&) const
  {
    if ( t.is_point() ) {
      return K2_Site_2( Base::operator()(t.point()) );
    }

    // t is a segment
    return K2_Site_2( Base::operator()(t.point(0)),
		      Base::operator()(t.point(1)) );    
  }

public:
  K2_Site_2
  operator()(const K1_Site_2& t) const
  {
    return convert_site(t, intersections_tag());
  }

#if defined(CGAL_CFG_USING_BASE_MEMBER_BUG) || defined(_MSC_VER)
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


CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_KERNEL_WRAPPER_2_H
