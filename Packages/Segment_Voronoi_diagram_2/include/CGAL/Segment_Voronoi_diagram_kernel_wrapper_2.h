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


#include <CGAL/Segment_Voronoi_diagram_site_2.h>



CGAL_BEGIN_NAMESPACE

template<class Kernel_base_2>
class Segment_Voronoi_diagram_kernel_wrapper_2
  : public Kernel_base_2
{
private:
  typedef Segment_Voronoi_diagram_kernel_wrapper_2<Kernel_base_2> Self;

public:
  //  typedef CGAL::Segment_Voronoi_diagram_site_2<Self>  Site_2;
  typedef CGAL::Segment_Voronoi_diagram_site_2<Kernel_base_2>  Site_2;
  //  typedef typename Site_2::Point_2    Point_2;
  //  typedef typename Site_2::Segment_2  Segment_2;

  //  typedef typename Kernel_base_2::Point_2    Point_2;
  //  typedef typename Kernel_base_2::Segment_2  Segment_2;
};




template<class K1, class K2, class Converter >
class Svd_cartesian_converter
  : public Converter
{
private:
  //  typedef Segment_Voronoi_diagram_kernel_wrapper_2<K1> K1W;
  //  typedef Segment_Voronoi_diagram_kernel_wrapper_2<K2> K2W;

  typedef typename K2::Site_2     K2_Site_2;
  typedef typename K2::Point_2    K2_Point_2;
  typedef typename K2::Segment_2  K2_Segment_2;

  typedef Converter               Base;

public:
  K2_Site_2
  operator()(const typename K1::Site_2& t) const
  {
    if ( t.is_point() ) {
      if ( t.is_exact() ) {
	return K2_Site_2( Base::operator()(t.point()) );
      } else {
	K2_Segment_2 s1 = Base::operator()(t.supporting_segment(0));
	K2_Segment_2 s2 = Base::operator()(t.supporting_segment(1));
	return K2_Site_2(s1, s2);
      }
    }

    if ( t.is_exact() ) {
      return K2_Site_2( Base::operator()(t.segment()) );
    } else {
      K2_Segment_2 supp = Base::operator()(t.supporting_segment());
      if ( t.is_exact(0) ) {
	K2_Segment_2 cs = Base::operator()(t.crossing_segment(1));
	return K2_Site_2(supp, cs, true);
      } else if ( t.is_exact(1) ) {
	K2_Segment_2 cs = Base::operator()(t.crossing_segment(0));
	return K2_Site_2(supp, cs, false);
      } else {
	K2_Segment_2 cs1 = Base::operator()(t.crossing_segment(0));
	K2_Segment_2 cs2 = Base::operator()(t.crossing_segment(1));
	return K2_Site_2(supp, cs1, cs2);
      }
    }
  }


  K2_Point_2
  operator()(const typename K1::Point_2& p) const
  {
    return Base::operator()(p);
  }


  K2_Segment_2
  operator()(const typename K1::Segment_2& s) const
  {
    return  Base::operator()(s);
  }


  Sign
  operator()(const Sign& s) const
  {
    return s;
  }

};


CGAL_END_NAMESPACE


#endif // CGAL_SEGMENT_VORONOI_DIAGRAM_KERNEL_WRAPPER_2_H
