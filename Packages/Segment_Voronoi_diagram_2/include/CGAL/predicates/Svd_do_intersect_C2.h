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




#ifndef CGAL_SVD_DO_INTERSECT_C2_H
#define CGAL_SVD_DO_INTERSECT_C2_H

#include <CGAL/enum.h>
#include <CGAL/determinant.h>
#include <CGAL/predicates/Svd_basic_predicates_C2.h>
#include <CGAL/predicates/Svd_are_same_points_C2.h>

CGAL_BEGIN_NAMESPACE

//---------------------------------------------------------------------

template<class RT>
std::pair<int,int>
svd_do_intersect_C2(const RT& x1, const RT& y1,
		    const RT& x2, const RT& y2,
		    const RT& x3, const RT& y3,
		    const RT& x4, const RT& y4)
{
  RT delta = -det2x2_by_formula(x2 - x1, x4 - x3, y2 - y1, y4 - y3);

  Sign s = CGAL::sign( delta );
  if ( s != CGAL::ZERO ) {
    return svd_do_intersect_non_parallel_C2(x1, y1, x2, y2,
					    x3, y3, x4, y4, delta);
  } else {
    return svd_do_intersect_parallel_C2(x1, y1, x2, y2,
					x3, y3, x4, y4);
  }
}


//---------------------------------------------------------------------


template<class RT>
std::pair<int,int>
svd_do_intersect_non_parallel_C2(const RT& x1, const RT& y1,
				 const RT& x2, const RT& y2,
				 const RT& x3, const RT& y3,
				 const RT& x4, const RT& y4,
				 const RT& D)
{
  RT Dt = -det2x2_by_formula(x3 - x1, x4 - x3,
			     y3 - y1, y4 - y3);

  RT Ds = det2x2_by_formula(x2 - x1, x3 - x1,
			    y2 - y1, y3 - y1);

  Sign s_D = CGAL::sign( D );
  Sign s_Dt = CGAL::sign( Dt );
  Sign s_Ds = CGAL::sign( Ds );

  Sign s_tdiff = CGAL::sign(Dt - D);
  Sign s_sdiff = CGAL::sign(Ds - D);

  Sign s_t = Sign(s_Dt * s_D);
  Sign s_s = Sign(s_Ds * s_D);

  Sign s_t_minus_1 = Sign(s_tdiff * s_D);
  Sign s_s_minus_1 = Sign(s_sdiff * s_D);

  if ( s_t == CGAL::NEGATIVE || s_t_minus_1 == CGAL::POSITIVE ||
       s_s == CGAL::NEGATIVE || s_s_minus_1 == CGAL::POSITIVE ) {
    //  t < 0 or t > 1 or s < 0 or s > 1
    return std::pair<int,int>(3,3);
  }

  int it(0), is(0);
  if ( s_t == CGAL::ZERO ) {
    it = 0;
  } else if ( s_t_minus_1 == CGAL::ZERO ) {
    it = 1;
  } else {
    it = 2;
  }
  if ( s_s == CGAL::ZERO ) {
    is = 0;
  } else if ( s_s_minus_1 == CGAL::ZERO ) {
    is = 1;
  } else {
    is = 2;
  }
  return std::pair<int,int>(it, is);
}




//---------------------------------------------------------------------

template<class RT>
std::pair<int,int>
svd_do_intersect_parallel_C2(const RT& x1, const RT& y1,
			     const RT& x2, const RT& y2,
			     const RT& x3, const RT& y3,
			     const RT& x4, const RT& y4)
{
  RT D1 = det2x2_by_formula(x2 - x1, x3 - x1,
			    y2 - y1, y3 - y1);

  if ( CGAL::sign( D1 ) != CGAL::ZERO ) {
    return std::pair<int,int>(3,3);
  }

  RT Dt3, Dt4, Dt;
  if ( CGAL::compare(x2, x1) != CGAL::EQUAL ) {
    Dt  = x2 - x1;
    Dt3 = x3 - x1;
    Dt4 = x4 - x1;
  } else {
    Dt  = y2 - y1;
    Dt3 = y3 - y1;
    Dt4 = y4 - y1;
  }

  Sign s_Dt = CGAL::sign( Dt );
  Sign s_Dt3 = CGAL::sign( Dt3 );
  Sign s_Dt4 = CGAL::sign( Dt4 );

  Sign s_t3 = Sign(s_Dt3 * s_Dt);
  Sign s_t4 = Sign(s_Dt4 * s_Dt);

  Sign s_t3diff = CGAL::sign( Dt3 - Dt );
  Sign s_t4diff = CGAL::sign( Dt4 - Dt );

  Sign s_t3_minus_1 = Sign(s_t3diff * s_Dt);
  Sign s_t4_minus_1 = Sign(s_t4diff * s_Dt);

  if ( (s_t3 == CGAL::NEGATIVE || s_t3_minus_1 == CGAL::POSITIVE) &&
       (s_t4 == CGAL::NEGATIVE || s_t4_minus_1 == CGAL::POSITIVE) ) {
    //  (t3 < 0 or t3 > 1) and (t4 < 0 or t4 > 1)
    return std::pair<int,int>(3,3); // no intersection
  }

  int it3(0), it4(0);
  if ( s_t3 == CGAL::ZERO ) { // t3 == 0
    it3 = 0;
  } else if ( s_t3_minus_1 == CGAL::ZERO ) { // t3 == 1
    it3 = 1;
  } else if ( s_t3 == CGAL::POSITIVE &&
	      s_t3_minus_1 == CGAL::NEGATIVE ) { // 0 < t3 < 1
    it3 = 2;
  } else { // t3 < 0 or t3 > 1
    it3 = 3;
  }

  if ( s_t4 == CGAL::ZERO ) { // t4 == 0
    it4 = 0;
  } else if ( s_t4_minus_1 == CGAL::ZERO ) { // t4 == 1
    it4 = 1;
  } else if ( s_t4 == CGAL::POSITIVE &&
	      s_t4_minus_1 == CGAL::NEGATIVE ) { // 0 < t4 < 1
    it4 = 2;
  } else { // t4 < 0 or t4 > 1
    it4 = 3;
  }

  if ( it3 < 2 && it4 < 2 ) { // segments are identical
    return std::pair<int,int>(4,4);
  } else if ( it3 < 2 && it4 == 3 ) { // segments intersect at p1 or p2
    return std::pair<int,int>(it3,0);
  } else if ( it3 == 3 && it4 < 2 ) { // segments intersect at p1 or p2
    return std::pair<int,int>(it4,1);
  } else {
    // MK: this case has to be further investigating to produce finer
    //     answers wrt the exact configuration
    return std::pair<int,int>(5,5);
  }
  
}


//---------------------------------------------------------------------
//---------------------------------------------------------------------

template<class K>
class Svd_do_intersect_C2
  : public Svd_basic_predicates_C2<K>
{
public:
  //  typedef std::pair<Intersection_type,Intersection_type>  result_type;
  typedef std::pair<int,int>       result_type;

  enum
    { FIRST_ENDPOINT = 0, SECOND_ENDPOINT, INTERIOR_POINT,
      NO_INTERSECTION, IDENTICAL, SUBSEGMENT
    } Intersection_type;

private:
  typedef Svd_basic_predicates_C2<K>  Base;

  typedef typename Base::Point_2        Point_2;
  typedef typename Base::Segment_2      Segment_2;
  typedef typename Base::Site_2         Site_2;
  typedef typename Base::Line_2         Line_2;

  typedef typename Base::FT             FT;
  typedef typename Base::RT             RT;

  typedef typename K::Orientation_2     Orientation_2;

  typedef Svd_are_same_points_C2<K>  Are_same_points_C2;

private:
  bool same_segments(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_segment() && q.is_segment() );
    return
      ( are_same(p.source_site(), q.source_site()) &&
	are_same(p.target_site(), q.target_site()) ) ||
      ( are_same(p.source_site(), q.target_site()) &&
	are_same(p.target_site(), q.source_site()) );
  }

  //------------------------------------------------------------------------

  result_type
  do_intersect_same_point(const Site_2& p, const Site_2& q,
			  unsigned int ip, unsigned int iq) const
  {
    CGAL_precondition( ip < 2 && iq < 2 );

    if ( same_segments(p.supporting_site(), q.supporting_site()) ) {
      Line_2 l = compute_supporting_line(p.supporting_segment());
      Line_2 lp;

      if ( ip == 0 ) {
	lp = compute_perpendicular(l, p.segment().source());
      } else {
	lp = compute_perpendicular(l, p.segment().target());
	lp = opposite_line(lp);
      }

      Oriented_side os;

      if ( iq == 0 ) {
	os = oriented_side_of_line(lp, q.segment().target());
      } else {
	os = oriented_side_of_line(lp, q.segment().source());
      }

      CGAL_assertion( os != ON_ORIENTED_BOUNDARY );

      if ( os == ON_POSITIVE_SIDE ) {
	return std::pair<int,int>(ip,iq);
      } else {
	return std::pair<int,int>(5,5);
      }
    }

    Point_2 p1 = p.supporting_segment().source();
    Point_2 p2 = p.supporting_segment().target();
    Point_2 p3;

    if ( iq == 0 ) {
      p3 = q.supporting_segment().target();
    } else {
      p3 = q.supporting_segment().source();
    }

    if ( Orientation_2()(p1, p2, p3) != COLLINEAR ) {
      return std::pair<int,int>(ip,iq);
    } else {
      Segment_2 s1 = p.segment();
      Segment_2 s2 = q.segment();

      return
	svd_do_intersect_parallel_C2( s1.source().x(), s1.source().y(),
				      s1.target().x(), s1.target().y(),
				      s2.source().x(), s2.source().y(),
				      s2.target().x(), s2.target().y() );
    }
  }


public:
  result_type
  operator()(const Site_2& p, const Site_2& q) const
  {
    CGAL_precondition( p.is_segment() && q.is_segment() );

    if ( same_segments(p, q) ) {
      return std::pair<int,int>(4,4);
    }

    if ( are_same(p.source_site(), q.source_site()) ) {
      return do_intersect_same_point(p, q, 0, 0);
    } if ( are_same(p.source_site(), q.target_site()) ) {
      return do_intersect_same_point(p, q, 0, 1);
    } else if ( are_same(p.target_site(), q.source_site()) ) {
      return do_intersect_same_point(p, q, 1, 0);
    } else if ( are_same(p.target_site(), q.target_site()) ) {
      return do_intersect_same_point(p, q, 1, 1);
    }
    

    Segment_2 s1 = p.segment();
    Segment_2 s2 = q.segment();

    std::pair<int,int> res =
      svd_do_intersect_C2( s1.source().x(), s1.source().y(),
			   s1.target().x(), s1.target().y(),
			   s2.source().x(), s2.source().y(),
			   s2.target().x(), s2.target().y() );

    return res;

    //    Intersection_type it1 = res.first;
    //    Intersection_type it2 = res.second;

    //    return result_type(it1, it2);
  }

private:
  Are_same_points_C2  are_same;

};

//---------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_DO_INTERSECT_C2_H
