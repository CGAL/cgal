#ifndef CGAL_SVD_ORIENTED_SIDE_OF_BISECTOR_C2_H
#define CGAL_SVD_ORIENTED_SIDE_OF_BISECTOR_C2_H

#include <CGAL/predicates/Svd_are_same_points_C2.h>

CGAL_BEGIN_NAMESPACE

template<class FT, class Method_tag>
Comparison_result
svd_compare_distance_ftC2(const FT& qx, const FT& qy,
			  const FT& sx, const FT& sy,
			  const FT& tx, const FT& ty,
			  const FT& px, const FT& py, Method_tag)
{
  // first check if (qx,qy) is inside, the boundary or at the exterior
  // of the band of the segment s

  //  must_be_filtered(qx);

  FT dx = sx - tx;
  FT dy = sy - ty;


  FT d1x = qx - sx;
  FT d1y = qy - sy;
  FT d2x = qx - tx;
  FT d2y = qy - ty;
#if 0
  std::cout << "inside svd_compare_distance_ftC2" << std::endl;
  std::cout << "px: " << px << std::endl;
  std::cout << "py: " << py << std::endl;
  std::cout << "sx: " << sx << std::endl;
  std::cout << "sy: " << sy << std::endl;
  std::cout << "tx: " << tx << std::endl;
  std::cout << "ty: " << ty << std::endl;
#endif

  if ( px == sx && py == sy ) {
#if 0
    std::cout << "point is same as segment's endpoint 1" << std::endl;
#endif
    FT o = dx * d1x + dy * d1y;
    if ( o >= FT(0) ) {
      return LARGER;
    }
  }

  if ( px == tx && py == ty ) {
#if 0
    std::cout << "point is same as segment's endpoint 2" << std::endl;
#endif
    FT o = dx * d2x + dy * d2y;
    if ( o <= FT(0) ) {
      return LARGER;
    }
  }


  FT d2_from_p = CGAL::square (qx-px) + CGAL::square(qy-py);

  FT dot1 = dx * d1x + dy * d1y;
  if ( dot1 >= FT(0) ) {
    // q is outside (or the boundary of) the band on the side of s.
    FT d2_from_s = CGAL::square(d1x) + CGAL::square(d1y);
    return CGAL::compare(d2_from_s, d2_from_p);
  }

  FT dot2 = dx * d2x + dy * d2y;
  if ( dot2 <= FT(0) ) {
    // q is outside (or the boundary of) the band on the side of t.
    FT d2_from_t = CGAL::square(d2x) + CGAL::square(d2y);
    return CGAL::compare(d2_from_t, d2_from_p);
  }

  // q is strictly in the interior of the band, so I have to compare
  // its distance from the supporting line.
  FT c = sx * ty - sy * tx;
  FT n = CGAL::square(dx) + CGAL::square(dy);
  FT d2_from_l = CGAL::square(qx * dy - qy * dx + c);
  return CGAL::compare(d2_from_l, d2_from_p * n);
}

template<class FT, class Method_tag>
Comparison_result
svd_compare_distance_ftC2(const FT& qx, const FT& qy,
			  const FT& s1x, const FT& s1y,
			  const FT& t1x, const FT& t1y,
			  const FT& s2x, const FT& s2y,
			  const FT& t2x, const FT& t2y, Method_tag)
{
  // first check if (qx,qy) is inside, the boundary or at the exterior
  // of the band of the segments s1, s2

  must_be_filtered(qx);

  FT d1x = s1x - t1x;
  FT d1y = s1y - t1y;

  FT d2x = s2x - t2x;
  FT d2y = s2y - t2y;

  FT dqs1x = qx - s1x;
  FT dqs1y = qy - s1y;
  FT dqt1x = qx - t1x;
  FT dqt1y = qy - t1y;

  FT dqs2x = qx - s2x;
  FT dqs2y = qy - s2y;
  FT dqt2x = qx - t2x;
  FT dqt2y = qy - t2y;

  FT dot_s1 = d1x * dqs1x + d1y * dqs1y;
  int idx1; // 0 for s1, 1 for interior of 1, 2 for t1;

  if ( qx == s1x && qy == s1y ) {
    idx1 = 0;
  } else if ( qx == t1x && qy == t1y ) {
    idx1 = 2;
  } else if ( dot_s1 >= FT(0) ) {
    // q is outside (or the boundary of) the band of 1 on the side of s1.
    idx1 = 0;
  } else {
    FT dot_t1 = d1x * dqt1x + d1y * dqt1y;
    if ( dot_t1 <= FT(0) ) {
      // q is outside (or the boundary of) the band of 1 on the side of t1.
      idx1 = 2;
    } else {
      idx1 = 1;
    }
  }

  FT dot_s2 = d2x * dqs2x + d2y * dqs2y;
  int idx2; // 0 for s2, 1 for interior of 2, 2 for t2;

  if ( qx == s2x && qy == s2y ) {
    idx2 = 0;
  } else if ( qx == t2x && qy == t2y ) {
    idx2 = 2;
  } else if ( dot_s2 >= FT(0) ) {
    // q is outside (or the boundary of) the band of 2 on the side of s2.
    idx2 = 0;
  } else {
    FT dot_t2 = d2x * dqt2x + d2y * dqt2y;
    if ( dot_t2 <= FT(0) ) {
      // q is outside (or the boundary of) the band of 2 on the side of t2.
      idx2 = 2;
    } else {
      idx2 = 1;
    }
  }

  if ( idx1 == 0 ) {
    FT d2_from_s1 = CGAL::square(dqs1x) + CGAL::square(dqs1y);
    //    if ( qx == s1x && qy == s1y ) { d2_from_s1 = FT(0); }
    if ( idx2 == 0 ) {
      FT d2_from_s2 = CGAL::square(dqs2x) + CGAL::square(dqs2y);
      //      if ( qx == s2x && qy == s2y ) { d2_from_s2 = FT(0); }

      if ( s1x == s2x && s1y == s2y ) { return EQUAL; }

      return CGAL::compare(d2_from_s1, d2_from_s2);
    } else if ( idx2 == 2 ) {

      FT d2_from_t2 = CGAL::square(dqt2x) + CGAL::square(dqt2y);
      //      if ( qx == t2x && qy == t2y ) { d2_from_t2 = FT(0); }

      if ( s1x == t2x && s1y == t2y ) { return EQUAL; }

      return CGAL::compare(d2_from_s1, d2_from_t2);
    } else { // idx2 == 1
      FT c2 = s2x * t2y - s2y * t2x;
      FT n2 = CGAL::square(d2x) + CGAL::square(d2y);
      FT d2_from_l2 = CGAL::square(qx * d2y - qy * d2x + c2);

      return CGAL::compare(d2_from_s1 * n2, d2_from_l2);
    }
  } else if ( idx1 == 2 ) {
     FT d2_from_t1 = CGAL::square(dqt1x) + CGAL::square(dqt1y);
     //     if ( qx == t1x && qy == t1y ) { d2_from_t1 = FT(0); }

     if ( idx2 == 0 ) {
       FT d2_from_s2 = CGAL::square(dqs2x) + CGAL::square(dqs2y);
       //       if ( qx == s2x && qy == s2y ) { d2_from_s2 = FT(0); }

       if ( t1x == s2x && t1y == s2y ) { return EQUAL; }

       return CGAL::compare(d2_from_t1, d2_from_s2);
     } else if ( idx2 == 2 ) {
       FT d2_from_t2 = CGAL::square(dqt2x) + CGAL::square(dqt2y);
       //       if ( qx == t2x && qy == t2y ) { d2_from_t2 = FT(0); }

       if ( t1x == t2x && t1y == t2y ) { return EQUAL; }

       return CGAL::compare(d2_from_t1, d2_from_t2);
    } else { // idx2 == 1
      FT c2 = s2x * t2y - s2y * t2x;
      FT n2 = CGAL::square(d2x) + CGAL::square(d2y);
      FT d2_from_l2 = CGAL::square(qx * d2y - qy * d2x + c2);

      return CGAL::compare(d2_from_t1 * n2, d2_from_l2);
    }
  } else { // idx1 == 1
    FT c1 = s1x * t1y - s1y * t1x;
    FT n1 = CGAL::square(d1x) + CGAL::square(d1y);
    FT d2_from_l1 = CGAL::square(qx * d1y - qy * d1x + c1);
    if ( idx2 == 0 ) {
      FT d2_from_s2 = CGAL::square(dqs2x) + CGAL::square(dqs2y);
      //      if ( qx == s2x && qy == s2y ) { d2_from_s2 = FT(0); }

      return CGAL::compare(d2_from_l1, d2_from_s2 * n1);
    } else if ( idx2 == 2 ) {
      FT d2_from_t2 = CGAL::square(dqt2x) + CGAL::square(dqt2y);
      //      if ( qx == t2x && qy == t2y ) { d2_from_t2 = FT(0); }

      return CGAL::compare(d2_from_l1, d2_from_t2 * n1);
    } else { // idx2 == 1
      FT c2 = s2x * t2y - s2y * t2x;
      FT n2 = CGAL::square(d2x) + CGAL::square(d2y);
      FT d2_from_l2 = CGAL::square(qx * d2y - qy * d2x + c2);

      return CGAL::compare(d2_from_l1 * n2, d2_from_l2 * n1);
    }
  }
}


template<class K, class Method_tag>
inline Comparison_result
svd_compare_distance_2(const typename K::Point_2& q,
		       const typename K::Segment_2& s,
		       const typename K::Point_2& p,
		       Cartesian_tag, Method_tag method_tag)
{
  return svd_compare_distance_ftC2(q.x(), q.y(),
				   s[0].x(), s[0].y(),
				   s[1].x(), s[1].y(),
				   p.x(), p.y(), method_tag);
}

template<class K, class Method_tag>
inline Comparison_result
svd_compare_distance_2(const typename K::Point_2& q,
		       const typename K::Segment_2& s,
		       const typename K::Point_2& p,
		       Homogeneous_tag, Method_tag method_tag)
{
  return svd_compare_distanceH2(q.hx(), q.hy(), q.hw(),
				s[0].hx(), s[0].hy(), s[0].hw(),
				s[1].hx(), s[1].hy(), s[1].hw(),
				p.hx(), p.hy(), p.hw(), method_tag);
}


template<class K, class Method_tag>
inline Comparison_result
svd_compare_distance_2(const typename K::Point_2& q,
		       const typename K::Segment_2& s1,
		       const typename K::Segment_2& s2,
		       Cartesian_tag, Method_tag method_tag)
{
  return svd_compare_distance_ftC2(q.x(), q.y(),
				   s1[0].x(), s1[0].y(),
				   s1[1].x(), s1[1].y(),
				   s2[0].x(), s2[0].y(),
				   s2[1].x(), s2[1].y(), method_tag);
}

template<class K, class Method_tag>
inline Comparison_result
svd_compare_distance_2(const typename K::Point_2& q,
		       const typename K::Segment_2& s1,
		       const typename K::Segment_2& s2,
		       Homogeneous_tag, Method_tag method_tag)
{
  return svd_compare_distanceH2(q.hx(), q.hy(), q.hw(),
				s1[0].hx(), s1[0].hy(), s1[0].hw(),
				s1[1].hx(), s1[1].hy(), s1[1].hw(),
				s2[0].hx(), s2[0].hy(), s2[0].hw(),
				s2[1].hx(), s2[1].hy(), s2[1].hw(),
				method_tag);
}



template<class K, class Method_tag>
class Svd_oriented_side_of_bisector_C2
{
public:
  typedef typename K::Site_2     Site_2;
  typedef Oriented_side          result_type;

  typedef typename K::Point_2    Point_2;
  typedef typename K::Segment_2  Segment_2;
  typedef typename K::Rep_tag    Rep_tag;
  typedef Svd_are_same_points_C2<K>  Are_same_points_C2;

private:
  Comparison_result
  operator()(const Point_2& q,
	     const Point_2& p1, const Point_2& p2) const
  {
    CGAL_precondition( p1 != p2 );

    if ( are_same(q, p1) ) { return SMALLER; }
    if ( are_same(q, p2) ) { return LARGER; }
    
    return compare_distance_to_point(q, p1, p2);
  }

  Comparison_result
  operator()(const Point_2& q,
	     const Point_2& p, const Segment_2& s) const
  {
    CGAL_precondition( !s.is_degenerate() );

    return
      opposite( svd_compare_distance_2<K>(q, s, p, Rep_tag(), Method_tag()) );
  }

  Comparison_result
  operator()(const Point_2& q,
	     const Segment_2& s, const Point_2& p) const
  {
    if ( are_same(q, p) ) { return LARGER; }
    if ( are_same(q, s.source()) ) { return SMALLER; }
    if ( are_same(q, s.target()) ) { return SMALLER; }

    return svd_compare_distance_2<K>(q, s, p, Rep_tag(), Method_tag());
  }

  Comparison_result
  operator()(const Point_2& q,
	     const Segment_2& s1, const Segment_2& s2) const
  {
    CGAL_precondition( !s1.is_degenerate() );
    CGAL_precondition( !s2.is_degenerate() );

    if (  ( are_same(q, s1.source()) && are_same(q, s2.source()) ) ||
	  ( are_same(q, s1.source()) && are_same(q, s2.target()) ) ||
	  ( are_same(q, s1.target()) && are_same(q, s2.source()) ) ||
	  ( are_same(q, s1.target()) && are_same(q, s2.target()) )  ) {
      return EQUAL;
    }

    if (  ( are_same(q, s1.source()) || are_same(q, s1.target()) ) &&
	  ( !are_same(q, s2.source()) && !are_same(q, s2.target()) )  ) {
      return SMALLER;
    }

    if (  ( are_same(q, s2.source()) || are_same(q, s2.target()) ) &&
	  ( !are_same(q, s1.source()) && !are_same(q, s1.target()) )  ) {
      return LARGER;
    }

    if (  ( are_same(s1.source(), s2.source()) &&
	    are_same(s1.target(), s2.target()) ) ||
	  ( are_same(s1.source(), s2.target()) &&
	    are_same(s1.target(), s2.source()) )  ) {
      return EQUAL;
    }

    return svd_compare_distance_2<K>(q, s1, s2, Rep_tag(), Method_tag());
  }

protected:
  Oriented_side
  operator()(const Site_2& t1, const Site_2& t2,
	     const Point_2& q) const
  {
    Comparison_result r;

    if ( t1.is_point() && t2.is_point() ) {
      r = operator()(q, t1.point(), t2.point());
      if ( r == LARGER ) { return ON_NEGATIVE_SIDE; }
      if ( r == SMALLER ) { return ON_POSITIVE_SIDE; }
      return ON_ORIENTED_BOUNDARY;
    } else if ( t1.is_segment() && t2.is_point() ) {
      r = operator()(q, t1.segment(), t2.point());
      if ( r == LARGER ) { return ON_NEGATIVE_SIDE; }
      if ( r == SMALLER ) { return ON_POSITIVE_SIDE; }
      return ON_NEGATIVE_SIDE;
    } else if ( t1.is_point() && t2.is_segment() ) {
      r = operator()(q, t1.point(), t2.segment());
      if ( r == LARGER ) { return ON_NEGATIVE_SIDE; }
      if ( r == SMALLER ) { return ON_POSITIVE_SIDE; }
      return ON_POSITIVE_SIDE;
    } else {
      r = operator()(q, t1.segment(), t2.segment());
      if ( r == LARGER ) { return ON_NEGATIVE_SIDE; }
      if ( r == SMALLER ) { return ON_POSITIVE_SIDE; }
      return ON_ORIENTED_BOUNDARY;
    }
  }

public:
  Oriented_side
  operator()(const Site_2& t1, const Site_2& t2,
	     const Site_2& q) const
  {
#if 1
    std::cout << "inside oriented side of bisector top "
	      << "level operator()" << std::endl;
    std::cout << "t1: " << t1 << std::endl;
    std::cout << "t2: " << t2 << std::endl;
    std::cout << "q: " << q << std::endl;
#endif
    CGAL_precondition( q.is_point() );
    return operator()(t1, t2, q.point());
  }


private:
  Are_same_points_C2  are_same;
};


CGAL_END_NAMESPACE

#ifdef CGAL_ARITHMETIC_FILTER_H
#ifndef CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H
#include <CGAL/Arithmetic_filter/predicates/svd_predicates_ftC2.h>
#endif // CGAL_ARITHMETIC_FILTER_SVD_PREDICATES_FTC2_H
#endif

#endif // CGAL_SVD_ORIENTED_SIDE_OF_BISECTOR_C2_H
