#ifndef CGAL_SVD_INCIRCLE_2_H
#define CGAL_SVD_INCIRCLE_2_H

#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h>


CGAL_BEGIN_NAMESPACE

//---------------------------------------------------------------------

template<class K, class Method_tag>
class Svd_incircle_2
{
private:
  typedef typename K::Point_2                       Point_2;
  typedef typename K::Segment_2                     Segment_2;
  typedef typename K::Site_2                        Site_2;
  typedef CGAL::Svd_voronoi_vertex_2<K,Method_tag>  Voronoi_vertex_2;

  typedef typename K::FT                            FT;
  typedef typename K::RT                            RT;

private:
  Sign incircle(const Site_2& p, const Site_2& q,
		const Point_2& t) const
  {
    if ( p.is_point() && q.is_point() ) {
      Orientation o = orientation(p.point(), q.point(), t);

      if ( o != COLLINEAR ) {
	return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
      }

      
      RT dtpx = p.point().x() - t.x();
      RT dtpy = p.point().y() - t.y();
      RT dtqx = q.point().x() - t.x();
      RT minus_dtqy = -q.point().y() + t.y();
      
      Sign s = sign_of_determinant2x2(dtpx, dtpy, minus_dtqy, dtqx);

      CGAL_assertion( s != ZERO );

      return s;
    }

    CGAL_assertion( p.is_point() || q.is_point() );

    Orientation o;
    if ( p.is_point() && q.is_segment() ) {
      Point_2 pq = (p.point() == q.source()) ? q.target() : q.source();
      o = orientation(p.point(), pq, t);
    } else { // p is a segment and q is a point
      Point_2 pp = (q.point() == p.source()) ? p.target() : p.source();
      o = orientation(pp, q.point(), t);
    }
    return ( o == RIGHT_TURN ) ? NEGATIVE : POSITIVE;
  }

  //-----------------------------------------------------------------------


  Sign incircle(const Point_2& p, const Point_2& q,
		const Segment_2& t) const
  {    
    if ( (p == t.source() || p == t.target()) &&
	 (q == t.source() || q == t.target()) ) {
	// if t is the segment joining p and q then t must be a vertex
	// on the convex hull
	return NEGATIVE;
    } else if ( p == t.source() || p == t.target() ) {
      // p is an endpoint of t
      // in this case the p,q,oo vertex is destroyed only if the
      // other endpoint of t is beyond
      Point_2 pt = (p == t.source()) ? t.target() : t.source();
      Orientation o = orientation(p, q, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else if ( q == t.source() || q == t.target() ) {
      Point_2 pt = (q == t.source()) ? t.target() : t.source();
      Orientation o = orientation(p, q, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // maybe here I should immediately return POSITIVE;
      // since we insert endpoints of segments first, p and q cannot
      // be consecutive points on the convex hull if one of the
      // endpoints of t is to the right of the line pq.
      Orientation o1 = orientation(p, q, t.source());
      Orientation o2 = orientation(p, q, t.target());

      if ( o1 == RIGHT_TURN || o2 == RIGHT_TURN ) {
	return NEGATIVE;
      }
      return POSITIVE;
    }
  }


  Sign incircle(const Segment_2& p, const Point_2& q,
		const Segment_2& t) const
  {
    if ( q == t.source() && q == t.target() ) {
      Point_2 pp = (q == p.source()) ? p.target() : p.source();
      Point_2 pt = (q == t.source()) ? t.target() : t.source();

      Orientation o = orientation(pp, q, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      return POSITIVE;
    }
  }


  Sign incircle(const Point_2& p, const Segment_2& q,
		const Segment_2& t) const
  {
    if ( p == t.source() || p == t.target() ) {
      Point_2 pq = (p == q.source()) ? q.target() : q.source();
      Point_2 pt = (p == t.source()) ? t.target() : t.source();

      Orientation o = orientation(p, pq, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // if p is not an endpoint of t, then either p and q should
      // not be on the convex hull or t does not affect the vertex
      // of p and q.
      return POSITIVE;
    }
  }


  Sign incircle(const Site_2& p, const Site_2& q,
		const Segment_2& t) const
  {
    if ( p.is_point() && q.is_point() ) {
      return incircle(p.point(), q.point(), t);
    } else if ( p.is_point() && q.is_segment() ) {
      return incircle(p.point(), q.segment(), t);
    } else { // p is a segment and q is a point
      return incircle(p.segment(), q.point(), t);
    }
  }


public:

  Sign operator()(const Site_2& p, const Site_2& q,
		  const Site_2& r, const Site_2& t) const
  {
    Voronoi_vertex_2 v(p, q, r);

    return v.incircle(t);
  }


  

  Sign operator()(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_assertion( !(p.is_segment() && q.is_segment()) );

    if ( p.is_point() && q.is_segment() ) {
      // p must be an endpoint of q
      CGAL_assertion( (p.point() == q.source()) ||
		      (p.point() == q.target()) );
    } else if ( p.is_segment() && q.is_point() ) {
      // q must be an endpoint of p
      CGAL_assertion( (p.source() == q.point()) ||
		      (p.target() == q.point()) );
    }

    if ( t.is_point() ) {
      return incircle(p, q, t.point());
    }

    return incircle(p, q, t.segment());
  }


};

//---------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_INCIRCLE_2_H
