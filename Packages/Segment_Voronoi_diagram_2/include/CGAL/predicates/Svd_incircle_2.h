#ifndef CGAL_SVD_INCIRCLE_2_H
#define CGAL_SVD_INCIRCLE_2_H

#include <CGAL/predicates/Segment_Voronoi_diagram_vertex_2.h>
#include <CGAL/predicates/Svd_are_same_points_C2.h>

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

  typedef Svd_are_same_points_C2<K>  Are_same_points_C2;

private:
  Are_same_points_C2  are_same;

private:
  Sign incircle_p(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_precondition( t.is_point() );

    if ( p.is_point() && q.is_point() ) {
      Point_2 pp = p.point(), qp = q.point(), tp = t.point();

      Orientation o = orientation(pp, qp, tp);

      if ( o != COLLINEAR ) {
	return (o == LEFT_TURN) ? POSITIVE : NEGATIVE;
      }

      RT dtpx = pp.x() - tp.x();
      RT dtpy = pp.y() - tp.y();
      RT dtqx = qp.x() - tp.x();
      RT minus_dtqy = -qp.y() + tp.y();
      
      Sign s = sign_of_determinant2x2(dtpx, dtpy, minus_dtqy, dtqx);

      CGAL_assertion( s != ZERO );

      return s;
    }

    CGAL_assertion( p.is_point() || q.is_point() );

    Orientation o;
    if ( p.is_point() && q.is_segment() ) {
      Point_2 pq = are_same(p, q.source_site()) ? q.target() : q.source();
      o = orientation(p.point(), pq, t.point());
    } else { // p is a segment and q is a point
      Point_2 pp = are_same(q, p.source_site()) ? p.target() : p.source();
      o = orientation(pp, q.point(), t.point());
    }
    return ( o == RIGHT_TURN ) ? NEGATIVE : POSITIVE;
  }

  //-----------------------------------------------------------------------


  Sign incircle_pps(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_point() && q.is_point() );
    //    CGAL_precondition( t.is_segment() );

    bool is_p_tsrc = are_same(p, t.source_site());
    bool is_p_ttrg = are_same(p, t.target_site());

    bool is_q_tsrc = are_same(q, t.source_site());
    bool is_q_ttrg = are_same(q, t.target_site());

    bool is_p_on_t = is_p_tsrc || is_p_ttrg;
    bool is_q_on_t = is_q_tsrc || is_q_ttrg;

    if ( is_p_on_t && is_q_on_t ) {
	// if t is the segment joining p and q then t must be a vertex
	// on the convex hull
	return NEGATIVE;
    } else if ( is_p_on_t ) {
      // p is an endpoint of t
      // in this case the p,q,oo vertex is destroyed only if the
      // other endpoint of t is beyond
      Point_2 pt = is_p_tsrc ? t.target() : t.source();
      Orientation o = orientation(p.point(), q.point(), pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else if ( is_q_on_t ) {
      Point_2 pt = is_q_tsrc ? t.target() : t.source();
      Orientation o = orientation(p.point(), q.point(), pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // maybe here I should immediately return POSITIVE;
      // since we insert endpoints of segments first, p and q cannot
      // be consecutive points on the convex hull if one of the
      // endpoints of t is to the right of the line pq.
      Point_2 pp = p.point(), qq = q.point();
      Orientation o1 = orientation(pp, qq, t.source());
      Orientation o2 = orientation(pp, qq, t.target());

      if ( o1 == RIGHT_TURN || o2 == RIGHT_TURN ) {
	return NEGATIVE;
      }
      return POSITIVE;
    }
  }


  Sign incircle_sps(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_segment() && q.is_point() );
    //    CGAL_precondition( t.is_segment() );

    bool is_q_tsrc = are_same(q, t.source_site());
    bool is_q_ttrg = are_same(q, t.target_site());

    bool is_q_on_t = is_q_tsrc || is_q_ttrg;

    //    if ( q == t.source() && q == t.target() ) {
    if ( is_q_on_t ) {
      Point_2 pp = are_same(q, p.source_site()) ? p.target() : p.source();
      Point_2 pt = is_q_tsrc ? t.target() : t.source();

      Orientation o = orientation(pp, q.point(), pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      return POSITIVE;
    }
  }


  Sign incircle_pss(const Site_2& p, const Site_2& q,
		    const Site_2& t) const
  {
    CGAL_precondition( p.is_point() && q.is_segment() );
    //    CGAL_precondition( t.is_segment() );

    bool is_p_tsrc = are_same(p, t.source_site());
    bool is_p_ttrg = are_same(p, t.target_site());

    bool is_p_on_t = is_p_tsrc || is_p_ttrg;

    if ( is_p_on_t ) {
      Point_2 pq = are_same(p, q.source_site()) ? q.target() : q.source();
      Point_2 pt = is_p_tsrc ? t.target() : t.source();

      Orientation o = orientation(p.point(), pq, pt);

      return (o == RIGHT_TURN) ? NEGATIVE : POSITIVE;
    } else {
      // if p is not an endpoint of t, then either p and q should
      // not be on the convex hull or t does not affect the vertex
      // of p and q.
      return POSITIVE;
    }
  }


  Sign incircle_s(const Site_2& p, const Site_2& q,
		  const Site_2& t) const
  {
    CGAL_precondition( t.is_segment() );

    if ( p.is_point() && q.is_point() ) {
      return incircle_pps(p, q, t);
    } else if ( p.is_point() && q.is_segment() ) {
      return incircle_pss(p, q, t);
    } else { // p is a segment and q is a point
      return incircle_sps(p, q, t);
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
      CGAL_assertion( are_same(p, q.source_site()) ||
		      are_same(p, q.target_site()) );
    } else if ( p.is_segment() && q.is_point() ) {
      // q must be an endpoint of p
      CGAL_assertion( are_same(p.source_site(), q) ||
		      are_same(p.target_site(), q) );
    }

    if ( t.is_point() ) {
      return incircle_p(p, q, t);
    }

    return incircle_s(p, q, t.segment());
  }


};

//---------------------------------------------------------------------

CGAL_END_NAMESPACE

#endif // CGAL_SVD_INCIRCLE_2_H
