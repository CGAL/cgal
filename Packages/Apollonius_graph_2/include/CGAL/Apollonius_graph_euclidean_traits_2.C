// ======================================================================
//
// Copyright (c) 1997 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       : 
// release_date  : 
//
// file          : include/CGAL/Apollonius_graph_euclidean_traits_2.C
// package       : Apollonius_graph_2
// source        : $RCSfile$
// revision      : $Revision$
// revision_date : $Date$
// author(s)     : Menelaos Karavelas <mkaravel@cse.nd.edu>
//
// coordinator   : Mariette Yvinec  <Mariette Yvinec@sophia.inria.fr>
//
// ======================================================================



// class implementation continued
//=================================

CGAL_BEGIN_NAMESPACE

//***********************************************************************
//***********************************************************************
//                            CONSTRUCTIONS
//***********************************************************************
//***********************************************************************


//-----------------------------------------------------------------------
//                        Apollonius vertex
//-----------------------------------------------------------------------
template < class Point, class We >
inline
Point_2< typename Point::R >
ad_circumcenter_2(const Weighted_point< Point,We >& p,
		  const Weighted_point< Point,We >& q,
		  const Weighted_point< Point,We >& r,
		  Cartesian_tag )
{
  typedef typename Kernel_traits<Point>::Kernel Rep;
  typename Rep::FT x,y;
  ad_circumcenterC2(p.x(),p.y(),p.weight(),
		    q.x(),q.y(),q.weight(),
		    r.x(),r.y(),r.weight(),x,y);
  return Point_2< Rep >(x,y);
}

template < class Point, class We >
inline
Point_2< typename Point::R >
ad_circumcenter_2(const Weighted_point< Point, We >& p,
		  const Weighted_point< Point, We >& q,
		  const Weighted_point< Point, We >& r,
		  Homogeneous_tag )
{
  typedef typename Kernel_traits<Point>::Kernel Rep;
  typename Rep::RT x,y,w;
  ad_circumcenterH2(p.hx(),p.hy(),p.hw(),p.weight(),
		    q.hx(),q.hy(),q.hw(),q.weight(),
		    r.hx(),r.hy(),r.hw(),r.weight(),
		    x,y,w);
  return Point_2< Rep >(x,y,w);
}

template < class Point, class We >
inline
Point_2< typename Point::R >
ad_circumcenter_2(const Weighted_point< Point, We >& p,
		  const Weighted_point< Point, We >& q,
		  const Weighted_point< Point, We >& r)
{
  typedef typename Kernel_traits<Point>::Kernel::Rep_tag Tag;
  return ad_circumcenter_2< Point, We >(p, q, r, Tag()); 
}


template < class R >
class Construct_Apollonius_vertex_2
{
public:
  typedef typename R::Point_2                Point;
  typedef typename R::RT                     Weight;
  typedef Weighted_point < Point, Weight >   Weighted_point;

  inline
  Point operator() ( const Weighted_point& p,
		     const Weighted_point& q,
		     const Weighted_point& r) const
  {
    //      CGAL_triangulation_precondition( ! collinear(p, q, r) );
    return ad_circumcenter_2< Point, Weight >(p,q,r);
  }
};


//-----------------------------------------------------------------------
//                     Apollonius weighted point
//-----------------------------------------------------------------------
template < class Point, class We >
inline
Weighted_point< Point, We >
ad_circumcircle_2(const Weighted_point< Point,We >& p,
		  const Weighted_point< Point,We >& q,
		  const Weighted_point< Point,We >& r,
		  Cartesian_tag )
{
  typename Kernel_traits<Point>::Kernel::FT x,y;
  We wt;
  ad_circumcircleC2(p.x(),p.y(),p.weight(),
		    q.x(),q.y(),q.weight(),
		    r.x(),r.y(),r.weight(),x,y,wt);
  return Weighted_point< Point, We >(Point(x,y), wt);
}

template < class Point, class We >
inline
Weighted_point< Point, We >
ad_circumcircle_2(const Weighted_point< Point, We >& p,
		  const Weighted_point< Point, We >& q,
		  const Weighted_point< Point, We >& r,
		  Homogeneous_tag )
{
  typename Kernel_traits<Point>::Kernel::RT x,y,w;
  We wt;
  ad_circumcircleH2(p.hx(),p.hy(),p.hw(),p.weight(),
		    q.hx(),q.hy(),q.hw(),q.weight(),
		    r.hx(),r.hy(),r.hw(),r.weight(),
		    x,y,w,wt);
  return Weighted_point< Point, We >(Point(x,y,w), wt);
}

template < class Point, class We >
inline
Weighted_point< Point, We >
ad_circumcircle_2(const Weighted_point< Point, We >& p,
		  const Weighted_point< Point, We >& q,
		  const Weighted_point< Point, We >& r)
{
  typedef typename Kernel_traits<Point>::Kernel::Rep_tag Tag;
  return ad_circumcircle_2< Point, We >(p, q, r, Tag()); 
}

template < class Point, class We, class Line >
inline
Line
ad_left_bitangent_line_2(const Weighted_point< Point,We >& p,
			 const Weighted_point< Point,We >& q,
			 Cartesian_tag )
{
  typename Kernel_traits<Point>::Kernel::FT a, b, c;
  ad_left_bitangent_lineC2(p.x(),p.y(),p.weight(),
			   q.x(),q.y(),q.weight(),
			   a,b,c);
  return Line(a, b, c);
}

template < class Point, class We, class Line >
inline
Line
ad_left_bitangent_line_2(const Weighted_point< Point, We >& p,
			 const Weighted_point< Point, We >& q,
			 Homogeneous_tag )
{
  typename Kernel_traits<Point>::Kernel::RT a, b, c;
  ad_left_bitangent_lineH2(p.hx(),p.hy(),p.hw(),p.weight(),
			   q.hx(),q.hy(),q.hw(),q.weight(),
			   a, b, c);
  return Line(a, b, c);
}

template < class Point, class We, class Line >
inline
Line
ad_left_bitangent_line_2(const Weighted_point< Point, We >& p,
			 const Weighted_point< Point, We >& q)
{
  typedef typename Kernel_traits<Point>::Kernel::Rep_tag Tag;
  return ad_left_bitangent_line_2< Point, We, Line >(p, q, Tag()); 
}


template < class R >
class Construct_Apollonius_weighted_point_2
{
public:
  typedef typename R::Line_2                 Line;
  typedef typename R::RT                     Weight;
  typedef typename R::Point_2                Point;
  typedef Weighted_point < Point, Weight >   Weighted_point;

  inline Weighted_point operator() ( const Weighted_point& p,
				     const Weighted_point& q,
				     const Weighted_point& r) const
  {
    //      CGAL_triangulation_precondition( ! collinear(p, q, r) );
    return ad_circumcircle_2< Point, Weight >(p,q,r);
  }

  inline Line operator()(const Weighted_point &p,
			 const Weighted_point &q) const
  {
    return ad_left_bitangent_line_2< Point, Weight, Line >(p, q);
  }
};


//-----------------------------------------------------------------------
//                        Apollonius bisector
//-----------------------------------------------------------------------


template< class R >
class Construct_Apollonius_bisector_2
{
public:
  typedef typename R::Point_2               Point;
  typedef typename R::RT                    Weight;
  typedef typename R::Line_2                Line;
  typedef Weighted_point< Point, Weight >   Weighted_point;
  typedef Hyperbola_2<Point, Weight>        Hyperbola;

  inline Object operator() (const Weighted_point& p,
			    const Weighted_point& q) const {
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Line l1(p.point(), q.point());
      Line l = l1.perpendicular(midpoint(p, q));
      return make_object(l);
    }

    Hyperbola h(p, q);
    return make_object(h);
  }
};

//-----------------------------------------------------------------------
//                      Apollonius bisector ray
//-----------------------------------------------------------------------


template<class R>
class Construct_Apollonius_bisector_ray_2
{
public:
  typedef typename R::Point_2               Point;
  typedef typename R::RT                    Weight;
  typedef typename R::Line_2                Line;
  typedef typename R::Ray_2                 Ray;
  typedef Weighted_point< Point, Weight >   Weighted_point;
  typedef Hyperbola_ray_2<Point, Weight>    Hyperbola_ray;
  typedef Sign                              Hyperbola_direction;
  typedef Construct_Apollonius_vertex_2<R>  Apollonius_vertex;

  inline Object
  operator() (const Weighted_point& p,
	      const Weighted_point& q,
	      const Point& r,
	      const Hyperbola_direction& direction) const {
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Line l1(q, p);
      Line l = l1.perpendicular(midpoint(p.point(), q.point()));
      Ray ray(r, l.direction());
      return make_object(ray);
    }
    Hyperbola_ray hr(p, q, r, direction);
    return make_object(hr);
  }

  inline Object
  operator() (const Weighted_point& p,
	      const Weighted_point& q,
	      const Weighted_point& r) const {
    Point c = Apollonius_vertex()(p, q, r);
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Line l1(q, p);
      Line l = l1.perpendicular(midpoint(p.point(), q.point()));
      Ray ray(c, l.direction());
      return make_object(ray);
    }
    Hyperbola_ray hr(p, q, c, NEGATIVE);
    return make_object(hr);
  }
};

//-----------------------------------------------------------------------
//                    Apollonius bisector segment
//-----------------------------------------------------------------------

template<class R>
class Construct_Apollonius_bisector_segment_2
{
public:
  typedef typename R::Point_2                 Point;
  typedef typename R::RT                      Weight;
  typedef typename R::Segment_2               Segment;
  typedef Weighted_point< Point, Weight >     Weighted_point;
  typedef Hyperbola_segment_2<Point, Weight>  Hyperbola_segment;
  typedef Construct_Apollonius_vertex_2<R>    Apollonius_vertex;

  inline Object operator() (const Weighted_point& p,
			    const Weighted_point& q,
			    const Point& r, const Point& s) const {
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Segment seg(r, s);
      return make_object(seg);
    }
    Hyperbola_segment hs(p, q, r, s);
    return make_object(hs);
  }

  inline Object operator() (const Weighted_point& p,
			    const Weighted_point& q,
			    const Weighted_point& r,
			    const Weighted_point& s) const {
    Apollonius_vertex apollonius_vertex;
    Point c_pqr = apollonius_vertex(p,q,r);
    Point c_qps = apollonius_vertex(q,p,s);
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Segment seg(c_pqr, c_qps);
      return make_object(seg);
    }
    Hyperbola_segment hs(p, q, c_pqr, c_qps);
    return make_object(hs);
  }

};

//-----------------------------------------------------------------------
//                    Apollonius primal ray
//-----------------------------------------------------------------------


template<class R>
class Construct_Apollonius_primal_ray_2
{
public:
  typedef typename R::Point_2                        Point;
  typedef typename R::RT                             Weight;
  typedef typename R::RT                             RT;
  typedef typename R::Line_2                         Line;
  typedef typename R::Ray_2                          Ray;
  typedef Weighted_point< Point, Weight >            Weighted_point;
  typedef Construct_Apollonius_weighted_point_2<R>   Apollonius_circle;

  inline Ray operator() (const Weighted_point& p,
			 const Weighted_point& r,
			 const Weighted_point& s) const {
    //
    Apollonius_circle apollonius_circle;
    Line l1 = apollonius_circle(r, p);
    Line l2 = apollonius_circle(p, s);

    RT d1 = CGAL_NTS sqrt( CGAL_NTS square(l1.a()) +
			   CGAL_NTS square(l1.b()) );
    RT d2 = CGAL_NTS sqrt( CGAL_NTS square(l2.a()) +
			   CGAL_NTS square(l2.b()) );
    RT a = l1.a() / d1 - l2.a() / d2;
    RT b = l1.b() / d1 - l2.b() / d2;
    Point c(p.x() + b, p.y() - a);
    return Ray(p, c);
  }
};

//-----------------------------------------------------------------------
//                    Apollonius primal segment
//-----------------------------------------------------------------------

template<class R>
class Construct_Apollonius_primal_segment_2
{
public:
  typedef typename R::Point_2                         Point;
  typedef typename R::RT                              Weight;
  typedef typename R::Line_2                          Line;
  typedef typename R::Segment_2                       Segment;
  typedef Weighted_point< Point, Weight >             Weighted_point;
  typedef Hyperbola_segment_2<Point,Weight>           Hyperbola_segment;
  typedef Parabola_segment_2<Point,Weight,Line>       Parabola_segment;
  typedef Construct_Apollonius_weighted_point_2<R>    Apollonius_circle;

  inline Segment
  operator() (const Weighted_point& p,
	      const Weighted_point& q) const {
    //
    return Segment(p.point(), q.point());
  }

  inline Object
  operator() (const Weighted_point& p,
	      const Weighted_point& q,
	      const Point& r, const Point& s) const {
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Segment seg(r, s);
      return make_object(seg);
    }
    Hyperbola_segment hs(p, q, r, s);
    return make_object(hs);
  }

  inline Object
  operator() (const Weighted_point& p,
	      const Weighted_point& q,
	      const Weighted_point& r,
	      const Weighted_point& s) const {
    Apollonius_circle apollonius_circle;
    Weighted_point c_pqr = apollonius_circle(p, q, r);
    Weighted_point c_qps = apollonius_circle(q, p, s);
    //
    Comparison_result cr = CGAL_NTS compare(c_pqr.weight(), c_qps.weight());
    if ( cr == EQUAL ) {
      Segment seg(p.point(), q.point());
      return make_object(seg);
    }
    Hyperbola_segment hs(c_pqr, c_qps, p.point(), q.point());
    return make_object(hs);
  }

#if 0
  inline Object
  operator() (const Weighted_point& p,
	      const Weighted_point& q,
	      const Weighted_point& c) const {
    //
    Circumcircle_object circumcircle;
    Line l = circumcircle(p, q);
    Parabola_segment ps(c, l, p.point(), q.point());
    return make_object(ps);
  }
#endif

  inline Parabola_segment
  operator() (const Weighted_point& p,
	      const Weighted_point& q,
	      const Weighted_point& r) const {
    //
    Apollonius_circle apollonius_circle;
    Weighted_point c = apollonius_circle(p, q, r);
    Line l = apollonius_circle(q, p);
    return Parabola_segment(c, l, q.point(), p.point());
  }
};


//***********************************************************************
//***********************************************************************
//                              PREDICATES
//***********************************************************************
//***********************************************************************

//-----------------------------------------------------------------------
//                        Compare weight
//-----------------------------------------------------------------------

template < class R >
class Compare_weight_2
{
public:
  typedef typename R::Point_2                Point;
  typedef typename R::RT                     Weight;
  typedef Weighted_point< Point, Weight >    Weighted_point;
  typedef Comparison_result                  result_type;

  inline
  Comparison_result operator()(const Weighted_point& p,
			       const Weighted_point& q) const
  {
    return CGAL_NTS compare(p.weight(), q.weight());
  }
};

//-----------------------------------------------------------------------
//                        Is hidden
//-----------------------------------------------------------------------

template < class Point, class We >
inline
bool
ad_is_hidden_test_2(const Weighted_point< Point, We >& p,
		    const Weighted_point< Point, We >& q,
		    Cartesian_tag, Sqrt_field_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::FT  FT;
  return ad_is_hidden_test_sqrtf_C2(p.x(), p.y(), FT(p.weight()),
				    q.x(), q.y(), FT(q.weight()));
}


template < class Point, class We >
inline
bool
ad_is_hidden_test_2(const Weighted_point< Point, We >& p,
		    const Weighted_point< Point, We >& q,
		    Cartesian_tag, Ring_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::FT  FT;
  return ad_is_hidden_test_ring_C2(p.x(), p.y(), FT(p.weight()),
				    q.x(), q.y(), FT(q.weight()));
}



template < class Point, class We, class Method_tag >
inline
bool
ad_is_hidden_test_2(const Weighted_point< Point, We >& p,
		    const Weighted_point< Point, We >& q,
		    Homogeneous_tag)
{
  typedef typename Kernel_traits<Point>::Kernel::RT  RT;
  Sign s = sign_of_ad_distance2_testH2(p.hx(), p.hy(), p.hw(), 
				       RT(p.weight()),
				       q.hx(), q.hy(), q.hw(),
				       RT(q.weight())
				       );
  if ( s == POSITIVE ) { return false; }
  return (CGAL_NTS compare(p.weight(), q.weight()) != SMALLER);
}

template < class Point, class We, class Method_tag >
inline
bool
ad_is_hidden_test_2(const Weighted_point< Point, We >& p,
		    const Weighted_point< Point, We >& q,
		    Cartesian_tag tag)
{
  return ad_is_hidden_test_2(p, q, tag, Method_tag());
}


template< class R, class Method_tag >
class Is_hidden_2
{
public:
  typedef typename R::Point_2                Point;
  typedef typename R::RT                     Weight;
  typedef Weighted_point < Point, Weight >   Weighted_point;
  typedef bool                               result_type;

  inline bool operator()(const Weighted_point &p,
			 const Weighted_point &q) const
  {
    typedef typename Kernel_traits<Point>::Kernel::Rep_tag Tag;
    return
      ad_is_hidden_test_2<Point,Weight,Method_tag>(p, q, Tag());
  }
};


//-----------------------------------------------------------------------
//                    Oriented side of bisector
//-----------------------------------------------------------------------


template < class Point, class We >
inline
Comparison_result
ad_distances_test_2(const Weighted_point< Point, We >& p1,
		    const Weighted_point< Point, We >& p2,
		    const Point& p, Cartesian_tag, Sqrt_field_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::FT  FT;
  return
    compare_ad_distances_test_sqrtf_C2(p1.x(), p1.y(), FT(p1.weight()),
				       p2.x(), p2.y(), FT(p2.weight()),
				       p.x(),  p.y());
}


template < class Point, class We >
inline
Comparison_result
ad_distances_test_2(const Weighted_point< Point, We >& p1,
		    const Weighted_point< Point, We >& p2,
		    const Point& p, Cartesian_tag, Ring_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::FT  FT;
  return
    compare_ad_distances_test_ring_C2(p1.x(), p1.y(), FT(p1.weight()),
				      p2.x(), p2.y(), FT(p2.weight()),
				      p.x(),  p.y());
}



template < class Point, class We, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const Weighted_point< Point, We >& p1,
		    const Weighted_point< Point, We >& p2,
		    const Point& p, Cartesian_tag tag)
{
  return
    ad_distances_test_2(p1, p2, p, tag, Method_tag());
}



template < class Point, class We, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const Weighted_point< Point, We >& p1,
		    const Weighted_point< Point, We >& p2,
		    const Point& p, Homogeneous_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::RT  RT;
  return compare_ad_distances_testH2(p1.hx(), p1.hy(), p1.hw(),
				     RT(p1.weight()),
				     p2.hx(), p2.hy(), p2.hw(),
				     RT(p2.weight()),
				     p.hx(), p.hy(), p.hw());
}





template < class Point, class We, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const Weighted_point< Point, We >& p1,
		    const Weighted_point< Point, We >& p2,
		    const Point& p)
{
  typedef typename Kernel_traits<Point>::Kernel::Rep_tag Tag;
  return ad_distances_test_2<Point,We,Method_tag>(p1, p2, p, Tag());
}



template< class R, class Method_tag >
class Oriented_side_of_bisector_2
{
public:
  typedef typename R::Point_2                Point;
  typedef typename R::RT                     Weight;
  typedef Weighted_point < Point, Weight >   Weighted_point;
  typedef Oriented_side                      result_type;

  inline Oriented_side operator()(const Weighted_point &p1,
				  const Weighted_point &p2,
				  const Point &p) const
  {
    Comparison_result r =
      ad_distances_test_2<Point,Weight,Method_tag>(p1, p2, p);
    if ( r == EQUAL ) { return ON_ORIENTED_BOUNDARY; }
    return ( r == LARGER ) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
  }
};



//-----------------------------------------------------------------------
//                        Vertex conflict
//-----------------------------------------------------------------------


template < class Point, class We >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag, Sqrt_field_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::FT  FT;
  return ad_incircle_test_sqrtf_C2(p1.x(), p1.y(), FT(p1.weight()),
				   p2.x(), p2.y(), FT(p2.weight()),
				    q.x(),  q.y(), FT( q.weight()));
}


template < class Point, class We >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag, Ring_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::FT  FT;
  return ad_incircle_test_ring_C2(p1.x(), p1.y(), FT(p1.weight()),
				  p2.x(), p2.y(), FT(p2.weight()),
				   q.x(),  q.y(), FT( q.weight()));
}


template < class Point, class We, class Method_tag >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag tag)
{
  return ad_incircle_test_2(p1, p2, q, tag, Method_tag());
}


template < class Point, class We, class Method_tag >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >& q,
		   Homogeneous_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::RT  RT;
  return 
    ad_incircle_testH2(p1.hx(), p1.hy(), p1.hw(), RT(p1.weight()),
		       p2.hx(), p2.hy(), p2.hw(), RT(p2.weight()),
		        q.hx(),  q.hy(),  q.hw(), RT( q.weight()));
}


//-----------------------------------------------------------------------


template < class Point, class We >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >& p3,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag, Sqrt_field_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::FT  FT;
  return ad_incircle_test_sqrtf_C2(p1.x(), p1.y(), FT(p1.weight()),
				   p2.x(), p2.y(), FT(p2.weight()),
				   p3.x(), p3.y(), FT(p3.weight()),
				    q.x(),  q.y(), FT( q.weight()));
}


template < class Point, class We >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >& p3,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag, Ring_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::FT  FT;
  return ad_incircle_test_ring_C2(p1.x(), p1.y(), FT(p1.weight()),
				  p2.x(), p2.y(), FT(p2.weight()),
				  p3.x(), p3.y(), FT(p3.weight()),
				  q.x(),  q.y(), FT( q.weight()));
}



template < class Point, class We, class Method_tag >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >& p3,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag tag)
{
  return ad_incircle_test_2(p1, p2, p3, q, tag, Method_tag());
}


template < class Point, class We, class Method_tag >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >& p3,
		   const Weighted_point< Point, We >& q,
		   Homogeneous_tag )
{
  typedef typename Kernel_traits<Point>::Kernel::RT  RT;
  return 
    ad_incircle_testH2(p1.hx(), p1.hy(), p1.hw(), RT(p1.weight()),
		       p2.hx(), p2.hy(), p2.hw(), RT(p2.weight()),
		       p3.hx(), p3.hy(), p3.hw(), RT(p3.weight()),
		        q.hx(),  q.hy(),  q.hw(), RT( q.weight()));
}


template < class R, class Method_tag >
class Vertex_conflict_2
{
public:
  typedef typename R::Point_2                Point;
  typedef typename R::RT                     Weight;
  typedef Weighted_point< Point, Weight >    Weighted_point;
  typedef Sign                               result_type;

  inline
  Sign operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& q) const
  {
    typedef typename Kernel_traits<Point>::Kernel::Rep_tag Tag;
    return
      ad_incircle_test_2<Point,Weight,Method_tag>(p1, p2, p3, q, Tag());
  }


  inline
  Sign operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& q) const
  {
    typedef typename Kernel_traits<Point>::Kernel::Rep_tag Tag;
    return
      ad_incircle_test_2<Point,Weight,Method_tag>(p1, p2, q, Tag());
  }
 


};

//-----------------------------------------------------------------------
//                    Finite edge interior conflict
//-----------------------------------------------------------------------

template < class Pt, class We>
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return
    ad_finite_edge_test_degenerated_sqrtf_C2(p1.x(), p1.y(),
					     FT(p1.weight()),
					     p2.x(), p2.y(),
					     FT(p2.weight()),
					     q.x(),  q.y(),
					     FT( q.weight()), b);
}

template < class Pt, class We>
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return
    ad_finite_edge_test_degenerated_ring_C2(p1.x(), p1.y(),
					    FT(p1.weight()),
					    p2.x(), p2.y(),
					    FT(p2.weight()),
					    q.x(),  q.y(),
					    FT( q.weight()), b);
}


template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag tag)
{
  return
    ad_finite_edge_test_2(p1, p2, q, b, tag, Method_tag());
}

template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& q,
		      bool b, Homogeneous_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::RT  RT;
  return
    ad_finite_edge_test_degeneratedH2(p1.hx(), p1.hy(),
				      p1.hw(),
				      RT(p1.weight()),
				      p2.hx(), p2.hy(),
				      p2.hw(),
				      RT(p2.weight()),
				      q.hx(),  q.hy(),
				      q.hw(),
				      RT( q.weight()), b);
}

template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& q, bool b)
{
  typedef typename Kernel_traits<Pt>::Kernel::Rep_tag Tag;
  return ad_finite_edge_test_2<Pt,We,Method_tag>
    (p1, p2, q, b, Tag());
}

//-----------------------------------------------------------------------

template < class Pt, class We >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return ad_finite_edge_test_degenerated_sqrtf_C2(p1.x(), p1.y(),
						  FT(p1.weight()),
						  p2.x(), p2.y(),
						  FT(p2.weight()),
						  p3.x(), p3.y(),
						  FT(p3.weight()),
						  q.x(),  q.y(),
						  FT( q.weight()), b);
}

template < class Pt, class We >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return ad_finite_edge_test_degenerated_ring_C2(p1.x(), p1.y(),
						 FT(p1.weight()),
						 p2.x(), p2.y(),
						 FT(p2.weight()),
						 p3.x(), p3.y(),
						 FT(p3.weight()),
						 q.x(),  q.y(),
						 FT( q.weight()), b);
}


template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag tag)
{
  return
    ad_finite_edge_test_2(p1, p2, p3, q, b, tag, Method_tag());
}

template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& q,
		      bool b, Homogeneous_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::RT  RT;
  return
    ad_finite_edge_test_degeneratedH2(p1.hx(), p1.hy(),
				      p1.hw(),
				      RT(p1.weight()),
				      p2.hx(), p2.hy(),
				      p2.hw(),
				      RT(p2.weight()),
				      p3.hx(), p3.hy(),
				      p3.hw(),
				      RT(p3.weight()),
				      q.hx(),  q.hy(),
				      q.hw(),
				      RT( q.weight()), b);
}

template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& q, bool b)
{
  typedef typename Kernel_traits<Pt>::Kernel::Rep_tag Tag;
  return
    ad_finite_edge_test_2<Pt,We,Method_tag>(p1, p2, p3, q, b, Tag());
}

//-----------------------------------------------------------------------


template < class Pt, class We >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& p4,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return
    ad_finite_edge_test_sqrtf_C2(p1.x(), p1.y(), FT(p1.weight()),
				 p2.x(), p2.y(), FT(p2.weight()),
				 p3.x(), p3.y(), FT(p3.weight()),
				 p4.x(), p4.y(), FT(p4.weight()),
				 q.x(),  q.y(), FT(q.weight()), b);
}

template < class Pt, class We >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& p4,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return
    ad_finite_edge_test_ring_C2(p1.x(), p1.y(), FT(p1.weight()),
				p2.x(), p2.y(), FT(p2.weight()),
				p3.x(), p3.y(), FT(p3.weight()),
				p4.x(), p4.y(), FT(p4.weight()),
				q.x(),  q.y(), FT(q.weight()), b);
}


template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& p4,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag tag)
{
  return ad_finite_edge_test_2(p1, p2, p3, p4, q, b,
			       tag, Method_tag());
}

template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& p4,
		      const Weighted_point< Pt, We >& q,
		      bool b, Homogeneous_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::RT  RT;
  return
    aw_Voronoi_diagram_finite_edge_testH2(p1.hx(), p1.hy(), p1.hw(),
					 RT(p1.weight()),
					 p2.hx(), p2.hy(), p2.hw(),
					 RT(p2.weight()),
					 p3.hx(), p3.hy(), p3.hw(),
					 RT(p3.weight()),
					 p4.hx(), p4.hy(), p4.hw(),
					 RT(p4.weight()),
					 q.hx(),  q.hy(), q.hw(),
					 RT( q.weight()), b);
}

template < class Pt, class We, class Method_tag >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& p4,
		      const Weighted_point< Pt, We >& q,
		      bool b)
{
  typedef typename Kernel_traits<Pt>::Kernel::Rep_tag Tag;
  return ad_finite_edge_test_2<Pt,We,Method_tag>
    (p1, p2, p3, p4, q, b, Tag());
}




template < class R, class Method_tag >
class Finite_edge_interior_conflict_2
{
public:
  typedef typename R::Point_2              Point;
  typedef typename R::RT                   Weight;
  typedef Weighted_point< Point, Weight >  Weighted_point;
  typedef bool                             result_type;

  inline
  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& q, bool b) const
  {
    return
      ad_finite_edge_test_2<Point,Weight,Method_tag>(p1, p2, p3, q, b);
  }

  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& q, bool b) const
  {
    return
      ad_finite_edge_test_2<Point,Weight,Method_tag>(p1, p2, q, b);
  }



  inline
  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& p4,
		  const Weighted_point& q,
		  bool b) const
  {
    return ad_finite_edge_test_2<Point,Weight,Method_tag>
      (p1, p2, p3, p4, q, b);
  }
};


//-----------------------------------------------------------------------
//                   Infinite edge interior conflict
//-----------------------------------------------------------------------

template < class Pt, class We >
inline
bool
ad_infinite_edge_test_2(const Weighted_point< Pt, We >& p2,
			const Weighted_point< Pt, We >& p3,
			const Weighted_point< Pt, We >& p4,
			const Weighted_point< Pt, We >& q,
			bool b, Cartesian_tag, Sqrt_field_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return
    ad_infinite_edge_test_sqrtf_C2(p2.x(), p2.y(), FT(p2.weight()),
				   p3.x(), p3.y(), FT(p3.weight()),
				   p4.x(), p4.y(), FT(p4.weight()),
				   q.x(),  q.y(), FT( q.weight()), b);
}


template < class Pt, class We >
inline
bool
ad_infinite_edge_test_2(const Weighted_point< Pt, We >& p2,
			const Weighted_point< Pt, We >& p3,
			const Weighted_point< Pt, We >& p4,
			const Weighted_point< Pt, We >& q,
			bool b, Cartesian_tag, Ring_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return
    ad_infinite_edge_test_ring_C2(p2.x(), p2.y(), FT(p2.weight()),
				  p3.x(), p3.y(), FT(p3.weight()),
				  p4.x(), p4.y(), FT(p4.weight()),
				  q.x(),  q.y(), FT( q.weight()), b);
}


template < class Pt, class We, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const Weighted_point< Pt, We >& p2,
			const Weighted_point< Pt, We >& p3,
			const Weighted_point< Pt, We >& p4,
			const Weighted_point< Pt, We >& q,
			bool b, Cartesian_tag tag)
{
  return
    ad_infinite_edge_test_2(p2, p3, p4, q, b, tag, Method_tag());
}



template < class Pt, class We, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const Weighted_point< Pt, We >& p2,
			const Weighted_point< Pt, We >& p3,
			const Weighted_point< Pt, We >& p4,
			const Weighted_point< Pt, We >& q,
			bool b, Homogeneous_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::RT  RT;
  return
    ad_infinite_edge_testH2(p2.hx(), p2.hy(), p2.hw(),
			    RT(p2.weight()),
			    p3.hx(), p3.hy(), p3.hw(),
			    RT(p3.weight()),
			    p4.hx(), p4.hy(), p4.hw(),
			    RT(p4.weight()),
			    q.hx(),  q.hy(), q.hw(),
			    RT( q.weight()), b);
}

template < class Pt, class We, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const Weighted_point< Pt, We >& p2,
			const Weighted_point< Pt, We >& p3,
			const Weighted_point< Pt, We >& p4,
			const Weighted_point< Pt, We >& q, bool b)
{
  typedef typename Kernel_traits<Pt>::Kernel::Rep_tag Tag;
  return ad_infinite_edge_test_2<Pt,We,Method_tag>
    (p2, p3, p4, q, b, Tag());
}


template < class R, class Method_tag >
class Infinite_edge_interior_conflict_2
{
public:
  typedef typename R::Point_2              Point;
  typedef typename R::RT                   Weight;
  typedef Weighted_point< Point, Weight >  Weighted_point;
  typedef bool                             result_type;

  inline
  bool operator()(const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& p4,
		  const Weighted_point& q, bool b) const
  {
    return ad_infinite_edge_test_2<Point,Weight,Method_tag>
      (p2, p3, p4, q, b);
  }
};


//-----------------------------------------------------------------------
//                          Is degenerate
//-----------------------------------------------------------------------



template < class Pt, class We >
inline
bool
ad_is_degenerate_edge_test_2(const Weighted_point< Pt, We >& p1,
			     const Weighted_point< Pt, We >& p2,
			     const Weighted_point< Pt, We >& p3,
			     const Weighted_point< Pt, We >& p4,
			     Cartesian_tag, Sqrt_field_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return
    ad_is_degenerate_edge_test_sqrtf_C2(p1.x(), p1.y(), FT(p1.weight()),
					p2.x(), p2.y(), FT(p2.weight()),
					p3.x(), p3.y(), FT(p3.weight()),
					p4.x(), p4.y(), FT(p4.weight()));
}

template < class Pt, class We >
inline
bool
ad_is_degenerate_edge_test_2(const Weighted_point< Pt, We >& p1,
			     const Weighted_point< Pt, We >& p2,
			     const Weighted_point< Pt, We >& p3,
			     const Weighted_point< Pt, We >& p4,
			     Cartesian_tag, Ring_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::FT  FT;
  return
    ad_is_degenerate_edge_test_ring_C2(p1.x(), p1.y(), FT(p1.weight()),
				       p2.x(), p2.y(), FT(p2.weight()),
				       p3.x(), p3.y(), FT(p3.weight()),
				       p4.x(), p4.y(), FT(p4.weight()));
}


template < class Pt, class We, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const Weighted_point< Pt, We >& p1,
			     const Weighted_point< Pt, We >& p2,
			     const Weighted_point< Pt, We >& p3,
			     const Weighted_point< Pt, We >& p4,
			     Cartesian_tag tag)
{
  return
    ad_is_degenerate_edge_test_2(p1, p2, p3, p4, tag, Method_tag());
}


template < class Pt, class We, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const Weighted_point< Pt, We >& p1,
			     const Weighted_point< Pt, We >& p2,
			     const Weighted_point< Pt, We >& p3,
			     const Weighted_point< Pt, We >& p4,
			     Homogeneous_tag)
{
  typedef typename Kernel_traits<Pt>::Kernel::RT  RT;
  return
    ad_is_degenerate_edge_testH2(p1.hx(), p1.hy(), p1.hw(),
				 RT(p1.weight()),
				 p2.hx(), p2.hy(), p2.hw(),
				 RT(p2.weight()),
				 p3.hx(), p3.hy(), p3.hw(),
				 RT(p3.weight()),
				 p4.hx(), p4.hy(), p4.hw(),
				 RT(p4.weight()));
}

template < class Pt, class We, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const Weighted_point< Pt, We >& p1,
			     const Weighted_point< Pt, We >& p2,
			     const Weighted_point< Pt, We >& p3,
			     const Weighted_point< Pt, We >& p4)
{
  typedef typename Kernel_traits<Pt>::Kernel::Rep_tag Tag;
  return ad_is_degenerate_edge_test_2<Pt,We,Method_tag>
    (p1, p2, p3, p4, Tag());
}



template < class R, class Method_tag >
class Is_degenerate_edge_2
{
public:
  typedef typename R::Point_2              Point;
  typedef typename R::RT                   Weight;
  typedef Weighted_point< Point, Weight >  Weighted_point;
  typedef bool                             result_type;

  inline
  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& p4) const
  {
    return ad_is_degenerate_edge_test_2<Point,Weight,Method_tag>
      (p1, p2, p3, p4);
  }
};




CGAL_END_NAMESPACE
