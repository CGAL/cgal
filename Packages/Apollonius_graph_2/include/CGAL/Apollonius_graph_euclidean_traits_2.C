// class implementation continued
//=================================

CGAL_BEGIN_NAMESPACE

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
//                        CONSTRUCTIONS
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------

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
  typedef typename Point::R Rep;
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
  typedef typename Point::R Rep;
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
  typedef typename Point::R::Rep_tag Tag;
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
  Point_2< R > operator() ( const Weighted_point& p,
			    const Weighted_point& q,
			    const Weighted_point& r)
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
  typename Point::FT x,y;
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
  typename Point::RT x,y,w;
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
  typedef typename Point::R::Rep_tag Tag;
  return ad_circumcircle_2< Point, We >(p, q, r, Tag()); 
}

//-----------------------------------------------------------------------
// the functions for the construction of the additively weighted
// Delaunay graph left bitangent
//-----------------------------------------------------------------------
template < class Point, class We, class Line >
inline
Line
ad_left_bitangent_line_2(const Weighted_point< Point,We >& p,
			 const Weighted_point< Point,We >& q,
			 Cartesian_tag )
{
  typename Point::FT a, b, c;
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
  typename Point::RT a, b, c;
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
  typedef typename Point::R::Rep_tag Tag;
  return ad_left_bitangent_line_2< Point, We, Line >(p, q, Tag()); 
}


//======================================================================

template < class Point, class Weight, class Line >
class Construct_Apollonius_weighted_point_2
{
public:
  typedef Weighted_point < Point, Weight >   Weighted_point;

  inline Weighted_point operator() ( const Weighted_point& p,
				     const Weighted_point& q,
				     const Weighted_point& r)
  {
    //      CGAL_triangulation_precondition( ! collinear(p, q, r) );
    return ad_circumcircle_2< Point, Weight >(p,q,r);
  }

  inline Line operator()(const Weighted_point &p,
			 const Weighted_point &q)
  {
    return ad_left_bitangent_line_2< Point, Weight, Line >(p, q);
  }
};



//----------------------------------------------------------------------------

template<class Point, class Weight, class Line, class Ray>
class Construct_Apollonius_primal_ray_2
{
public:
  typedef typename Point::R::RT                           RT;
  typedef Weighted_point< Point, Weight >                 Weighted_point;
  typedef Construct_Apollonius_weighted_point_2<Point,Weight,Line>
                                                     Circumcircle_object;

  inline Object operator() (const Weighted_point& p,
			    const Weighted_point& r,
			    const Weighted_point& s) const {
    //
    Circumcircle_object circumcircle;
    Line l1 = circumcircle(r, p);
    Line l2 = circumcircle(p, s);

    RT d1 = CGAL_NTS sqrt( CGAL_NTS square(l1.a()) +
			   CGAL_NTS square(l1.b()) );
    RT d2 = CGAL_NTS sqrt( CGAL_NTS square(l2.a()) +
			   CGAL_NTS square(l2.b()) );
    RT a = l1.a() / d1 - l2.a() / d2;
    RT b = l1.b() / d1 - l2.b() / d2;
    Point c(p.x() + b, p.y() - a);
    Ray ray(p, c);
    return make_object(ray);
  }
};


template<class Point, class Weight, class Line, class Segment>
class Construct_Apollonius_primal_segment_2
{
public:
  typedef Weighted_point< Point, Weight >                 Weighted_point;
  typedef Hyperbola_segment_2<Point,Weight>               Hyperbola_segment;
  typedef Parabola_segment_2<Point,Weight,Line>           Parabola_segment;
  typedef Construct_Apollonius_weighted_point_2<Point,Weight,Line> 
                                                        Circumcircle_object;

  inline Object
  operator() (const Weighted_point& p,
	      const Weighted_point& q) const {
    //
    Segment s(p.point(), q.point());
    return make_object(s);
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
	      const Weighted_point& c) const {
    //
    Circumcircle_object circumcircle;
    Line l = circumcircle(p, q);
    Parabola_segment ps(c, l, p.point(), q.point());
    return make_object(ps);
  }
};


template<class Point, class Weight, class Line, class Ray>
class Construct_Apollonius_bisector_ray_2
{
public:
  typedef Weighted_point< Point, Weight >   Weighted_point;
  typedef Hyperbola_ray_2<Point, Weight>    Hyperbola_ray;
  typedef Sign                              Hyperbola_direction;

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
};

template<class Point, class Weight, class Segment>
class Construct_Apollonius_bisector_segment_2
{
public:
  typedef Weighted_point< Point, Weight >     Weighted_point;
  typedef Hyperbola_segment_2<Point, Weight>  Hyperbola_segment;

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
};


template< class Point, class Weight, class Line >
class Construct_Apollonius_bisector_2
{
public:
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




//***********************************************************************
//***********************************************************************
//************ THE PREDICATES **************
//***********************************************************************
//***********************************************************************


template < class Point, class We >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag, Naive_tag )
{
  typedef typename Point::FT  FT;
  return ad_incircle_test_naive_C2(p1.x(), p1.y(), FT(p1.weight()),
				   p2.x(), p2.y(), FT(p2.weight()),
				    q.x(),  q.y(), FT( q.weight()));
}


template < class Point, class We >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag, Algebraic1_tag )
{
  typedef typename Point::FT  FT;
  return ad_incircle_test_alg1_C2(p1.x(), p1.y(), FT(p1.weight()),
				  p2.x(), p2.y(), FT(p2.weight()),
				   q.x(),  q.y(), FT( q.weight()));
}

template < class Point, class We >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag, Algebraic2_tag )
{
  typedef typename Point::FT  FT;
  return ad_incircle_test_alg2_C2(p1.x(), p1.y(), FT(p1.weight()),
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
  typedef typename Point::RT  RT;
  return 
    ad_incircle_testH2(p1.hx(), p1.hy(), p1.hw(), RT(p1.weight()),
		       p2.hx(), p2.hy(), p2.hw(), RT(p2.weight()),
		        q.hx(),  q.hy(),  q.hw(), RT( q.weight()));
}



//------------------------------------

template < class Point, class We >
inline
Sign
ad_incircle_test_2(const Weighted_point< Point, We >& p1,
		   const Weighted_point< Point, We >& p2,
		   const Weighted_point< Point, We >& p3,
		   const Weighted_point< Point, We >&  q,
		   Cartesian_tag, Naive_tag )
{
  typedef typename Point::FT  FT;
  return ad_incircle_test_naive_C2(p1.x(), p1.y(), FT(p1.weight()),
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
		   Cartesian_tag, Algebraic1_tag )
{
  typedef typename Point::FT  FT;
  return ad_incircle_test_alg1_C2(p1.x(), p1.y(), FT(p1.weight()),
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
		   Cartesian_tag, Algebraic2_tag )
{
  typedef typename Point::FT  FT;
  return ad_incircle_test_alg2_C2(p1.x(), p1.y(), FT(p1.weight()),
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
  typedef typename Point::RT  RT;
  return 
    ad_incircle_testH2(p1.hx(), p1.hy(), p1.hw(), RT(p1.weight()),
		       p2.hx(), p2.hy(), p2.hw(), RT(p2.weight()),
		       p3.hx(), p3.hy(), p3.hw(), RT(p3.weight()),
		        q.hx(),  q.hy(),  q.hw(), RT( q.weight()));
}



template < class Point, class Weight, class Method_tag >
class Vertex_conflict_2
{
public:
  typedef Weighted_point< Point, Weight >    Weighted_point;

  inline
  Sign operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& q)
  {
    typedef typename Point::R::Rep_tag Tag;
    return
      ad_incircle_test_2<Point,Weight,Method_tag>(p1, p2, p3, q, Tag());
  }


  inline
  Sign operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& q)
  {
    typedef typename Point::R::Rep_tag Tag;
    return
      ad_incircle_test_2<Point,Weight,Method_tag>(p1, p2, q, Tag());
  }
 


};

//-----------------------------------------------------------------------
// the functions for the
// Additively_weighted_Voronoi_diagram_compare_weight_test_2 predicate
//-----------------------------------------------------------------------
template < class Point, class Weight >
class Compare_weight_2
{
public:
  typedef Weighted_point< Point, Weight >    Weighted_point;
  inline
  Comparison_result operator()(const Weighted_point& p,
			       const Weighted_point& q)
  {
    return CGAL_NTS compare(p.weight(), q.weight());
  }
};


//-----------------------------------------------------------------------
// the functions for the
// Additively_weighted_Voronoi_diagram_compare_distances_test_2 predicate
//-----------------------------------------------------------------------


template < class Point, class We >
inline
Comparison_result
ad_distances_test_2(const Weighted_point< Point, We >& p1,
		    const Weighted_point< Point, We >& p2,
		    const Point& p, Cartesian_tag, Naive_tag )
{
  typedef typename Point::FT  FT;
  return
    compare_ad_distances_test_naive_C2(p1.x(), p1.y(), FT(p1.weight()),
				       p2.x(), p2.y(), FT(p2.weight()),
				       p.x(),  p.y());
}


template < class Point, class We >
inline
Comparison_result
ad_distances_test_2(const Weighted_point< Point, We >& p1,
		    const Weighted_point< Point, We >& p2,
		    const Point& p, Cartesian_tag, Algebraic1_tag )
{
  typedef typename Point::FT  FT;
  return
    compare_ad_distances_test_alg1_C2(p1.x(), p1.y(), FT(p1.weight()),
				      p2.x(), p2.y(), FT(p2.weight()),
				      p.x(),  p.y());
}

template < class Point, class We >
inline
Comparison_result
ad_distances_test_2(const Weighted_point< Point, We >& p1,
		    const Weighted_point< Point, We >& p2,
		    const Point& p, Cartesian_tag, Algebraic2_tag )
{
  typedef typename Point::FT  FT;
  return
    compare_ad_distances_test_alg2_C2(p1.x(), p1.y(), FT(p1.weight()),
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
  typedef typename Point::RT  RT;
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
  typedef typename Point::R::Rep_tag Tag;
  return ad_distances_test_2<Point,We,Method_tag>(p1, p2, p, Tag());
}



template< class Point, class Weight, class Method_tag >
class Oriented_side_of_bisector_2
{
public:
  typedef Weighted_point < Point, Weight >   Weighted_point;

  inline Oriented_side operator()(const Weighted_point &p1,
				  const Weighted_point &p2,
				  const Point &p)
  {
    Comparison_result r =
      ad_distances_test_2<Point,Weight,Method_tag>(p1, p2, p);
    if ( r == EQUAL ) { return ON_ORIENTED_BOUNDARY; }
    return ( r == LARGER ) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
  }
};


//-----------------------------------------------------------------------


template < class Point, class We >
inline
bool
ad_is_trivial_test_2(const Weighted_point< Point, We >& p,
		     const Weighted_point< Point, We >& q,
		     Cartesian_tag, Naive_tag )
{
  typedef typename Point::FT  FT;
  return ad_is_trivial_test_naive_C2(p.x(), p.y(), FT(p.weight()),
				     q.x(), q.y(), FT(q.weight()));
}


template < class Point, class We, class Algebraic_tag >
inline
bool
ad_is_trivial_test_2(const Weighted_point< Point, We >& p,
		     const Weighted_point< Point, We >& q,
		     Cartesian_tag, Algebraic_tag )
{
  typedef typename Point::FT  FT;
  return ad_is_trivial_test_alg_C2(p.x(), p.y(), FT(p.weight()),
				   q.x(), q.y(), FT(q.weight()));
}



template < class Point, class We, class Method_tag >
inline
bool
ad_is_trivial_test_2(const Weighted_point< Point, We >& p,
		     const Weighted_point< Point, We >& q,
		     Homogeneous_tag)
{
  typedef typename Point::RT  RT;
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
ad_is_trivial_test_2(const Weighted_point< Point, We >& p,
		     const Weighted_point< Point, We >& q,
		     Cartesian_tag tag)
{
  return ad_is_trivial_test_2(p, q, tag, Method_tag());
}




template< class Point, class Weight, class Method_tag >
class Is_trivial_2
{
public:
  typedef Weighted_point < Point, Weight >   Weighted_point;

  inline bool operator()(const Weighted_point &p,
			 const Weighted_point &q)
  {
    typedef typename Point::R::Rep_tag Tag;
    return
      ad_is_trivial_test_2<Point,Weight,Method_tag>(p, q, Tag());
#if 0
    Sign s = ad_distance2_test_2<Point,Weight,Method_tag>(p, q);
    if ( s == POSITIVE ) { return false; }
    return (CGAL_NTS compare(p.weight(), q.weight()) != SMALLER);
#endif
  }
};




//-----------------------------------------------------------------------
// the functions for the Aw_Delaunay_infinite_edge_test_2 predicate
//-----------------------------------------------------------------------

template < class Pt, class We >
inline
bool
ad_infinite_edge_test_2(const Weighted_point< Pt, We >& p2,
			const Weighted_point< Pt, We >& p3,
			const Weighted_point< Pt, We >& p4,
			const Weighted_point< Pt, We >& q,
			bool b, Cartesian_tag, Naive_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_infinite_edge_test_naive_C2(p2.x(), p2.y(), FT(p2.weight()),
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
			bool b, Cartesian_tag, Algebraic1_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_infinite_edge_test_alg1_C2(p2.x(), p2.y(), FT(p2.weight()),
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
			bool b, Cartesian_tag, Algebraic2_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_infinite_edge_test_alg2_C2(p2.x(), p2.y(), FT(p2.weight()),
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
  typedef typename Pt::RT  RT;
  return
    ad_infinite_edge_testH2(p2.hx(), p2.hy(), p2.hw(),
			    RT(p1.weight()),
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
  typedef typename Pt::R::Rep_tag Tag;
  return ad_infinite_edge_test_2<Pt,We,Method_tag>
    (p2, p3, p4, q, b, Tag());
}

template < class Point, class Weight, class Method_tag >
class Infinite_edge_interior_conflict_2
{
public:
  typedef Weighted_point< Point, Weight >  Weighted_point;
  inline
  bool operator()(const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& p4,
		  const Weighted_point& q, bool b)
  {
    return ad_infinite_edge_test_2<Point,Weight,Method_tag>
      (p2, p3, p4, q, b);
  }
};


//-----------------------------------------------------------------------
// the functions for the Aw_Delaunay_finite_edge_test_denegenerated_2
// predicate
//-----------------------------------------------------------------------

// Below are the functions for the case where the degenerated edge has
// two infinite vertices
template < class Pt, class We>
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Naive_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_finite_edge_test_degenerated_naive_C2(p1.x(), p1.y(),
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
		      bool b, Cartesian_tag, Algebraic1_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_finite_edge_test_degenerated_alg1_C2(p1.x(), p1.y(),
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
		      bool b, Cartesian_tag, Algebraic2_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_finite_edge_test_degenerated_alg2_C2(p1.x(), p1.y(),
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
  typedef typename Pt::RT  RT;
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
  typedef typename Pt::R::Rep_tag Tag;
  return ad_finite_edge_test_2<Pt,We,Method_tag>
    (p1, p2, q, b, Tag());
}

// Below are the functions for the case where the degenerated edge has
// one finite vertex and one infinite
template < class Pt, class We >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Naive_tag)
{
  typedef typename Pt::FT  FT;
  return ad_finite_edge_test_degenerated_naive_C2(p1.x(), p1.y(),
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
		      bool b, Cartesian_tag, Algebraic1_tag)
{
  typedef typename Pt::FT  FT;
  return ad_finite_edge_test_degenerated_alg1_C2(p1.x(), p1.y(),
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
		      bool b, Cartesian_tag, Algebraic2_tag)
{
  typedef typename Pt::FT  FT;
  return ad_finite_edge_test_degenerated_alg2_C2(p1.x(), p1.y(),
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
  typedef typename Pt::RT  RT;
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
  typedef typename Pt::R::Rep_tag Tag;
  return
    ad_finite_edge_test_2<Pt,We,Method_tag>(p1, p2, p3, q, b, Tag());
}

//-----------------------------------------------------------------------
// the functions for the Aw_Delaunay_finite_edge_test_2 predicate
//-----------------------------------------------------------------------



template < class Pt, class We >
inline
bool
ad_finite_edge_test_2(const Weighted_point< Pt, We >& p1,
		      const Weighted_point< Pt, We >& p2,
		      const Weighted_point< Pt, We >& p3,
		      const Weighted_point< Pt, We >& p4,
		      const Weighted_point< Pt, We >& q,
		      bool b, Cartesian_tag, Naive_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_finite_edge_test_naive_C2(p1.x(), p1.y(), FT(p1.weight()),
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
		      bool b, Cartesian_tag, Algebraic1_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_finite_edge_test_alg1_C2(p1.x(), p1.y(), FT(p1.weight()),
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
		      bool b, Cartesian_tag, Algebraic2_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_finite_edge_test_alg2_C2(p1.x(), p1.y(), FT(p1.weight()),
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
  typedef typename Pt::RT  RT;
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
  typedef typename Pt::R::Rep_tag Tag;
  return ad_finite_edge_test_2<Pt,We,Method_tag>
    (p1, p2, p3, p4, q, b, Tag());
}




template < class Point, class Weight, class Method_tag >
class Finite_edge_interior_conflict_2
{
public:
  typedef Weighted_point< Point, Weight >  Weighted_point;




  inline
  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& q, bool b)
  {
    return
      ad_finite_edge_test_2<Point,Weight,Method_tag>(p1, p2, p3, q, b);
  }

  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& q, bool b)
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
		  bool b)
  {
    return ad_finite_edge_test_2<Point,Weight,Method_tag>
      (p1, p2, p3, p4, q, b);
  }
};


//-----------------------------------------------------------------------
//-----------------------------------------------------------------------



template < class Pt, class We >
inline
bool
ad_is_degenerate_edge_test_2(const Weighted_point< Pt, We >& p1,
			     const Weighted_point< Pt, We >& p2,
			     const Weighted_point< Pt, We >& p3,
			     const Weighted_point< Pt, We >& p4,
			     Cartesian_tag, Naive_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_is_degenerate_edge_test_naive_C2(p1.x(), p1.y(), FT(p1.weight()),
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
			     Cartesian_tag, Algebraic1_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_is_degenerate_edge_test_alg1_C2(p1.x(), p1.y(), FT(p1.weight()),
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
			     Cartesian_tag, Algebraic2_tag)
{
  typedef typename Pt::FT  FT;
  return
    ad_is_degenerate_edge_test_alg2_C2(p1.x(), p1.y(), FT(p1.weight()),
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
  typedef typename Pt::RT  RT;
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
  typedef typename Pt::R::Rep_tag Tag;
  return ad_is_degenerate_edge_test_2<Pt,We,Method_tag>
    (p1, p2, p3, p4, Tag());
}



template < class Point, class Weight, class Method_tag >
class Is_degenerate_edge_2
{
public:
  typedef Weighted_point< Point, Weight >  Weighted_point;


  inline
  bool operator()(const Weighted_point& p1,
		  const Weighted_point& p2,
		  const Weighted_point& p3,
		  const Weighted_point& p4)
  {
    return ad_is_degenerate_edge_test_2<Point,Weight,Method_tag>
      (p1, p2, p3, p4);
  }
};




CGAL_END_NAMESPACE
