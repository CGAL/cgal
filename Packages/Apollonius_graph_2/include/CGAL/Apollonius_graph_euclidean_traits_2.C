// ======================================================================
//
// Copyright (c) 2003 The CGAL Consortium
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
// coordinator   :
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
template < class K >
inline
typename K::Point_2
ad_circumcenter_2(const typename K::Site_2& p,
		  const typename K::Site_2& q,
		  const typename K::Site_2& r,
		  Cartesian_tag )
{
  typename K::FT x,y;
  ad_circumcenterC2(p.x(),p.y(),p.weight(),
		    q.x(),q.y(),q.weight(),
		    r.x(),r.y(),r.weight(),x,y);
  return typename K::Point_2(x,y);
}

template < class K >
inline
typename K::Point_2
ad_circumcenter_2(const typename K::Site_2& p,
		  const typename K::Site_2& q,
		  const typename K::Site_2& r,
		  Homogeneous_tag )
{
  typename K::RT x,y,w;
  ad_circumcenterH2(p.hx(),p.hy(),p.hw(),p.weight(),
		    q.hx(),q.hy(),q.hw(),q.weight(),
		    r.hx(),r.hy(),r.hw(),r.weight(),
		    x,y,w);
  return typename K::Point_2(x,y,w);
}

template < class K >
inline
typename K::Point_2
ad_circumcenter_2(const typename K::Site_2& p,
		  const typename K::Site_2& q,
		  const typename K::Site_2& r)
{
  typedef typename K::Rep_tag Tag;
  return ad_circumcenter_2< K >(p, q, r, Tag()); 
}


template < class K >
class Construct_Apollonius_vertex_2
{
public:
  typedef typename K::Point_2              Point_2;
  typedef typename K::Site_2    Site_2;

  inline
  Point_2 operator() (const Site_2& p,
		      const Site_2& q,
		      const Site_2& r) const
  {
    //      CGAL_triangulation_precondition( ! collinear(p, q, r) );
    return ad_circumcenter_2< K >(p,q,r);
  }
};


//-----------------------------------------------------------------------
//                     Apollonius weighted point
//-----------------------------------------------------------------------
template < class K >
inline
typename K::Site_2
ad_circumcircle_2(const typename K::Site_2& p,
		  const typename K::Site_2& q,
		  const typename K::Site_2& r,
		  Cartesian_tag )
{
  typename K::FT x, y, wt;
  ad_circumcircleC2(p.x(),p.y(),p.weight(),
		    q.x(),q.y(),q.weight(),
		    r.x(),r.y(),r.weight(),x,y,wt);
  return typename K::Site_2(typename K::Point_2(x,y), wt);
}

template < class K >
inline
typename K::Site_2
ad_circumcircle_2(const typename K::Site_2& p,
		  const typename K::Site_2& q,
		  const typename K::Site_2& r,
		  Homogeneous_tag )
{
  typename K::RT x, y, w, wt;
  ad_circumcircleH2(p.hx(),p.hy(),p.hw(),p.weight(),
		    q.hx(),q.hy(),q.hw(),q.weight(),
		    r.hx(),r.hy(),r.hw(),r.weight(),
		    x,y,w,wt);
  return typename K::Site_2(typename K::Point_2(x,y,w), wt);
}

template < class K >
inline
typename K::Site_2
ad_circumcircle_2(const typename K::Site_2& p,
		  const typename K::Site_2& q,
		  const typename K::Site_2& r)
{
  typedef typename K::Rep_tag Tag;
  return ad_circumcircle_2< K >(p, q, r, Tag()); 
}

template < class K >
inline
typename K::Line_2
ad_left_bitangent_line_2(const typename K::Site_2& p,
			 const typename K::Site_2& q,
			 Cartesian_tag )
{
  typename K::FT a, b, c;
  ad_left_bitangent_lineC2(p.x(),p.y(),p.weight(),
			   q.x(),q.y(),q.weight(),
			   a,b,c);
  return typename K::Line_2(a, b, c);
}

template < class K >
inline
typename K::Line_2
ad_left_bitangent_line_2(const typename K::Site_2& p,
			 const typename K::Site_2& q,
			 Homogeneous_tag )
{
  typename K::RT a, b, c;
  ad_left_bitangent_lineH2(p.hx(),p.hy(),p.hw(),p.weight(),
			   q.hx(),q.hy(),q.hw(),q.weight(),
			   a, b, c);
  return typename K::Line_2(a, b, c);
}

template < class K >
inline
typename K::Line_2
ad_left_bitangent_line_2(const typename K::Site_2& p,
			 const typename K::Site_2& q)
{
  typedef typename K::Rep_tag Tag;
  return ad_left_bitangent_line_2< K >(p, q, Tag()); 
}


template < class K >
class Construct_Apollonius_weighted_point_2
{
public:
  typedef typename K::Line_2             Line_2;
  typedef typename K::Point_2            Point_2;
  typedef typename K::Site_2  Site_2;

  inline Site_2 operator()(const Site_2& p,
				      const Site_2& q,
				      const Site_2& r) const
  {
    //      CGAL_triangulation_precondition( ! collinear(p, q, r) );
    return ad_circumcircle_2< K >(p,q,r);
  }

  inline Line_2 operator()(const Site_2 &p,
			   const Site_2 &q) const
  {
    return ad_left_bitangent_line_2< K >(p, q);
  }
};


//-----------------------------------------------------------------------
//                        Apollonius bisector
//-----------------------------------------------------------------------


template< class K >
class Construct_Apollonius_bisector_2
{
public:
  typedef typename K::Point_2                Point_2;
  typedef typename K::RT                     Weight;
  typedef typename K::Line_2                 Line_2;
  typedef typename K::Site_2      Site_2;
  typedef typename K::Object_2               Object_2;
  typedef typename K::Construct_object_2     Construct_object_2;
  typedef CGAL::Hyperbola_2<Point_2, Weight> Hyperbola_2;

private:
  template<class T>
  Object_2 make_object(const T& t) const
  {
    return Construct_object_2()(t);
  }

public:
  inline Object_2 operator() (const Site_2& p,
			      const Site_2& q) const {
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Line_2 l1(p.point(), q.point());
      Line_2 l = l1.perpendicular(midpoint(p.point(), q.point()));
      return make_object(l);
    }

    Hyperbola_2 h(p, q);
    return make_object(h);
  }
};

//-----------------------------------------------------------------------
//                      Apollonius bisector ray
//-----------------------------------------------------------------------


template< class K >
class Construct_Apollonius_bisector_ray_2
{
public:
  typedef typename K::Point_2               Point_2;
  typedef typename K::RT                    Weight;
  typedef typename K::Line_2                Line_2;
  typedef typename K::Ray_2                 Ray_2;
  typedef typename K::Site_2     Site_2;
  typedef typename K::Object_2              Object_2;
  typedef typename K::Construct_object_2    Construct_object_2;
  typedef CGAL::Hyperbola_ray_2<Point_2, Weight>  Hyperbola_ray_2;
  typedef CGAL::Sign                              Hyperbola_direction;
  typedef CGAL::Construct_Apollonius_vertex_2<K>  Apollonius_vertex_2;

private:
  template<class T>
  Object_2 make_object(const T& t) const
  {
    return Construct_object_2()(t);
  }

public:
  inline Object_2
  operator() (const Site_2& p,
	      const Site_2& q,
	      const Point_2& r,
	      const Hyperbola_direction& direction) const {
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Line_2 l1(q, p);
      Line_2 l = l1.perpendicular(midpoint(p.point(), q.point()));
      Ray_2 ray(r, l.direction());
      return make_object(ray);
    }
    Hyperbola_ray_2 hr(p, q, r, direction);
    return make_object(hr);
  }

  inline Object_2
  operator() (const Site_2& p,
	      const Site_2& q,
	      const Site_2& r) const {
    Point_2 c = Apollonius_vertex_2()(p, q, r);
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Line_2 l1(q.point(), p.point());
      Line_2 l = l1.perpendicular(midpoint(p.point(), q.point()));
      Ray_2 ray(c, l.direction());
      return make_object(ray);
    }
    Hyperbola_ray_2 hr(p, q, c, NEGATIVE);
    return make_object(hr);
  }
};

//-----------------------------------------------------------------------
//                    Apollonius bisector segment
//-----------------------------------------------------------------------

template< class K >
class Construct_Apollonius_bisector_segment_2
{
public:
  typedef typename K::Point_2                 Point_2;
  typedef typename K::RT                      Weight;
  typedef typename K::Segment_2               Segment_2;
  typedef typename K::Site_2       Site_2;
  typedef typename K::Object_2                Object_2;
  typedef typename K::Construct_object_2      Construct_object_2;
  typedef CGAL::Hyperbola_segment_2<Point_2, Weight>  Hyperbola_segment_2;
  typedef CGAL::Construct_Apollonius_vertex_2<K>      Apollonius_vertex_2;

private:
  template<class T>
  Object_2 make_object(const T& t) const
  {
    return Construct_object_2()(t);
  } 

public:
  inline Object_2 operator() (const Site_2& p,
			      const Site_2& q,
			      const Point_2& r, const Point_2& s) const {
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Segment_2 seg(r, s);
      return make_object(seg);
    }
    Hyperbola_segment_2 hs(p, q, r, s);
    return make_object(hs);
  }

  inline Object_2 operator() (const Site_2& p,
			      const Site_2& q,
			      const Site_2& r,
			      const Site_2& s) const {
    Apollonius_vertex_2 apollonius_vertex_2;
    Point_2 c_pqr = apollonius_vertex_2(p,q,r);
    Point_2 c_qps = apollonius_vertex_2(q,p,s);
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Segment_2 seg(c_pqr, c_qps);
      return make_object(seg);
    }
    Hyperbola_segment_2 hs(p, q, c_pqr, c_qps);
    return make_object(hs);
  }

};

//-----------------------------------------------------------------------
//                    Apollonius primal ray
//-----------------------------------------------------------------------


template< class K >
class Construct_Apollonius_primal_ray_2
{
public:
  typedef typename K::Point_2                        Point_2;
  typedef typename K::RT                             Weight;
  typedef typename K::RT                             RT;
  typedef typename K::Line_2                         Line_2;
  typedef typename K::Ray_2                          Ray_2;
  typedef typename K::Site_2              Site_2;
  typedef CGAL::Construct_Apollonius_weighted_point_2<K> Apollonius_circle_2;

  inline Ray_2 operator() (const Site_2& p,
			   const Site_2& r,
			   const Site_2& s) const {
    //
    Apollonius_circle_2 apollonius_circle_2;
    Line_2 l1 = apollonius_circle_2(r, p);
    Line_2 l2 = apollonius_circle_2(p, s);

    RT d1 = CGAL_NTS sqrt( CGAL_NTS square(l1.a()) +
			   CGAL_NTS square(l1.b()) );
    RT d2 = CGAL_NTS sqrt( CGAL_NTS square(l2.a()) +
			   CGAL_NTS square(l2.b()) );
    RT a = l1.a() / d1 - l2.a() / d2;
    RT b = l1.b() / d1 - l2.b() / d2;
    Point_2 c(p.x() + b, p.y() - a);
    return Ray_2(p.point(), c);
  }
};

//-----------------------------------------------------------------------
//                    Apollonius primal segment
//-----------------------------------------------------------------------

template< class K >
class Construct_Apollonius_primal_segment_2
{
public:
  typedef typename K::Point_2                         Point_2;
  typedef typename K::RT                              Weight;
  typedef typename K::Line_2                          Line_2;
  typedef typename K::Segment_2                       Segment_2;
  typedef typename K::Site_2               Site_2;
  typedef typename K::Object_2                        Object_2;
  typedef typename K::Construct_object_2              Construct_object_2;
  typedef CGAL::Hyperbola_segment_2<Point_2,Weight>      Hyperbola_segment_2;
  typedef CGAL::Parabola_segment_2<Point_2,Weight,Line_2> 
  /*                                                   */ Parabola_segment_2;
  typedef CGAL::Construct_Apollonius_weighted_point_2<K> Apollonius_circle_2;


private:
  template<class T>
  Object_2 make_object(const T& t) const
  {
    return Construct_object_2()(t);
  }

public:

  inline Segment_2
  operator() (const Site_2& p,
	      const Site_2& q) const {
    //
    return Segment_2(p.point(), q.point());
  }

  inline Object_2
  operator() (const Site_2& p,
	      const Site_2& q,
	      const Point_2& r, const Point_2& s) const {
    //
    Comparison_result cr = CGAL_NTS compare(p.weight(), q.weight());
    if ( cr == EQUAL ) {
      Segment_2 seg(r, s);
      return make_object(seg);
    }
    Hyperbola_segment_2 hs(p, q, r, s);
    return make_object(hs);
  }

  inline Object_2
  operator() (const Site_2& p,
	      const Site_2& q,
	      const Site_2& r,
	      const Site_2& s) const {
    Apollonius_circle_2 apollonius_circle_2;
    Site_2 c_pqr = apollonius_circle_2(p, q, r);
    Site_2 c_qps = apollonius_circle_2(q, p, s);
    //
    Comparison_result cr = CGAL_NTS compare(c_pqr.weight(), c_qps.weight());
    if ( cr == EQUAL ) {
      Segment_2 seg(p.point(), q.point());
      return make_object(seg);
    }
    Hyperbola_segment_2 hs(c_pqr, c_qps, p.point(), q.point());
    return make_object(hs);
  }

  inline Parabola_segment_2
  operator() (const Site_2& p,
	      const Site_2& q,
	      const Site_2& r) const {
    //
    Apollonius_circle_2 apollonius_circle_2;
    Site_2 c = apollonius_circle_2(p, q, r);
    Line_2 l = apollonius_circle_2(q, p);
    return Parabola_segment_2(c, l, q.point(), p.point());
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

template < class K >
class Compare_weight_2
{
public:
  typedef typename K::Site_2  Site_2;
  typedef Comparison_result             result_type;

  inline
  Comparison_result operator()(const Site_2& p,
			       const Site_2& q) const
  {
    return CGAL_NTS compare(p.weight(), q.weight());
  }
};

//-----------------------------------------------------------------------
//                        Is hidden
//-----------------------------------------------------------------------

template < class K >
inline
bool
ad_is_hidden_test_2(const typename K::Site_2& p,
		    const typename K::Site_2& q,
		    Cartesian_tag, Sqrt_field_tag )
{
  return ad_is_hidden_test_sqrtf_C2(p.x(), p.y(), p.weight(),
				    q.x(), q.y(), q.weight());
}


template < class K >
inline
bool
ad_is_hidden_test_2(const typename K::Site_2& p,
		    const typename K::Site_2& q,
		    Cartesian_tag, Ring_tag )
{
  return ad_is_hidden_test_ring_C2(p.x(), p.y(), p.weight(),
				   q.x(), q.y(), q.weight());
}



template < class K, class Method_tag >
inline
bool
ad_is_hidden_test_2(const typename K::Site_2& p,
		    const typename K::Site_2& q,
		    Homogeneous_tag)
{
  Sign s = sign_of_ad_distance2_testH2(p.hx(), p.hy(), p.hw(), 
				       p.weight(),
				       q.hx(), q.hy(), q.hw(),
				       q.weight());
  if ( s == POSITIVE ) { return false; }
  return (CGAL_NTS compare(p.weight(), q.weight()) != SMALLER);
}

template < class K, class Method_tag >
inline
bool
ad_is_hidden_test_2(const typename K::Site_2& p,
		    const typename K::Site_2& q,
		    Cartesian_tag tag)
{
  return ad_is_hidden_test_2< K >(p, q, tag, Method_tag());
}


template< class K, class Method_tag >
class Is_hidden_2
{
public:
  typedef typename K::Site_2   Site_2;
  typedef bool                           result_type;

  inline bool operator()(const Site_2 &p,
			 const Site_2 &q) const
  {
    typedef typename K::Rep_tag Tag;
    return ad_is_hidden_test_2<K,Method_tag>(p, q, Tag());
  }
};


//-----------------------------------------------------------------------
//                    Oriented side of bisector
//-----------------------------------------------------------------------


template < class K >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p,
		    Cartesian_tag, Sqrt_field_tag )
{
  return
    compare_ad_distances_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				       p2.x(), p2.y(), p2.weight(),
				       p.x(),  p.y());
}


template < class K >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p,
		    Cartesian_tag, Ring_tag)
{
  return compare_ad_distances_test_ring_C2(p1.x(), p1.y(), p1.weight(),
					   p2.x(), p2.y(), p2.weight(),
					   p.x(),  p.y());
}



template < class K, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p, Cartesian_tag tag)
{
  return ad_distances_test_2< K >(p1, p2, p, tag, Method_tag());
}



template < class K, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p, Homogeneous_tag )
{
  return compare_ad_distances_testH2(p1.hx(), p1.hy(), p1.hw(),
				     p1.weight(),
				     p2.hx(), p2.hy(), p2.hw(),
				     p2.weight(),
				     p.hx(), p.hy(), p.hw());
}





template < class K, class Method_tag >
inline
Comparison_result
ad_distances_test_2(const typename K::Site_2& p1,
		    const typename K::Site_2& p2,
		    const typename K::Point_2& p)
{
  typedef typename K::Rep_tag Tag;
  return ad_distances_test_2<K,Method_tag>(p1, p2, p, Tag());
}



template< class K, class Method_tag >
class Oriented_side_of_bisector_2
{
public:
  typedef typename K::Point_2             Point_2;
  typedef typename K::Site_2    Site_2;
  typedef Oriented_side                   result_type;

  inline Oriented_side operator()(const Site_2& p1,
				  const Site_2& p2,
				  const Point_2 &p) const
  {
    Comparison_result r = ad_distances_test_2<K,Method_tag>(p1, p2, p);

    if ( r == EQUAL ) { return ON_ORIENTED_BOUNDARY; }
    return ( r == LARGER ) ? ON_NEGATIVE_SIDE : ON_POSITIVE_SIDE;
  }
};



//-----------------------------------------------------------------------
//                        Vertex conflict
//-----------------------------------------------------------------------


template < class K >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2&  q,
		   Cartesian_tag, Sqrt_field_tag )
{
  return ad_incircle_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				   p2.x(), p2.y(), p2.weight(),
				    q.x(),  q.y(),  q.weight());
}


template < class K >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2&  q,
		   Cartesian_tag, Ring_tag )
{
  return ad_incircle_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				  p2.x(), p2.y(), p2.weight(),
				   q.x(),  q.y(),  q.weight());
}


template < class K, class Method_tag >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2&  q,
		   Cartesian_tag tag)
{
  return ad_incircle_test_2< K >(p1, p2, q, tag, Method_tag());
}


template < class K, class Method_tag >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& q,
		   Homogeneous_tag )
{
  return 
    ad_incircle_testH2(p1.hx(), p1.hy(), p1.hw(), p1.weight(),
		       p2.hx(), p2.hy(), p2.hw(), p2.weight(),
		        q.hx(),  q.hy(),  q.hw(),  q.weight());
}


//-----------------------------------------------------------------------


template < class K >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& p3,
		   const typename K::Site_2&  q,
		   Cartesian_tag, Sqrt_field_tag )
{
  return ad_incircle_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				   p2.x(), p2.y(), p2.weight(),
				   p3.x(), p3.y(), p3.weight(),
				    q.x(),  q.y(),  q.weight());
}


template < class K >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& p3,
		   const typename K::Site_2&  q,
		   Cartesian_tag, Ring_tag )
{
  return ad_incircle_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				  p2.x(), p2.y(), p2.weight(),
				  p3.x(), p3.y(), p3.weight(),
				  q.x(),  q.y(),   q.weight());
}



template < class K, class Method_tag >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& p3,
		   const typename K::Site_2&  q,
		   Cartesian_tag tag)
{
  return ad_incircle_test_2< K >(p1, p2, p3, q, tag, Method_tag());
}


template < class K, class Method_tag >
inline
Sign
ad_incircle_test_2(const typename K::Site_2& p1,
		   const typename K::Site_2& p2,
		   const typename K::Site_2& p3,
		   const typename K::Site_2& q,
		   Homogeneous_tag )
{
  return 
    ad_incircle_testH2(p1.hx(), p1.hy(), p1.hw(), p1.weight(),
		       p2.hx(), p2.hy(), p2.hw(), p2.weight(),
		       p3.hx(), p3.hy(), p3.hw(), p3.weight(),
		        q.hx(),  q.hy(),  q.hw(),  q.weight());
}


template < class K, class Method_tag >
class Vertex_conflict_2
{
public:
  typedef typename K::Site_2      Site_2;
  typedef Sign                               result_type;

  inline
  Sign operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& q) const
  {
    typedef typename K::Rep_tag Tag;
    return ad_incircle_test_2<K,Method_tag>(p1, p2, p3, q, Tag());
  }


  inline
  Sign operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& q) const
  {
    typedef typename K::Rep_tag Tag;
    return ad_incircle_test_2<K,Method_tag>(p1, p2, q, Tag());
  }
 


};

//-----------------------------------------------------------------------
//                    Finite edge interior conflict
//-----------------------------------------------------------------------

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  return
    ad_finite_edge_test_degenerated_sqrtf_C2(p1.x(), p1.y(),
					     p1.weight(),
					     p2.x(), p2.y(),
					     p2.weight(),
					     q.x(),  q.y(),
					     q.weight(), b);
}

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  return
    ad_finite_edge_test_degenerated_ring_C2(p1.x(), p1.y(),
					    p1.weight(),
					    p2.x(), p2.y(),
					    p2.weight(),
					    q.x(),  q.y(),
					    q.weight(), b);
}


template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag tag)
{
  return
    ad_finite_edge_test_2< K >(p1, p2, q, b, tag, Method_tag());
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q,
		      bool b, Homogeneous_tag)
{
  return
    ad_finite_edge_test_degeneratedH2(p1.hx(), p1.hy(),
				      p1.hw(),
				      p1.weight(),
				      p2.hx(), p2.hy(),
				      p2.hw(),
				      p2.weight(),
				      q.hx(),  q.hy(),
				      q.hw(),
				      q.weight(), b);
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& q, bool b)
{
  typedef typename K::Rep_tag Tag;
  return ad_finite_edge_test_2<K,Method_tag>(p1, p2, q, b, Tag());
}

//-----------------------------------------------------------------------

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  return ad_finite_edge_test_degenerated_sqrtf_C2(p1.x(), p1.y(),
						  p1.weight(),
						  p2.x(), p2.y(),
						  p2.weight(),
						  p3.x(), p3.y(),
						  p3.weight(),
						  q.x(),  q.y(),
						  q.weight(), b);
}

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  return ad_finite_edge_test_degenerated_ring_C2(p1.x(), p1.y(),
						 p1.weight(),
						 p2.x(), p2.y(),
						 p2.weight(),
						 p3.x(), p3.y(),
						 p3.weight(),
						 q.x(),  q.y(),
						 q.weight(), b);
}


template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag tag)
{
  return
    ad_finite_edge_test_2< K >(p1, p2, p3, q, b, tag, Method_tag());
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q,
		      bool b, Homogeneous_tag)
{
  return
    ad_finite_edge_test_degeneratedH2(p1.hx(), p1.hy(),
				      p1.hw(),
				      p1.weight(),
				      p2.hx(), p2.hy(),
				      p2.hw(),
				      p2.weight(),
				      p3.hx(), p3.hy(),
				      p3.hw(),
				      p3.weight(),
				      q.hx(),  q.hy(),
				      q.hw(),
				      q.weight(), b);
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& q, bool b)
{
  typedef typename K::Rep_tag Tag;
  return ad_finite_edge_test_2<K,Method_tag>(p1, p2, p3, q, b, Tag());
}

//-----------------------------------------------------------------------


template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Sqrt_field_tag)
{
  return
    ad_finite_edge_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
				 p2.x(), p2.y(), p2.weight(),
				 p3.x(), p3.y(), p3.weight(),
				 p4.x(), p4.y(), p4.weight(),
				 q.x(),  q.y(), q.weight(), b);
}

template < class K >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag, Ring_tag)
{
  return
    ad_finite_edge_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				p2.x(), p2.y(), p2.weight(),
				p3.x(), p3.y(), p3.weight(),
				p4.x(), p4.y(), p4.weight(),
				 q.x(),  q.y(),  q.weight(), b);
}


template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b, Cartesian_tag tag)
{
  return ad_finite_edge_test_2< K >(p1, p2, p3, p4, q, b,
				    tag, Method_tag());
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b, Homogeneous_tag)
{
  return
    ad_Voronoi_diagram_finite_edge_testH2(p1.hx(), p1.hy(), p1.hw(),
					  p1.weight(),
					  p2.hx(), p2.hy(), p2.hw(),
					  p2.weight(),
					  p3.hx(), p3.hy(), p3.hw(),
					  p3.weight(),
					  p4.hx(), p4.hy(), p4.hw(),
					  p4.weight(),
					  q.hx(),  q.hy(), q.hw(),
					  q.weight(), b);
}

template < class K, class Method_tag >
inline
bool
ad_finite_edge_test_2(const typename K::Site_2& p1,
		      const typename K::Site_2& p2,
		      const typename K::Site_2& p3,
		      const typename K::Site_2& p4,
		      const typename K::Site_2& q,
		      bool b)
{
  typedef typename K::Rep_tag Tag;
  return ad_finite_edge_test_2<K,Method_tag>
    (p1, p2, p3, p4, q, b, Tag());
}




template < class K, class Method_tag >
class Finite_edge_interior_conflict_2
{
public:
  typedef typename K::Site_2  Site_2;
  typedef bool                           result_type;

  inline
  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& q, bool b) const
  {
    return ad_finite_edge_test_2<K,Method_tag>(p1, p2, p3, q, b);
  }

  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& q, bool b) const
  {
    return ad_finite_edge_test_2<K,Method_tag>(p1, p2, q, b);
  }



  inline
  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& p4,
		  const Site_2& q,
		  bool b) const
  {
    return ad_finite_edge_test_2<K,Method_tag>(p1, p2, p3, p4, q, b);
  }
};


//-----------------------------------------------------------------------
//                   Infinite edge interior conflict
//-----------------------------------------------------------------------

template < class K >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q,
			bool b, Cartesian_tag, Sqrt_field_tag)
{
  return
    ad_infinite_edge_test_sqrtf_C2(p2.x(), p2.y(), p2.weight(),
				   p3.x(), p3.y(), p3.weight(),
				   p4.x(), p4.y(), p4.weight(),
				    q.x(),  q.y(),  q.weight(), b);
}


template < class K >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q,
			bool b, Cartesian_tag, Ring_tag)
{
  return
    ad_infinite_edge_test_ring_C2(p2.x(), p2.y(), p2.weight(),
				  p3.x(), p3.y(), p3.weight(),
				  p4.x(), p4.y(), p4.weight(),
				   q.x(),  q.y(),  q.weight(), b);
}


template < class K, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q,
			bool b, Cartesian_tag tag)
{
  return
    ad_infinite_edge_test_2<K>(p2, p3, p4, q, b, tag, Method_tag());
}



template < class K, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q,
			bool b, Homogeneous_tag)
{
  return
    ad_infinite_edge_testH2(p2.hx(), p2.hy(), p2.hw(),
			    p2.weight(),
			    p3.hx(), p3.hy(), p3.hw(),
			    p3.weight(),
			    p4.hx(), p4.hy(), p4.hw(),
			    p4.weight(),
			    q.hx(),  q.hy(), q.hw(),
			    q.weight(), b);
}

template < class K, class Method_tag >
inline
bool
ad_infinite_edge_test_2(const typename K::Site_2& p2,
			const typename K::Site_2& p3,
			const typename K::Site_2& p4,
			const typename K::Site_2& q, bool b)
{
  typedef typename K::Rep_tag Tag;
  return ad_infinite_edge_test_2<K,Method_tag>(p2, p3, p4, q, b, Tag());
}


template < class K, class Method_tag >
class Infinite_edge_interior_conflict_2
{
public:
  typedef typename K::Site_2 Site_2;
  typedef bool                          result_type;

  inline
  bool operator()(const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& p4,
		  const Site_2& q, bool b) const
  {
    return ad_infinite_edge_test_2<K,Method_tag>(p2, p3, p4, q, b);
  }
};


//-----------------------------------------------------------------------
//                          Is degenerate
//-----------------------------------------------------------------------



template < class K >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4,
			     Cartesian_tag, Sqrt_field_tag)
{
  return
    ad_is_degenerate_edge_test_sqrtf_C2(p1.x(), p1.y(), p1.weight(),
					p2.x(), p2.y(), p2.weight(),
					p3.x(), p3.y(), p3.weight(),
					p4.x(), p4.y(), p4.weight());
}

template < class K >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4,
			     Cartesian_tag, Ring_tag)
{
  return
    ad_is_degenerate_edge_test_ring_C2(p1.x(), p1.y(), p1.weight(),
				       p2.x(), p2.y(), p2.weight(),
				       p3.x(), p3.y(), p3.weight(),
				       p4.x(), p4.y(), p4.weight());
}


template < class K, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4,
			     Cartesian_tag tag)
{
  return
    ad_is_degenerate_edge_test_2< K >(p1, p2, p3, p4, tag, Method_tag());
}


template < class K, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4,
			     Homogeneous_tag)
{
  return
    ad_is_degenerate_edge_testH2(p1.hx(), p1.hy(), p1.hw(),
				 p1.weight(),
				 p2.hx(), p2.hy(), p2.hw(),
				 p2.weight(),
				 p3.hx(), p3.hy(), p3.hw(),
				 p3.weight(),
				 p4.hx(), p4.hy(), p4.hw(),
				 p4.weight());
}

template < class K, class Method_tag >
inline
bool
ad_is_degenerate_edge_test_2(const typename K::Site_2& p1,
			     const typename K::Site_2& p2,
			     const typename K::Site_2& p3,
			     const typename K::Site_2& p4)
{
  typedef typename K::Rep_tag Tag;
  return
    ad_is_degenerate_edge_test_2<K,Method_tag>(p1, p2, p3, p4, Tag());
}



template < class K, class Method_tag >
class Is_degenerate_edge_2
{
public:
  typedef typename K::Site_2    Site_2;
  typedef bool                             result_type;

  inline
  bool operator()(const Site_2& p1,
		  const Site_2& p2,
		  const Site_2& p3,
		  const Site_2& p4) const
  {
    return
      ad_is_degenerate_edge_test_2<K,Method_tag>(p1, p2, p3, p4);
  }
};




CGAL_END_NAMESPACE
