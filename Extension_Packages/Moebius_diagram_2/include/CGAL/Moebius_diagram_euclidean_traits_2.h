#ifndef CHR_MOEBIUS_DIAGRAM_EUCLIDEAN_TRAITS_2_H
#define CHR_MOEBIUS_DIAGRAM_EUCLIDEAN_TRAITS_2_H


#include <CGAL/basic.h>
#include <CGAL/Point_2.h>
#include <CGAL/Line_2.h>
#include <CGAL/Circle_2.h>
#include <CGAL/Conic_2.h>
#include <CGAL/Object.h>


#include <CGAL/Moebius_utils.h>
#include <CGAL/Moebius_regular_euclidean_traits_3.h>
#include <CGAL/Moebius_point.h>


#include <CGAL/predicates/Moebius_diagram_ftC2.h>
#include <CGAL/constructions/Moebius_diagram_ftC2.h>


CGAL_BEGIN_NAMESPACE

// ----------------- predicates ---------------------

// ------- compare on oriented line ------

template <class R>
inline
bool
compare_on_oriented_line_2 (const Line_2<R> &l, const Point_2<R> &p, const Point_2<R> &q,
			    Cartesian_tag)
{
  typedef typename R::FT FT;
  return compare_on_oriented_lineC2 (l.a(),l.b(), p.x(),p.y(), q.x(),q.y());
}


template <class R>
inline
bool
compare_on_oriented_line_2 (const Line_2<R> &l, const Point_2<R> &p, const Point_2<R> &q,
			    Homogeneous_tag)
{
  typedef typename R::FT FT;
  return compare_on_oriented_lineC2 (l.a(),l.b(), p.x(),p.y(), q.x(),q.y());
}


template <class R>
inline
bool
compare_on_oriented_line_2 (const Line_2<R> &l, const Point_2<R> &p, const Point_2<R> &q)
{
  typedef typename R::Rep_tag Tag;
  return compare_on_oriented_line_2 (l, p, q, Tag ());
}

template <class R>
class Compare_on_oriented_line_2
{
 public:
  typedef Point_2<R> Point;
  typedef Line_2<R> Line;
  typedef /* blip blop */ bool result_type;

  result_type
    operator() (const Line &l, const Point &p, const Point &q) const
    { 
      return compare_on_oriented_line_2 (l, p, q);
    }
};

// ------- moebius orientation -------
template <class P, class W>
inline
CGAL::Orientation
moebius_orientation_2 (const Moebius_point<P,W> &p,
		      const Moebius_point<P,W> &q,
		      Cartesian_tag)
{
  typedef typename P::R::RT RT;
  return moebius_orientationC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
			       q.x(), q.y(), RT(q.lambda()), RT (q.mu()));
}
template <class P, class W>
inline
CGAL::Orientation
moebius_orientation_2 (const Moebius_point<P,W> &p,
		      const Moebius_point<P,W> &q,
		      Homogeneous_tag)
{
  typedef typename P::R::RT RT;
  return moebius_orientationC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
			       q.x(), q.y(), RT(q.lambda()), RT (q.mu()));
}
template <class P, class W>
inline
CGAL::Orientation
moebius_orientation_2 (const Moebius_point<P,W> &p,
		      const Moebius_point<P,W> &q)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_orientation_2 (p, q, Tag ());
}
template <class P, class W>
class Moebius_orientation_2
{
 public:
  typedef Moebius_point<P,W> Point;
  typedef CGAL::Orientation result_type;

  result_type
    operator() (const Point &p, const Point &q) const
    { 
      return moebius_orientation_2 (p, q);
    }
};

// ------- moebius has circle -------
template <class P, class W>
inline
bool
moebius_has_circle_2 (const Moebius_point<P,W> &p,
		     const Moebius_point<P,W> &q,
		     Cartesian_tag)
{
  typedef typename P::R::RT RT;
  return moebius_has_circleC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
			      q.x(), q.y(), RT(q.lambda()), RT (q.mu()));
}
template <class P, class W>
inline
bool
moebius_has_circle_2 (const Moebius_point<P,W> &p,
		     const Moebius_point<P,W> &q,
		     Homogeneous_tag)
{
  typedef typename P::R::RT RT;
  return moebius_has_circleC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
			      q.x(), q.y(), RT(q.lambda()), RT (q.mu()));
}
template <class P, class W>
inline
bool
moebius_has_circle_2 (const Moebius_point<P,W> &p,
		     const Moebius_point<P,W> &q)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_has_circle_2 (p, q, Tag ());
}
template <class P, class W>
class Moebius_has_circle_2
{
 public:
  typedef Moebius_point<P,W> Point;
  typedef bool result_type;

  result_type
    operator() (const Point &p, const Point &q) const
    { 
      return moebius_has_circle_2 (p, q);
    }
};

// ------- moebius circle cross line -------
template <class P, class W>
inline
bool
moebius_circle_cross_line_2 (const Moebius_point<P,W> &p,
			    const Moebius_point<P,W> &q,
			    const Moebius_point<P,W> &r,
			    Cartesian_tag)
{
  typedef typename P::R::RT RT;
  return moebius_circle_cross_lineC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
				     q.x(), q.y(), RT(q.lambda()), RT (q.mu()),
				     r.x(), r.y(), RT(r.lambda()), RT (r.mu()));
}
template <class P, class W>
inline
bool
moebius_circle_cross_line_2 (const Moebius_point<P,W> &p,
			    const Moebius_point<P,W> &q,
			    const Moebius_point<P,W> &r,
			    Homogeneous_tag)
{
  typedef typename P::R::RT RT;
  return moebius_circle_cross_lineC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
				     q.x(), q.y(), RT(q.lambda()), RT (q.mu()),
				     r.x(), r.y(), RT(r.lambda()), RT (r.mu()));
}
template <class P, class W>
inline
bool
moebius_circle_cross_line_2 (const Moebius_point<P,W> &p,
			    const Moebius_point<P,W> &q,
			    const Moebius_point<P,W> &r)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_circle_cross_line_2 (p, q, r, Tag ());
}
template <class P, class W>
class Moebius_circle_cross_line_2
{
 public:
  typedef Moebius_point<P,W> Point;
  typedef bool result_type;

  result_type
    operator() (const Point &p, const Point &q, const Point &r) const
    { 
      return moebius_circle_cross_line_2 (p, q, r);
    }
};

// ------- moebius circle side of center -------
template <class P, class W>
inline
CGAL::Oriented_side
moebius_circle_side_of_center_2 (const Moebius_point<P,W> &p,
				const Moebius_point<P,W> &q,
				const Moebius_point<P,W> &r,
				Cartesian_tag)
{
  typedef typename P::R::RT RT;
  return moebius_circle_side_of_centerC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
					 q.x(), q.y(), RT(q.lambda()), RT (q.mu()),
					 r.x(), r.y(), RT(r.lambda()), RT (r.mu()));
}
template <class P, class W>
inline
CGAL::Oriented_side
moebius_circle_side_of_center_2 (const Moebius_point<P,W> &p,
				const Moebius_point<P,W> &q,
				const Moebius_point<P,W> &r,
				Homogeneous_tag)
{
  typedef typename P::R::RT RT;
  return moebius_circle_side_of_centerC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
					 q.x(), q.y(), RT(q.lambda()), RT (q.mu()),
					 r.x(), r.y(), RT(r.lambda()), RT (r.mu()));
}
template <class P, class W>
inline
CGAL::Oriented_side
moebius_circle_side_of_center_2 (const Moebius_point<P,W> &p,
				const Moebius_point<P,W> &q,
				const Moebius_point<P,W> &r)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_circle_side_of_center_2 (p, q, r, Tag ());
}
template <class P, class W>
class Moebius_circle_side_of_center_2
{
 public:
  typedef Moebius_point<P,W> Point;
  typedef CGAL::Oriented_side result_type;

  result_type
    operator() (const Point &p, const Point &q, const Point &r) const
    { 
      return moebius_circle_side_of_center_2 (p, q, r);
    }
};

// ------- moebius side of vertex -------
template <class P, class W>
inline
CGAL::Bounded_side
moebius_side_of_vertex_2 (const Moebius_point<P,W> &p,
			 const Moebius_point<P,W> &q,
			 const Moebius_point<P,W> &r,
			 const Moebius_point<P,W> &s,
			 Cartesian_tag)
{
  typedef typename P::R::RT RT;
  return moebius_side_of_vertexC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
				  q.x(), q.y(), RT(q.lambda()), RT (q.mu()),
				  r.x(), r.y(), RT(r.lambda()), RT (r.mu()),
				  s.x(), s.y(), RT(s.lambda()), RT (s.mu()));
}
template <class P, class W>
inline
CGAL::Bounded_side
moebius_side_of_vertex_2 (const Moebius_point<P,W> &p,
			 const Moebius_point<P,W> &q,
			 const Moebius_point<P,W> &r,
			 const Moebius_point<P,W> &s,
			 Homogeneous_tag)
{
  typedef typename P::R::RT RT;
  return moebius_side_of_vertexC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
				  q.x(), q.y(), RT(q.lambda()), RT (q.mu()),
				  r.x(), r.y(), RT(r.lambda()), RT (r.mu()),
				  s.x(), s.y(), RT(s.lambda()), RT (s.mu()));
}
template <class P, class W>
inline
CGAL::Bounded_side
moebius_side_of_vertex_2 (const Moebius_point<P,W> &p,
			 const Moebius_point<P,W> &q,
			 const Moebius_point<P,W> &r,
			 const Moebius_point<P,W> &s)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_side_of_vertex_2 (p, q, r, s, Tag ());
}
template <class P, class W>
class Moebius_side_of_vertex_2
{
 public:
  typedef Moebius_point<P,W> Point;
  typedef CGAL::Bounded_side result_type;

  result_type
    operator() (const Point &p, const Point &q, const Point &r, const Point &s) const
    { 
      return moebius_side_of_vertex_2 (p, q, r, s);
    }
};

// ------- moebius side of vertex -------
template <class P, class W>
inline
CGAL::Oriented_side
moebius_circle_side_of_vertex_2 (const Moebius_point<P,W> &p,
				const Moebius_point<P,W> &q,
				const Moebius_point<P,W> &r,
				const Moebius_point<P,W> &s,
				Cartesian_tag)
{
  typedef typename P::R::RT RT;
  return moebius_circle_side_of_vertexC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
					 q.x(), q.y(), RT(q.lambda()), RT (q.mu()),
					 r.x(), r.y(), RT(r.lambda()), RT (r.mu()),
					 s.x(), s.y(), RT(s.lambda()), RT (s.mu()));
}
template <class P, class W>
inline
CGAL::Oriented_side
moebius_circle_side_of_vertex_2 (const Moebius_point<P,W> &p,
				const Moebius_point<P,W> &q,
				const Moebius_point<P,W> &r,
				const Moebius_point<P,W> &s,
				Homogeneous_tag)
{
  typedef typename P::R::RT RT;
  return moebius_circle_side_of_vertexC2 (p.x(), p.y(), RT(p.lambda()), RT (p.mu()),
					 q.x(), q.y(), RT(q.lambda()), RT (q.mu()),
					 r.x(), r.y(), RT(r.lambda()), RT (r.mu()),
					 s.x(), s.y(), RT(s.lambda()), RT (s.mu()));
}
template <class P, class W>
inline
CGAL::Oriented_side
moebius_circle_side_of_vertex_2 (const Moebius_point<P,W> &p,
				const Moebius_point<P,W> &q,
				const Moebius_point<P,W> &r,
				const Moebius_point<P,W> &s)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_circle_side_of_vertex_2 (p, q, r, s, Tag ());
}
template <class P, class W>
class Moebius_circle_side_of_vertex_2
{
 public:
  typedef Moebius_point<P,W> Point;
  typedef CGAL::Oriented_side result_type;

  result_type
    operator() (const Point &p, const Point &q, const Point &r, const Point &s) const
    { 
      return moebius_circle_side_of_vertex_2 (p, q, r, s);
    }
};


// ----------------- constructions ---------------------

// ------- moebius bisector -------
template <class P, class W>
inline
CGAL::Line_2<typename P::R>
moebius_line_2 (const Moebius_point<P,W> &p,
	       const Moebius_point<P,W> &q,
	       Cartesian_tag)
{
  typedef typename P::R::FT FT;
  typedef CGAL::Line_2<typename P::R> Line;
  FT a, b, c;
  moebius_lineC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
		 q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
		 a, b, c);
  return Line (a, b, c);
}

template <class P, class W>
inline
CGAL::Line_2<typename P::R>
moebius_line_2 (const Moebius_point<P,W> &p,
	       const Moebius_point<P,W> &q,
	       Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  typedef CGAL::Line_2<typename P::R> Line;
  FT a, b, c;
  moebius_lineC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
		 q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
		 a, b, c);
  return Line (a, b, c);
}

template <class P, class W>
inline
CGAL::Line_2<typename P::R>
moebius_line_2 (const Moebius_point<P,W> &p,
	       const Moebius_point<P,W> &q)
{
  typedef typename P::R::Rep_tag Tag;
  CGAL::Line_2<typename P::R> l = moebius_line_2 (p, q, Tag ());
  TRACE ("moebius_line_2 ("<<p<<","<<q<<") = "<<l<<"\n");
  return l;
}

template <class P, class W>
class Moebius_construct_line_2
{
 public:
  typedef CGAL::Moebius_point<P, W> Point;
  typedef CGAL::Line_2<typename P::R> Line;
  typedef /* blip blop */ Line result_type;

  result_type
    operator() (const Point &p,
		const Point &q) const
    { 
      return moebius_line_2 (p, q);
    }
};

template <class P, class W>
inline
CGAL::Circle_2<typename P::R>
moebius_circle_2 (const Moebius_point<P,W> &p,
		 const Moebius_point<P,W> &q,
		 Cartesian_tag)
{
  typedef typename P::R::FT FT;
  typedef CGAL::Circle_2<typename P::R> Circle;
  FT x, y, r;
  CGAL::Orientation o = moebius_circleC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
					 q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
					 x, y, r);
  return Circle (P (x, y), r, o);
}

template <class P, class W>
inline
CGAL::Circle_2<typename P::R>
moebius_circle_2 (const Moebius_point<P,W> &p,
		 const Moebius_point<P,W> &q,
		 Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  typedef CGAL::Circle_2<typename P::R> Circle;
  FT x, y, r;
  CGAL::Orientation o = moebius_circleC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
					 q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
					 x, y, r);
  return Circle (P (x, y), r, o);
}

template <class P, class W>
inline
CGAL::Circle_2<typename P::R>
moebius_circle_2 (const Moebius_point<P,W> &p,
		 const Moebius_point<P,W> &q)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_circle_2 (p, q, Tag ());
}

template <class P, class W>
class Moebius_construct_circle_2
{
 public:
  typedef CGAL::Moebius_point<P, W> Point;
  typedef CGAL::Circle_2<typename P::R> Circle;
  typedef /* blip blop */ Circle result_type;

  result_type
    operator() (const Point &p,
		const Point &q) const
    { 
      return moebius_circle_2<P,W> (p, q);
    }
};

// ------- custom intersections --------

template <class K, class Pair>
inline
Object
moebius_circle_2_line_2_intersection_2 (const Circle_2<K> &c,
				       const Line_2<K> &l,
				       Cartesian_tag)
{
  typedef typename K::FT FT;
  FT x1, y1, x2, y2;
  if (moebius_circle_2_line_2_intersectionC2 (c.center().x(), c.center().y(), c.squared_radius(),
					     l.a(), l.b(), l.c(),
					     x1, y1, x2, y2) == 2)
    return make_object (Pair (Point_2<K> (x1, y1), Point_2<K> (x2, y2)));
  return Object ();
}

template <class K, class Pair>
inline
Object
moebius_circle_2_line_2_intersection_2 (const Circle_2<K> &c,
				       const Line_2<K> &l,
				       Homogeneous_tag)
{
  typedef typename K::FT FT;
  FT x1, y1, x2, y2;
  if (moebius_circle_2_line_2_intersectionC2 (c.center().x(), c.center().y(), c.squared_radius(),
					     l.a(), l.b(), l.c(),
					     x1, y1, x2, y2) == 2)
    return make_object (Pair (Point_2<K> (x1, y1), Point_2<K> (x2, y2)));
  return make_object ();
}

template <class K, class Pair>
inline
Object
moebius_circle_2_line_2_intersection_2 (const Circle_2<K> &c,
				       const Line_2<K> &l)
{
  typedef typename K::Rep_tag Tag;
  return moebius_circle_2_line_2_intersection_2<K,Pair> (c, l, Tag ());
}

template <class K, class Pair>
class Moebius_circle_2_line_2_intersection_2
{
 public:
  typedef Object result_type;
  Object operator () (const Circle_2<K> &c, const Line_2<K> &l)
    { return moebius_circle_2_line_2_intersection_2<K,Pair> (c, l); }
};

template <class K, class Pair>
inline
Object
moebius_circle_2_circle_2_intersection_2 (const Circle_2<K> &c1,
					 const Circle_2<K> &c2,
					 Cartesian_tag)
{
  typedef typename K::FT FT;
  FT x1, y1, x2, y2;
  if (moebius_circle_2_circle_2_intersectionC2 (c1.center().x(), c1.center().y(), c1.squared_radius(),
					       c2.center().x(), c2.center().y(), c2.squared_radius(),
					       x1, y1, x2, y2) == 2)
    return make_object (Pair (Point_2<K> (x1, y1), Point_2<K> (x2, y2)));
  return Object ();
}

template <class K, class Pair>
inline
Object
moebius_circle_2_circle_2_intersection_2 (const Circle_2<K> &c1,
					 const Circle_2<K> &c2,
					 Homogeneous_tag)
{
  typedef typename K::FT FT;
  FT x1, y1, x2, y2;
  if (moebius_circle_2_circle_2_intersectionC2 (c1.center().x(), c1.center().y(), c1.squared_radius(),
					       c2.center().x(), c2.center().y(), c2.squared_radius(),
					       x1, y1, x2, y2) == 2)
    return make_object (Pair (Point_2<K> (x1, y1), Point_2<K> (x2, y2)));
  return make_object ();
}

template <class K, class Pair>
inline
Object
moebius_circle_2_circle_2_intersection_2 (const Circle_2<K> &c,
					 const Circle_2<K> &l)
{
  typedef typename K::Rep_tag Tag;
  return moebius_circle_2_circle_2_intersection_2<K,Pair> (c, l, Tag ());
}

template <class K, class Pair>
class Moebius_circle_2_circle_2_intersection_2
{
 public:
  typedef Object result_type;
  Object operator () (const Circle_2<K> &c, const Circle_2<K> &l)
    { return moebius_circle_2_circle_2_intersection_2<K,Pair> (c, l); }
};


// ------- moebius vertices -------


template <class P, class W, class Pair>
inline
Pair
moebius_circle_vertices_2 (const Moebius_point<P,W> &p,
			  const Moebius_point<P,W> &q,
			  const Moebius_point<P,W> &r,
			  Cartesian_tag)
{
  typedef typename P::R::FT FT;
  FT x1, y1, x2, y2;
  moebius_circle_verticesC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			    q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			    r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
			    x1, y1, x2, y2);
  return Pair (P (x1,y1), P (x2, y2));
}
 
template <class P, class W, class Pair>
inline
Pair
moebius_circle_vertices_2 (const Moebius_point<P,W> &p,
			  const Moebius_point<P,W> &q,
			  const Moebius_point<P,W> &r,
			  Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  FT x1, y1, x2, y2;
  moebius_circle_verticesC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			    q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			    r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
			    x1, y1, x2, y2);
  return Pair (P (x1,y1), P (x2, y2));
}

template <class P, class W, class Pair>
inline
Pair
moebius_circle_vertices_2 (const Moebius_point<P,W> &p,
			  const Moebius_point<P,W> &q,
			  const Moebius_point<P,W> &r)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_circle_vertices_2<P,W,Pair> (p, q, r, Tag ());
}


template <class P, class W, class Pair>
class Moebius_construct_circle_vertices_2
{
 public:
  typedef CGAL::Moebius_point<P, W> Point;
  typedef /* blip blop */ Pair result_type;

  result_type
    operator() (const Point &p,
		const Point &q,
    		const Point &r) const
    { 
      return moebius_circle_vertices_2<P,W,Pair> (p, q, r);
    }
};


template <class P, class W, class Pair>
inline
Pair
moebius_line_vertices_2 (const Moebius_point<P,W> &p,
			const Moebius_point<P,W> &q,
			const Moebius_point<P,W> &r,
			Cartesian_tag)
{
  typedef typename P::R::FT FT;
  FT x1, y1, x2, y2;
  moebius_line_verticesC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			  q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			  r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
			  x1, y1, x2, y2);
  return Pair (P (x1,y1), P (x2, y2));
}
 
template <class P, class W, class Pair>
inline
Pair
moebius_line_vertices_2 (const Moebius_point<P,W> &p,
			const Moebius_point<P,W> &q,
			const Moebius_point<P,W> &r,
			Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  FT x1, y1, x2, y2;
  moebius_line_verticesC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			  q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			  r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
			  x1, y1, x2, y2);
  return Pair (P (x1,y1), P (x2, y2));
}

template <class P, class W, class Pair>
inline
Pair
moebius_line_vertices_2 (const Moebius_point<P,W> &p,
			const Moebius_point<P,W> &q,
			const Moebius_point<P,W> &r)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_line_vertices_2<P,W,Pair> (p, q, r, Tag ());
}


template <class P, class W, class Pair>
class Moebius_construct_line_vertices_2
{
 public:
  typedef CGAL::Moebius_point<P, W> Point;
  typedef /* blip blop */ Pair result_type;

  result_type
    operator() (const Point &p,
		const Point &q,
    		const Point &r) const
    { 
      return moebius_line_vertices_2<P,W,Pair> (p, q, r);
    }
};


template <class P, class W>
inline
P
moebius_vertex_2 (const Moebius_point<P,W> &p,
		 const Moebius_point<P,W> &q,
		 const Moebius_point<P,W> &r,
		 Cartesian_tag)
{
  typedef typename P::R::FT FT;
  FT x, y;
  moebius_vertexC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
		   q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
		   r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
		   x, y);
  return P (x,y);
}
 
template <class P, class W>
inline
P
moebius_vertex_2 (const Moebius_point<P,W> &p,
		 const Moebius_point<P,W> &q,
		 const Moebius_point<P,W> &r,
		 Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  FT x, y;
  moebius_vertexC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
		   q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
		   r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
		   x, y);
  return P (x,y);
}

template <class P, class W>
inline
P
moebius_vertex_2 (const Moebius_point<P,W> &p,
		 const Moebius_point<P,W> &q,
		 const Moebius_point<P,W> &r)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_vertex_2 (p, q, r, Tag ());
}

template <class P, class W>
class Moebius_construct_vertex_2
{
 public:
  typedef CGAL::Moebius_point<P, W> Point;
  typedef /* blip blop */ P result_type;

  result_type
    operator() (const Point &p,
		const Point &q,
    		const Point &r) const
    { 
      return moebius_vertex_2 (p, q, r);
    }
};


// ------- regular circum center -------


template <class P, class W>
class Moebius_construct_regular_circum_center_2
{
 public:
  typedef CGAL::Moebius_point<P, W> Point;
  typedef /* blip blop */ P result_type;

  result_type
    operator() (const Point &p,
		const Point &q,
    		const Point &r,
    		const Point &s) const
    { 
      return moebius_regular_circumcenter_2 (p, q, r, s);
    }
};

template <class P, class W>
inline
P
moebius_regular_circumcenter_2 (const Moebius_point<P,W> &p,
			       const Moebius_point<P,W> &q,
			       const Moebius_point<P,W> &r,
			       const Moebius_point<P,W> &s,
			       Cartesian_tag)
{
  typedef typename P::R::FT FT;
  FT x, y;
  moebius_regular_circumcenterC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
				 q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
				 r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
				 s.x(), s.y(), FT (s.lambda()), FT (s.mu()),
				 x, y);
  return P (x, y);
}

template <class P, class W>
inline
P
moebius_regular_circumcenter_2 (const Moebius_point<P,W> &p,
			       const Moebius_point<P,W> &q,
			       const Moebius_point<P,W> &r,
			       const Moebius_point<P,W> &s,
			       Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  FT x, y;
  moebius_regular_circumcenterC2 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
				 q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
				 r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
				 s.x(), s.y(), FT (s.lambda()), FT (s.mu()),
				 x, y);
  return P (x, y);
}

template <class P, class W>
inline
P
moebius_regular_circumcenter_2 (const Moebius_point<P,W> &p,
			       const Moebius_point<P,W> &q,
			       const Moebius_point<P,W> &r,
			       const Moebius_point<P,W> &s)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_regular_circumcenter_2 (p, q, r, s, Tag ());
}

template <class K, class Weight_ = typename K::FT >
class Moebius_diagram_euclidean_traits_2
: public K
{
 public:
  typedef Weight_ Weight;

  typedef typename K::Point_2 Bare_point_2;
  typedef Moebius_point<Bare_point_2, Weight> Point_2;
  typedef CGAL::Conic_2<K> Conic_2;
  typedef CGAL::Segment_2<K> Segment_2;

  typedef CGAL::Moebius_regular_triangulation_euclidean_traits_3<K, Weight> Regular_traits;

  typedef CGAL::Compare_on_oriented_line_2<K> Compare_on_oriented_line_2;

  typedef CGAL::Moebius_orientation_2<Bare_point_2,Weight> Moebius_orientation_2;
  typedef CGAL::Moebius_has_circle_2<Bare_point_2,Weight> Moebius_has_circle_2;
  typedef CGAL::Moebius_circle_cross_line_2<Bare_point_2,Weight> Moebius_circle_cross_line_2;
  typedef CGAL::Moebius_circle_side_of_center_2<Bare_point_2,Weight> Moebius_circle_side_of_center_2;
  typedef CGAL::Moebius_side_of_vertex_2<Bare_point_2,Weight> Moebius_side_of_vertex_2;
  typedef CGAL::Moebius_circle_side_of_vertex_2<Bare_point_2,Weight> Moebius_circle_side_of_vertex_2;

  typedef CGAL::Moebius_construct_line_2<Bare_point_2,Weight> Moebius_construct_line_2;
  typedef CGAL::Moebius_construct_circle_2<Bare_point_2,Weight> Moebius_construct_circle_2;

  //typedef CGAL::Moebius_construct_regular_circum_center_2<Bare_point_2,Weight> Moebius_construct_regular_circum_center_2;
  typedef CGAL::Moebius_construct_circle_vertices_2<Bare_point_2,Weight,Segment_2> Moebius_construct_circle_vertices_2;
  typedef CGAL::Moebius_construct_line_vertices_2<Bare_point_2,Weight,Segment_2> Moebius_construct_line_vertices_2;
  typedef CGAL::Moebius_construct_vertex_2<Bare_point_2,Weight> Moebius_construct_vertex_2;

  Compare_on_oriented_line_2
    compare_on_oriented_line_2_object () const
    { return Compare_on_oriented_line_2 (); }

  Moebius_orientation_2
    moebius_orientation_2_object () const
    { return Moebius_orientation_2 (); }

  Moebius_has_circle_2
    has_circle_2_object () const
    { return Moebius_has_circle_2 (); }

  Moebius_circle_cross_line_2
    circle_cross_line_2_object () const
    { return Moebius_circle_cross_line_2 (); }

  Moebius_circle_side_of_center_2
    circle_side_of_center_2_object () const
    { return Moebius_circle_side_of_center_2 (); }

  Moebius_side_of_vertex_2
    side_of_vertex_2_object () const
    { return Moebius_side_of_vertex_2 (); }

  Moebius_circle_side_of_vertex_2
    circle_side_of_vertex_2_object () const
    { return Moebius_circle_side_of_vertex_2 (); }



  Moebius_construct_line_2
    construct_line_2_object () const
    { return Moebius_construct_line_2 (); }

  Moebius_construct_circle_2
    construct_circle_2_object () const
    { return Moebius_construct_circle_2 (); }

  //  Moebius_construct_regular_circum_center_2
  // construst_regular_circumcenter_2_object () const
  // { return Moebius_construct_regular_circum_center_2 (); }

  Moebius_construct_circle_vertices_2
    construct_circle_vertices_2_object () const
    { return Moebius_construct_circle_vertices_2 (); }

  Moebius_construct_line_vertices_2
    construct_line_vertices_2_object () const
    { return Moebius_construct_line_vertices_2 (); }

  Moebius_construct_vertex_2
    construct_vertex_2_object () const
    { return Moebius_construct_vertex_2 (); }
};




CGAL_END_NAMESPACE



#endif// CHR_MOEBIUS_DIAGRAM_EUCLIDEAN_TRAITS_2_H
