#ifndef CHR_MOEBIUS_REGULAR_EUCLIDEAN_TRAITS_3_H
#define CHR_MOEBIUS_REGULAR_EUCLIDEAN_TRAITS_3_H

#include <CGAL/basic.h>
#include <CGAL/Moebius_point.h>
#include <CGAL/Regular_triangulation_euclidean_traits_3.h>
#include <CGAL/predicates/Moebius_diagram_regularftC2.h>

CGAL_BEGIN_NAMESPACE

template <class P, class W>
inline
Orientation
moebius_orientation_3 (const Moebius_point<P, W> &p,
		      const Moebius_point<P, W> &q,
		      const Moebius_point<P, W> &r,
		      const Moebius_point<P, W> &s,
		      Cartesian_tag)
{
  typedef typename P::R::FT FT;
   return moebius_orientationC3 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
				q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
				r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
				s.x(), s.y(), FT (s.lambda()), FT (s.mu()));
}


template <class P, class W>
inline
Orientation
moebius_orientation_3 (const Moebius_point<P, W> &p,
		      const Moebius_point<P, W> &q,
		      const Moebius_point<P, W> &r,
		      const Moebius_point<P, W> &s,
		      Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  return moebius_orientationC3 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			       q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			       r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
			       s.x(), s.y(), FT (s.lambda()), FT (s.mu()));
}


template <class P, class W>
inline
Orientation
moebius_orientation_3 (const Moebius_point<P, W> &p,
		      const Moebius_point<P, W> &q,
		      const Moebius_point<P, W> &r,
		      const Moebius_point<P, W> &s)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_orientation_3 (p, q, r, s, Tag ());
}


template <class Point, class Weight>
class Moebius_orientation_3
{
public:
  typedef Orientation  result_type;
  typedef CGAL::Moebius_point <Point, Weight>        Weighted_point;

  Orientation operator() ( const Weighted_point & p,
			   const Weighted_point & q,
			   const Weighted_point & r,
			   const Weighted_point & s)
    {
      //std::cerr << "orientation "
      //	<< "(" << p << ")"
      //	<< "(" << q << ")"
      //	<< "(" << r << ")"
      //	<< "(" << s << ")\n";
      return moebius_orientation_3 (p, q, r, s);
    }
};


// ------------------------

template <class P, class W>
inline
Orientation
moebius_coplanar_orientation_3 (const Moebius_point<P, W> &p,
			       const Moebius_point<P, W> &q,
			       const Moebius_point<P, W> &r,
			       Cartesian_tag)
{
  typedef typename P::R::FT FT;
   return moebius_coplanar_orientationC3 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
					 q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
					 r.x(), r.y(), FT (r.lambda()), FT (r.mu()));
}


template <class P, class W>
inline
Orientation
moebius_coplanar_orientation_3 (const Moebius_point<P, W> &p,
		      const Moebius_point<P, W> &q,
		      const Moebius_point<P, W> &r,
		      Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  return moebius_coplanar_orientationC3 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
					q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
					r.x(), r.y(), FT (r.lambda()), FT (r.mu()));
}


template <class P, class W>
inline
Orientation
moebius_coplanar_orientation_3 (const Moebius_point<P, W> &p,
			       const Moebius_point<P, W> &q,
			       const Moebius_point<P, W> &r)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_coplanar_orientation_3 (p, q, r, Tag ());
}


template <class Point, class Weight>
class Moebius_coplanar_orientation_3
{
public:
  typedef Orientation  result_type;
  typedef CGAL::Moebius_point <Point, Weight>        Weighted_point;

  Orientation operator() ( const Weighted_point & p,
			   const Weighted_point & q,
			   const Weighted_point & r)
    {
      return moebius_coplanar_orientation_3 (p, q, r);
    }
};

// ------------------------

template <class P, class W>
inline
Comparison_result
moebius_compare_xyz_3 (const Moebius_point<P, W> &p,
		      const Moebius_point<P, W> &q,
		      Cartesian_tag)
{
  typedef typename P::R::FT FT;
   return moebius_compare_xyzC3 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
				q.x(), q.y(), FT (q.lambda()), FT (q.mu()));
}


template <class P, class W>
inline
Comparison_result
moebius_compare_xyz_3 (const Moebius_point<P, W> &p,
		      const Moebius_point<P, W> &q,
		      Homogeneous_tag)
{
  typedef typename P::R::FT FT;
  return moebius_compare_xyzC3 (p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			       q.x(), q.y(), FT (q.lambda()), FT (q.mu()));
}


template <class P, class W>
inline
Comparison_result
moebius_compare_xyz_3 (const Moebius_point<P, W> &p,
		      const Moebius_point<P, W> &q)
{
  typedef typename P::R::Rep_tag Tag;
  return moebius_compare_xyz_3 (p, q, Tag ());
}


template <class Point, class Weight>
class Moebius_compare_xyz_3
{
public:
  typedef Comparison_result  result_type;
  typedef CGAL::Moebius_point <Point, Weight>        Weighted_point;

  Comparison_result operator() ( const Weighted_point & p,
				 const Weighted_point & q)
    {
      return moebius_compare_xyz_3 (p, q);
    }
};

// ------------------------


template <class Point, class Weight>
class Moebius_power_test_3
{
public:
  typedef Oriented_side  result_type;
  typedef CGAL::Moebius_point <Point, Weight>        Weighted_point;

  Oriented_side operator() ( const Weighted_point & p,
			     const Weighted_point & q,
			     const Weighted_point & r,
			     const Weighted_point & s,
			     const Weighted_point & t) const
    {
      return CGAL::moebius_power_test(p,q,r,s,t);
    }
  Oriented_side operator() ( const Weighted_point & p,
			     const Weighted_point & q,
			     const Weighted_point & r,
			     const Weighted_point & s) const
    {
      return CGAL::moebius_power_test(p,q,r,s);
    }

  Oriented_side operator() ( const Weighted_point & p,
			     const Weighted_point & q,
			     const Weighted_point & r) const
    {
      return CGAL::moebius_power_test(p,q,r);
    }

  Oriented_side operator() ( const Weighted_point & p,
			     const Weighted_point & q) const
    {
      return CGAL::moebius_power_test(p,q);
    }
};

template <class K, class Weight = CGAL_TYPENAME_MSVC_NULL K::FT>
class Moebius_regular_triangulation_euclidean_traits_3 : public K
{
 public:
  typedef typename K::Point_2                        Bare_point;
  typedef CGAL::Moebius_point<Bare_point, Weight>     Weighted_point;
  typedef Weighted_point                             Point_3;

  // The next typedef is there for backward compatibility
  // Some users take their point type from the traits class.
  // Before this type was Point
  typedef Point_3 Point;

  typedef CGAL::Moebius_power_test_3<Bare_point, Weight> Power_test_3;
  typedef CGAL::Moebius_orientation_3<Bare_point, Weight> Orientation_3;
  typedef CGAL::Moebius_coplanar_orientation_3<Bare_point, Weight> Coplanar_orientation_3;
  typedef CGAL::Moebius_compare_xyz_3<Bare_point, Weight> Compare_xyz_3;

  Power_test_3 
  power_test_3_object() const
    {  return Power_test_3(); }

  Orientation_3
    orientation_3_object() const
    { return Orientation_3 (); }

  Coplanar_orientation_3
    coplanar_orientation_3_object() const
    { return Coplanar_orientation_3 (); }

  Compare_xyz_3
    compare_xyz_3_object() const
    { return Compare_xyz_3 (); }
};


// regular triangulation predicates
// Cartesian versions.
template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt, Weight> &p,
		  const Moebius_point<pt, Weight> &q,
		  const Moebius_point<pt, Weight> &r,
		  const Moebius_point<pt, Weight> &s,
		  const Moebius_point<pt, Weight> &t,
		  Cartesian_tag)
{
  typedef typename pt::R::FT FT;
  return moebius_power_testC3(p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			     q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			     r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
			     s.x(), s.y(), FT (s.lambda()), FT (s.mu()),
			     t.x(), t.y(), FT (t.lambda()), FT (t.mu()));
}

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt, Weight> &p,
		  const Moebius_point<pt, Weight> &q,
		  const Moebius_point<pt, Weight> &r,
		  const Moebius_point<pt, Weight> &t,
		  Cartesian_tag)
{
  typedef typename pt::R::FT FT;
  return moebius_power_testC3(p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			     q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			     r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
			     t.x(), t.y(), FT (t.lambda()), FT (t.mu()));
}

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt, Weight> &p,
		  const Moebius_point<pt, Weight> &q,
		  const Moebius_point<pt, Weight> &t,
		  Cartesian_tag)
{
  typedef typename pt::R::FT FT;
  return moebius_power_testC3(p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			     q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			     t.x(), t.y(), FT (t.lambda()), FT (t.mu()));
}

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt, Weight> &p,
		  const Moebius_point<pt, Weight> &q,
		  Cartesian_tag)
{
  typedef typename pt::R::FT FT;
  return moebius_power_testC3(p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			     q.x(), q.y(), FT (p.lambda()), FT (q.mu()));
}


// Homogeneous versions.
template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt, Weight> &p,
		  const Moebius_point<pt, Weight> &q,
		  const Moebius_point<pt, Weight> &r,
		  const Moebius_point<pt, Weight> &s,
		  const Moebius_point<pt, Weight> &t,
		  Homogeneous_tag)
{
  typedef typename pt::R::RT RT;
  return moebius_power_testH3(p.hx(), p.hy(), p.hw(), RT (p.lambda()), RT(p.mu()),
			     q.hx(), q.hy(), q.hw(), RT (q.lambda()), RT(q.mu()),
			     r.hx(), r.hy(), r.hw(), RT (r.lambda()), RT(r.mu()),
			     s.hx(), s.hy(), s.hw(), RT (s.lambda()), RT(s.mu()),
			     t.hx(), t.hy(), t.hw(), RT (t.lambda()), RT(t.mu()));
}

// The 2 following call the cartesian version over FT, because an homogeneous
// special version would be boring to write.

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt, Weight> &p,
		  const Moebius_point<pt, Weight> &q,
		  const Moebius_point<pt, Weight> &r,
		  const Moebius_point<pt, Weight> &t,
		  Homogeneous_tag)
{
  typedef typename pt::R::FT FT;
  return moebius_power_testC3(p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			     q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			     r.x(), r.y(), FT (r.lambda()), FT (r.mu()),
			     t.x(), t.y(), FT (t.lambda()), FT (t.mu()));
}

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt, Weight> &p,
		  const Moebius_point<pt, Weight> &q,
		  const Moebius_point<pt, Weight> &t,
		  Homogeneous_tag)
{
  typedef typename pt::R::FT FT;
  return moebius_power_testC3(p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			     q.x(), q.y(), FT (q.lambda()), FT (q.mu()),
			     t.x(), t.y(), FT (t.lambda()), FT (t.mu()));
}

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt, Weight> &p,
		  const Moebius_point<pt, Weight> &q,
		  Homogeneous_tag)
{
  typedef typename pt::R::FT FT;
  return moebius_power_testC3(p.x(), p.y(), FT (p.lambda()), FT (p.mu()),
			     q.x(), q.y(), FT (q.lambda()), FT (q.mu()));
}


// Kludges for M$.

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt,Weight> &p,
		  const Moebius_point<pt,Weight> &q,
		  const Moebius_point<pt,Weight> &r,
		  const Moebius_point<pt,Weight> &s,
		  const Moebius_point<pt,Weight> &t)
{
  typedef typename pt::R::Rep_tag Tag;
  return moebius_power_test(p,q,r,s,t, Tag());
}

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt,Weight> &p,
		  const Moebius_point<pt,Weight> &q,
		  const Moebius_point<pt,Weight> &r,
		  const Moebius_point<pt,Weight> &t)
{
  typedef typename pt::R::Rep_tag Tag;
  return moebius_power_test(p,q,r,t, Tag());
}

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt,Weight> &p,
		  const Moebius_point<pt,Weight> &q,
		  const Moebius_point<pt,Weight> &t)
{
  typedef typename pt::R::Rep_tag Tag;
  return moebius_power_test(p,q,t, Tag());
}

template < class pt, class Weight >
inline
Oriented_side
moebius_power_test(const Moebius_point<pt,Weight> &p,
		  const Moebius_point<pt,Weight> &q)
{
  typedef typename pt::R::Rep_tag Tag;
  return moebius_power_test(p,q, Tag());
}




CGAL_END_NAMESPACE
#endif //  CHR_MOEBIUS_REGULAR_EUCLIDEAN_TRAITS_3_H
