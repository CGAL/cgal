#ifndef CGAL_CONVEX_HULL_TRAITS_DUAL_3_H
#define CGAL_CONVEX_HULL_TRAITS_DUAL_3_H

// Dual predicates
#include <CGAL/Voronoi_covariance_3/predicates.h>

namespace CGAL
{
  namespace Voronoi_covariance_3
  {
    template < class R_ >
    class Point_triple_dual
    {
    protected:
      typedef typename R_::FT                   FT;
      typedef typename R_::Point_3              Point_3;
      typedef typename R_::Vector_3             Vector_3;

      Point_3  p_,  q_,  r_;
    public:
      typedef R_                                     R;

      Point_triple_dual() {}

      Point_triple_dual(const Point_3 &p, const Point_3 &q, const Point_3 &r)
	: p_(p), q_(q), r_(r)
      {}

      const Point_3& p() const { return p_; }
      const Point_3& q() const { return q_; }
      const Point_3& r() const { return r_; }

    };

    template <class K>
    class Point_triple_dual_has_on_positive_side_3 {

    public:
      typedef typename K::Point_3 Point_3;
      typedef typename K::Plane_3 Plane_3;
      bool
      operator()( const Plane_3& pl, const Point_3& p) const
      {
	typename K::Orientation_3 o;
	return ( o(pl.p(), pl.q(), pl.r(), p) == CGAL::POSITIVE );
      }
    };
    template <class K, class OldK>
    class Point_triple_dual_construct_orthogonal_vector_3
    {
    public:

      typedef typename K::Vector_3 Vector_3;
      typedef typename K::Plane_3 Plane_3;

      Vector_3 operator()(const Plane_3& plane) const
      {
	typename OldK::Construct_orthogonal_vector_3
	  construct_orthogonal_vector_3;
	return construct_orthogonal_vector_3(plane.p(), plane.q(), plane.r());
      }
    };


    template <class K>
    class Point_triple_dual_oriented_side_3
    {
    public:

      typedef typename K::Point_3 Point_3;
      typedef typename K::Plane_3 Plane_3;
      typedef Oriented_side    result_type;

      result_type
      operator()( const Plane_3& pl, const Point_3& p) const
      {
	typename K::Orientation_3 o;
	Orientation ori = o(pl.p(), pl.q(), pl.r(), p) ;
	if(ori > 0) return ON_POSITIVE_SIDE;
	if(ori < 0) return ON_NEGATIVE_SIDE;
	return ON_ORIENTED_BOUNDARY;
      }
    };

    template <typename K, typename OldK>
    class Point_triple_dual_less_signed_distance_to_plane_3
    {
    public:
      typedef typename K::Point_3 Point_3;
      typedef typename K::Plane_3 Plane_3;

      typedef bool             result_type;

      bool
      operator()( const Plane_3& h, const Point_3& p, const Point_3& q) const
      {
	const Point_3& hp = h.p();
	const Point_3& hq = h.q();
	const Point_3& hr = h.r();
	//typename OldK::Less_signed_distance_to_plane_3
	//	less_signed_distance_to_plane_3;
	// return less_signed_distance_to_plane_3(hp, hq, hr, p, q);
	return has_smaller_signed_dist_to_planeC3(hp.x(), hp.y(), hp.z(),
						  hq.x(), hq.y(), hq.z(),
						  hr.x(), hr.y(), hr.z(),
						  p.x(), p.y(), p.z(),
						  q.x(), q.y(), q.z());

      }
    };



    template <class T>
    class Max_coordinate_dual_3
    {
    public:

      int operator()(const T& v)
      {
	if (CGAL_NTS abs(v.x()) >= CGAL_NTS abs(v.y()))
	  {
	    if (CGAL_NTS abs(v.x()) >= CGAL_NTS abs(v.z())) return 0;
	    return 2;
	  }
	else
	  {
	    if (CGAL_NTS abs(v.y()) >= CGAL_NTS abs(v.z())) return 1;
	    return 2;
	  }
      }
    };

    template <typename GT>
    struct GT3_for_CH3_dual {
      typedef typename GT::Point_3 Point_2;
    };

    template <class R_>
    class Convex_hull_traits_dual_3
    {
    public:
      typedef R_                                     R;
      typedef Convex_hull_traits_dual_3<R>           Self;

      /* typedef typename R::Point_3                    Point_3; */
      /* typedef Point_triple_dual<R>                   Plane_3; */
      /* typedef typename R::Segment_3                  Segment_3; */
      /* typedef typename R::Triangle_3                 Triangle_3; */
      /* typedef typename R::Vector_3                   Vector_3; */

      // For ch_quickhull_polyhedron
      typedef typename R::Plane_3 Point_2;

      // Dual
      typedef typename R::Plane_3         Point_3;
      typedef Plane_dual<R>               Plane_3;
      typedef Segment_dual<R>             Segment_3;
      typedef Plane_dual<R>               Triangle_3;
      typedef typename R::Vector_3        Vector_3;

      // Construct objects
      class Construct_segment_3 {
      public:
	Segment_3 operator ()(const Point_3& p, const Point_3& q)
	{
	  return Segment_3(p, q);
	}
      };

      class Construct_triangle_3 {
      public:
	Triangle_3 operator ()(const Point_3& p, const Point_3& q, const Point_3& r)
	{
	  return Triangle_3(p, q, r);
	}
      };

      typedef typename R::Construct_vector_3         Construct_vector_3;
      typedef typename R::RT                         RT;

      class Construct_orthogonal_vector_3 {
      public:
	typedef typename R::Plane_3 Primal_plane_3;

	Vector_3 operator ()(const Plane_3& plane)
	{
	  Primal_plane_3 p1 = plane.p1;
	  Primal_plane_3 p2 = plane.p2;
	  Primal_plane_3 p3 = plane.p3;

	  RT alpha = (p1.d() * p2.b() - p2.d() * p1.b()) * (p1.d() * p3.c() - p3.d() * p1.c()) - (p1.d() * p2.d() - p2.d() * p1.c()) * (p1.d() * p3.b() - p3.d() * p1.b());
	  RT beta  = (p1.d() * p2.c() - p2.d() * p1.c()) * (p1.d() * p3.a() - p3.d() * p1.a()) - (p1.d() * p2.a() - p2.d() * p1.a()) * (p1.d() * p3.c() - p3.d() * p1.c());
	  RT gamma = (p1.d() * p2.a() - p2.d() * p1.a()) * (p1.d() * p3.b() - p3.d() * p1.b()) - (p1.d() * p2.b() - p2.d() * p1.b()) * (p1.d() * p3.a() - p3.d() * p1.a());

	  return Vector_3(alpha, beta, gamma);
	}
      };

      class Construct_plane_3 {
      public:
	Plane_3 operator ()(const Point_3& p, const Point_3& q, const Point_3& r)
	{
	  return Plane_3(p,q,r);
	}
      };

      /* typedef typename R::Construct_segment_3        Construct_segment_3; */
      /* typedef typename R::Construct_triangle_3       Construct_triangle_3; */
      /* typedef typename R::Construct_vector_3         Construct_vector_3; */

      // Predicates
      typedef Equal_3_dual_point<R>                         Equal_3;
      typedef Collinear_3_dual_point<R>                     Collinear_3;
      typedef Coplanar_3_dual_point<R>                      Coplanar_3;
      typedef Has_on_positive_side_3_dual_point<R>          Has_on_positive_side_3;
      typedef Less_distance_to_point_3_dual_point<R>        Less_distance_to_point_3;
      typedef Less_signed_distance_to_plane_3_dual_point<R> Less_signed_distance_to_plane_3;
      typedef Oriented_side_3_dual_point<R>                 Oriented_side_3;

      /* typedef typename R::Equal_3                    Equal_3; */
      /* typedef typename R::Collinear_3                Collinear_3; */
      /* typedef typename R::Coplanar_3                 Coplanar_3; */
      /* typedef Point_triple_dual_has_on_positive_side_3<Self>     Has_on_positive_side_3; */
      /* typedef typename R::Less_distance_to_point_3   Less_distance_to_point_3; */
      /* typedef  Point_triple_dual_less_signed_distance_to_plane_3<Self, R> */
      /*                                              Less_signed_distance_to_plane_3; */

      // Why is it needed ?
      /* typedef Point_triple_dual_construct_orthogonal_vector_3<Self, R> */
      /*                                                Construct_orthogonal_vector_3; */
      typedef typename R::Orientation_3              Orientation_3;

      // required for degenerate case of all points coplanar
      /* typedef CGAL::Max_coordinate_dual_3<Vector_3>             Max_coordinate_dual_3; */

      Construct_segment_3
      construct_segment_3_object() const
      { return Construct_segment_3(); }

      Construct_plane_3
      construct_plane_3_object() const
      { return Construct_plane_3(); }

      Construct_triangle_3
      construct_triangle_3_object() const
      { return Construct_triangle_3(); }

      Construct_vector_3
      construct_vector_3_object() const
      { return Construct_vector_3(); }

      Construct_orthogonal_vector_3
      construct_orthogonal_vector_3_object() const
      { return Construct_orthogonal_vector_3(); }

      Collinear_3
      collinear_3_object() const
      { return Collinear_3(); }

      Coplanar_3
      coplanar_3_object() const
      { return Coplanar_3(); }

      Less_distance_to_point_3
      less_distance_to_point_3_object() const
      { return Less_distance_to_point_3(); }

      Has_on_positive_side_3
      has_on_positive_side_3_object() const
      { return Has_on_positive_side_3(); }

      Equal_3
      equal_3_object() const
      { return Equal_3(); }

      Less_signed_distance_to_plane_3
      less_signed_distance_to_plane_3_object() const
      { return Less_signed_distance_to_plane_3(); }

      // TODO
      /* Max_coordinate_dual_3 */
      /* max_coordinate_3_object() const */
      /* { return Max_coordinate_dual_3(); } */
    };

  } // namespace CGAL
}

#endif // CGAL_CONVEX_HULL_TRAITS_DUAL_3_H
