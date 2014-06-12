#ifndef CGAL_CONVEX_HULL_TRAITS_DUAL_3_H
#define CGAL_CONVEX_HULL_TRAITS_DUAL_3_H

#include <CGAL/Voronoi_covariance_3/predicates.h>
#include <CGAL/Filtered_predicate.h>
#include <CGAL/Cartesian_converter.h>

namespace CGAL
{
  namespace Voronoi_covariance_3
  {
    // Base traits class for dual predicates
    template <class R_>
      class Convex_hull_traits_base_dual_3
      {
        private:
          // Origin
          typedef typename R_::Point_3 Primal_point_3;
          Primal_point_3 origin;

        public:

          Convex_hull_traits_base_dual_3 (Primal_point_3 o =
                                          Primal_point_3(0, 0, 0)) : origin(o)
          {}

          typedef R_                                     R;
          typedef Convex_hull_traits_base_dual_3<R>      Self;

          // Dual
          typedef typename R::Plane_3         Point_3;
          typedef Plane_dual<R>               Plane_3;
          typedef Segment_dual<R>             Segment_3;
          typedef Plane_dual<R>               Triangle_3;
          typedef Vector_dual<R>              Vector_3;

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
              Triangle_3 operator ()(const Point_3& p,
                                     const Point_3& q,
                                     const Point_3& r)
              {
                return Triangle_3(p, q, r);
              }
          };

          class Construct_vector_3 {
            public:
                Vector_3 operator ()(const Point_3& p,
                                     const Point_3& q)
              {
                return Vector_3(p, q);
              }

              Vector_3 operator ()(int x,
                                   int y,
                                   int z)
              {
                  // TODO
                  Point_3 p(1, 1, 1, 1);
                  Point_3 q(1, 1, 1, 1);

                  return Vector_3(p, q);
              }
          };

          typedef typename R::RT                         RT;

          class Construct_orthogonal_vector_3 {
              private:
                  // Origin
                  typedef typename R_::Point_3 Primal_point_3;
                  Primal_point_3 origin;

              public:
              typedef typename R::Plane_3 Primal_plane_3;

              Construct_orthogonal_vector_3 (Primal_point_3 o =
                                             Primal_point_3(0, 0, 0)) : origin(o)
              {}

              Vector_3 operator ()(const Plane_3& plane)
              {
                Primal_plane_3 p1 = plane.p1;
                Primal_plane_3 p2 = plane.p2;
                Primal_plane_3 p3 = plane.p3;

                RT dp1 = p1.d() + origin.x() * p1.a()
                    + origin.y() * p1.b() + origin.z() * p1.c();
                RT dp2 = p2.d() + origin.x() * p2.a()
                    + origin.y() * p2.b() + origin.z() * p2.c();
                RT dp3 = p3.d() + origin.x() * p3.a()
                    + origin.y() * p3.b() + origin.z() * p3.c();

                // Normal to the dual plane
                RT alpha = (dp1 * p2.b() - dp2 * p1.b()) *
                    (dp1 * p3.c() - dp3 * p1.c()) -
                    (dp1 * p2.c() - dp2 * p1.c()) *
                    (dp1 * p3.b() - dp3 * p1.b());

                RT beta  = (dp1 * p2.c() - dp2 * p1.c()) *
                    (dp1 * p3.a() - dp3 * p1.a()) -
                    (dp1 * p2.a() - dp2 * p1.a()) *
                    (dp1 * p3.c() - dp3 * p1.c());

                RT gamma = (dp1 * p2.a() - dp2 * p1.a()) *
                    (dp1 * p3.b() - dp3 * p1.b()) -
                    (dp1 * p2.b() - dp2 * p1.b()) *
                    (dp1 * p3.a() - dp3 * p1.a());

                return Vector_3(alpha, beta, gamma);
              }
          };

          class Construct_plane_3 {
              public:
                  Plane_3 operator ()(const Point_3& p,
                                      const Point_3& q,
                                      const Point_3& r)
                  {
                      return Plane_3(p,q,r);
                  }
          };

          // Predicates
          typedef Equal_3_dual_point<R>                         Equal_3;
          typedef Collinear_3_dual_point<R>                     Collinear_3;
          typedef Coplanar_3_dual_point<R>                      Coplanar_3;
          typedef Has_on_positive_side_3_dual_point<R>          Has_on_positive_side_3;
          typedef Less_distance_to_point_3_dual_point<R>        Less_distance_to_point_3;
          typedef Less_signed_distance_to_plane_3_dual_point<R> Less_signed_distance_to_plane_3;
          typedef Orientation_3_dual_point<R> Orientation_3;

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
              { return Construct_orthogonal_vector_3(origin); }

          Collinear_3
              collinear_3_object() const
              { return Collinear_3(origin); }

          Coplanar_3
              coplanar_3_object() const
              { return Coplanar_3(origin); }

          Less_distance_to_point_3
              less_distance_to_point_3_object() const
              { return Less_distance_to_point_3(origin); }

          Has_on_positive_side_3
              has_on_positive_side_3_object() const
              { return Has_on_positive_side_3(origin); }

          Equal_3
              equal_3_object() const
              { return Equal_3(origin); }

          Less_signed_distance_to_plane_3
              less_signed_distance_to_plane_3_object() const
              { return Less_signed_distance_to_plane_3(origin); }

          Orientation_3
              orientation_3_object() const
              { return Orientation_3(origin); }
      };

    // Non-filtered traits class
    template <class R_, bool Has_filtered_predicates = R_::Has_filtered_predicates >
        class Convex_hull_traits_dual_3
        : public Convex_hull_traits_base_dual_3<R_>
        {} ;

    // Converter for dual planes
    template <class K1, class K2>
        struct Cartesian_converter_dual : public Cartesian_converter<K1, K2>
    {
        using CGAL::Cartesian_converter<K1, K2>::operator();

        Plane_dual<K2> operator() (const Plane_dual<K1> &in) const
        {
            return Plane_dual<K2>(operator()(in.p1),
                                  operator()(in.p2),
                                  operator()(in.p3));
        }

        Vector_dual<K2> operator() (const Vector_dual<K1> &in) const
        {
            return Vector_dual<K2>(operator()(in.p),
                                   operator()(in.q));
        }
    };

    // Filtered traits
    template <typename R_>
        class Convex_hull_filtered_traits_dual_3
        : public Convex_hull_traits_base_dual_3<R_>
        {
            private:
                // Origin
                typedef typename R_::Point_3 Primal_point_3;
                Primal_point_3 origin;

            public:
                Convex_hull_filtered_traits_dual_3 (Primal_point_3 o =
                                                    Primal_point_3(0, 0, 0)) : origin(o)
                {}

                // Exact traits is based on the exact kernel.
                typedef Convex_hull_traits_dual_3<typename R_::Exact_kernel_rt>
                    Exact_traits;

                // Filtering traits is based on the filtering kernel.
                typedef Convex_hull_traits_dual_3<typename R_::Approximate_kernel>
                    Filtering_traits;

                // Converters
                typedef Cartesian_converter_dual<R_, typename R_::Exact_kernel_rt> Converter_exact_dual;
                typedef Cartesian_converter_dual<R_, typename R_::Approximate_kernel> Converter_approx_dual;

                // Filtered predicates
                typedef Filtered_predicate<
                    typename Exact_traits::Equal_3,
                             typename Filtering_traits::Equal_3,
                             Converter_exact_dual ,
                             Converter_approx_dual > Equal_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Collinear_3,
                             typename Filtering_traits::Collinear_3,
                             Converter_exact_dual,
                             Converter_approx_dual > Collinear_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Coplanar_3,
                             typename Filtering_traits::Coplanar_3,
                             Converter_exact_dual,
                             Converter_approx_dual > Coplanar_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Less_distance_to_point_3,
                             typename Filtering_traits::Less_distance_to_point_3,
                             Converter_exact_dual,
                             Converter_approx_dual > Less_distance_to_point_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Has_on_positive_side_3,
                             typename Filtering_traits::Has_on_positive_side_3,
                             Converter_exact_dual,
                             Converter_approx_dual > Has_on_positive_side_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Less_signed_distance_to_plane_3,
                             typename Filtering_traits::Less_signed_distance_to_plane_3,
                             Converter_exact_dual,
                             Converter_approx_dual > Less_signed_distance_to_plane_3;

                typedef Filtered_predicate<
                    typename Exact_traits::Orientation_3,
                             typename Filtering_traits::Orientation_3,
                             Converter_exact_dual,
                             Converter_approx_dual > Orientation_3;

                Collinear_3 collinear_3_object() const
                { return Collinear_3(origin); }

                Coplanar_3 coplanar_3_object() const
                { return Coplanar_3(origin); }

                Less_distance_to_point_3 less_distance_to_point_3_object() const
                { return Less_distance_to_point_3(origin); }

                Equal_3 equal_3_object() const
                { return Equal_3(origin); }

                Has_on_positive_side_3 has_on_positive_side_3_object() const
                { return Has_on_positive_side_3(origin); }

                Less_signed_distance_to_plane_3
                less_signed_distance_to_plane_3_object() const
                { return Less_signed_distance_to_plane_3(origin); }

                Orientation_3
                    orientation_3_object() const
                    { return Orientation_3(origin); }

                // Constructions are inherited
        };

    // Traits specialization
    template <typename R_>
        class Convex_hull_traits_dual_3<R_, true>
        : public Convex_hull_filtered_traits_dual_3<R_>
        {
            private:
                typedef typename R_::Point_3 Primal_point_3;

            public:
                Convex_hull_traits_dual_3 (Primal_point_3 o =
                                           Primal_point_3(0, 0, 0)) : Convex_hull_filtered_traits_dual_3<R_>(o)
                {}
        } ;
  } // namespace Voronoi_covariance_3
} // namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_DUAL_3_H

