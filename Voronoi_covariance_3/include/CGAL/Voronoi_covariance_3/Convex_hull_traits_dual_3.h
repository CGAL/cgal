#ifndef CGAL_CONVEX_HULL_TRAITS_DUAL_3_H
#define CGAL_CONVEX_HULL_TRAITS_DUAL_3_H

#include <CGAL/Voronoi_covariance_3/predicates.h>
#include <CGAL/Filtered_predicate.h>

namespace CGAL
{
  namespace Voronoi_covariance_3
  {
    // Base traits class for dual predicates
    template <class R_>
      class Convex_hull_traits_base_dual_3
      {
        public:
          typedef R_                                     R;
          typedef Convex_hull_traits_base_dual_3<R>      Self;

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
              Triangle_3 operator ()(const Point_3& p,
                                     const Point_3& q,
                                     const Point_3& r)
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

                RT alpha = (p1.d() * p2.b() - p2.d() * p1.b()) *
                    (p1.d() * p3.c() - p3.d() * p1.c()) -
                    (p1.d() * p2.d() - p2.d() * p1.c()) *
                    (p1.d() * p3.b() - p3.d() * p1.b());

                RT beta  = (p1.d() * p2.c() - p2.d() * p1.c()) *
                    (p1.d() * p3.a() - p3.d() * p1.a()) -
                    (p1.d() * p2.a() - p2.d() * p1.a()) *
                    (p1.d() * p3.c() - p3.d() * p1.c());

                RT gamma = (p1.d() * p2.a() - p2.d() * p1.a()) *
                    (p1.d() * p3.b() - p3.d() * p1.b()) -
                    (p1.d() * p2.b() - p2.d() * p1.b()) *
                    (p1.d() * p3.a() - p3.d() * p1.a());

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
      };

    // Non-filtered traits class
    template <class R_, bool Has_filtered_predicates = R_::Has_filtered_predicates >
        class Convex_hull_traits_dual_3
            : public Convex_hull_traits_base_dual_3<R_>
    {} ;

    // Converter for dual planes
    template <class K1, class K2>
        struct Cartesian_converter_dual : public CGAL::Cartesian_converter<K1, K2>
    {
        using CGAL::Cartesian_converter<K1, K2>::operator();

        Plane_dual<K2> operator() (const Plane_dual<K1> &in) const
        {
            return Plane_dual<K2>(operator()(in.p1),
                                  operator()(in.p2),
                                  operator()(in.p3));
        }
    };

    // Filtered traits
    template <typename R_>
        class Convex_hull_filtered_traits_dual_3
        : public Convex_hull_traits_base_dual_3<R_>
        {
            public:
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

                Collinear_3 collinear_3_object() const
                    { return Collinear_3(); }

                Coplanar_3 coplanar_3_object() const
                    { return Coplanar_3(); }

                Less_distance_to_point_3 less_distance_to_point_3_object() const
                    { return Less_distance_to_point_3(); }

                Equal_3 equal_3_object() const
                    { return Equal_3(); }

                Has_on_positive_side_3 has_on_positive_side_3_object() const
                    { return Has_on_positive_side_3(); }

                Less_signed_distance_to_plane_3
                    less_signed_distance_to_plane_3_object() const
                    { return Less_signed_distance_to_plane_3(); }

                // Constructions are inherited
        };

    // Traits specialization
    template <typename R_>
        class Convex_hull_traits_dual_3<R_, true>
        : public Convex_hull_filtered_traits_dual_3<R_>
        {} ;
  } // namespace Voronoi_covariance_3
} // namespace CGAL

#endif // CGAL_CONVEX_HULL_TRAITS_DUAL_3_H

