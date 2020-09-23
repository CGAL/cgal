// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stephane Tayeb
//
//******************************************************************************
// File Description :
//******************************************************************************

#ifndef CGAL_MESH_3_ROBUST_INTERSECTION_TRAITS_3_H
#define CGAL_MESH_3_ROBUST_INTERSECTION_TRAITS_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/Mesh_3/config.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Mesh_3/Profile_counter.h>
#include <CGAL/internal/Static_filters/tools.h>


namespace CGAL {

namespace Mesh_3 {

template <typename K>
struct Vector_plane_orientation_3 {
  typedef typename K::Point_3 Point_3;

  typedef CGAL::Orientation result_type;

  Orientation
  operator()(const Point_3& p, const Point_3& q,
             const Point_3& a, const Point_3& b, const Point_3& c) const
  {
    typedef typename K::Vector_3 Vector;

    Vector ab = b - a;
    Vector ac = c - a;
    Vector pq = q - p;

    return K().orientation_3_object()(ab, ac, pq);
  }
}; // end struct template Vector_plane_orientation_3

template <typename K>
struct Vector_plane_orientation_3_filtered_predicate
  : public CGAL::Filtered_predicate<
    Vector_plane_orientation_3<CGAL::Exact_predicates_exact_constructions_kernel>,
    Vector_plane_orientation_3<Simple_cartesian<Interval_nt_advanced> >,
    Cartesian_converter<K, CGAL::Exact_predicates_exact_constructions_kernel>,
    Cartesian_converter<K, Simple_cartesian<Interval_nt_advanced> >
  >
{}; // end struct template Vector_plane_orientation_3_filtered_predicate

template <typename K>
struct Vector_plane_orientation_3_static_filter :
    public Vector_plane_orientation_3_filtered_predicate<K>
{
  typedef typename K::Point_3 Point_3;
  typedef Vector_plane_orientation_3_filtered_predicate<K> Base;

  Orientation
  operator()(const Point_3& p, const Point_3& q,
             const Point_3& a, const Point_3& b, const Point_3& c) const
  {
    CGAL_BRANCH_PROFILER(std::string("semi-static failures/calls to   : ")
                         + std::string(CGAL_PRETTY_FUNCTION), tmp);

    CGAL::internal::Static_filters_predicates::Get_approx<Point_3> get_approx;
    // Identity functor for all points but lazy points.

    double px, py, pz, qx, qy, qz;
    double ax, ay, az, bx, by, bz, cx, cy, cz;

    using CGAL::internal::fit_in_double;

    if (fit_in_double(get_approx(p).x(), px) && fit_in_double(get_approx(p).y(), py) &&
        fit_in_double(get_approx(p).z(), pz) &&
        fit_in_double(get_approx(q).x(), qx) && fit_in_double(get_approx(q).y(), qy) &&
        fit_in_double(get_approx(q).z(), qz) &&
        fit_in_double(get_approx(a).x(), ax) && fit_in_double(get_approx(a).y(), ay) &&
        fit_in_double(get_approx(a).z(), az) &&
        fit_in_double(get_approx(b).x(), bx) && fit_in_double(get_approx(b).y(), by) &&
        fit_in_double(get_approx(b).z(), bz) &&
        fit_in_double(get_approx(c).x(), cx) && fit_in_double(get_approx(c).y(), cy) &&
        fit_in_double(get_approx(c).z(), cz))
    { // This bloc is not indented because it was added in a second step,
      // and one wants to avoid the reindentation of the whole code

    double abx = bx - ax;
    double acx = cx - ax;
    double pqx = qx - px;
    double aby = by - ay;
    double acy = cy - ay;
    double pqy = qy - py;
    double abz = bz - az;
    double acz = cz - az;
    double pqz = qz - pz;

    double maxx = CGAL::abs(abx);
    double maxy = CGAL::abs(aby);
    double maxz = CGAL::abs(abz);

    double aprx = CGAL::abs(acx);
    double apsx = CGAL::abs(pqx);

    double apry = CGAL::abs(acy);
    double apsy = CGAL::abs(pqy);

    double aprz = CGAL::abs(acz);
    double apsz = CGAL::abs(pqz);

#ifdef CGAL_USE_SSE2_MAX
    CGAL::Max<double> mmax;

    maxx = mmax(maxx, aprx, apsx);
    maxy = mmax(maxy, apry, apsy);
    maxz = mmax(maxz, aprz, apsz);
#else
    if (maxx < aprx) maxx = aprx;
    if (maxx < apsx) maxx = apsx;
    if (maxy < apry) maxy = apry;
    if (maxy < apsy) maxy = apsy;
    if (maxz < aprz) maxz = aprz;
    if (maxz < apsz) maxz = apsz;
#endif
    double det = CGAL::determinant(abx, aby, abz,
                                   acx, acy, acz,
                                   pqx, pqy, pqz);

    double eps = 5.1107127829973299e-15 * maxx * maxy * maxz;

#ifdef CGAL_USE_SSE2_MAX
#if 0
    CGAL::Min<double> mmin;
    double tmp = mmin(maxx, maxy, maxz);
    maxz = mmax(maxx, maxy, maxz);
    maxx = tmp;
#else
    sse2minmax(maxx,maxy,maxz);
    // maxy can contain ANY element
#endif
#else
    // Sort maxx < maxy < maxz.
    if (maxx > maxz)
      std::swap(maxx, maxz);
    if (maxy > maxz)
      std::swap(maxy, maxz);
    else if (maxy < maxx)
      std::swap(maxx, maxy);
#endif
    // Protect against underflow in the computation of eps.
    if (maxx < 1e-97) /* cbrt(min_double/eps) */ {
      if (maxx == 0)
        return ZERO;
    }
    // Protect against overflow in the computation of det.
    else if (maxz < 1e102) /* cbrt(max_double [hadamard]/4) */ {

      if (det > eps)  return POSITIVE;
      if (det < -eps) return NEGATIVE;
    }

    } // end of the unindented block
    CGAL_BRANCH_PROFILER_BRANCH(tmp);

    return Base::operator()(p, q, a, b, c);
  }
}; // end struct template Vector_plane_orientation_3_static_filter

// An exact line-plane intersection that returns always a point.
// The line it (pq), and the plane is (abc).
//
// Precondition: [pq] and [abc] are not coplanar.
template <typename K>
typename K::Point_3
lp_intersection(const typename K::Point_3& p, const typename K::Point_3& q,
                const typename K::Point_3& a, const typename K::Point_3& b,
                const typename K::Point_3& c,
                const K &)
{
  CGAL_MESH_3_PROFILER(std::string(CGAL_PRETTY_FUNCTION));
  typedef Exact_predicates_exact_constructions_kernel   EPEC;
  typedef Simple_cartesian<CGAL::internal::Exact_field_selector<double>::Type> EK;
  typedef Simple_cartesian<Interval_nt_advanced>        FK; // filtering kernel
  typedef Cartesian_converter<typename K::Kernel, EK>    To_exact;
  typedef Cartesian_converter<typename K::Kernel, FK>    To_filtering;
  typedef Cartesian_converter<EK, typename K::Kernel>    Back_from_exact;

  To_filtering to_filtering;

  CGAL_BRANCH_PROFILER(std::string(" failures/calls to   : ") + std::string(CGAL_PRETTY_FUNCTION), tmp);

  { // First, one computes the intersection using a filtering kernel
    // (intervals as number type).
    Protect_FPU_rounding<true> protect;
    const typename FK::Point_3 fp = to_filtering(p), fq = to_filtering(q);
    const typename FK::Point_3
      fa = to_filtering(a), fb = to_filtering(b), fc = to_filtering(c);

    const typename FK::Vector_3 ap = fp - fa;
    const typename FK::Vector_3 ab = fb - fa;
    const typename FK::Vector_3 ac = fc - fa;
    const typename FK::Vector_3 qp = fq - fp;
    const typename FK::Vector_3 abac = cross_product(ab, ac);
    const FK::FT den = (qp * abac);
    // [Note: den is equal to det(qp, ab, ac).]

    // If the denominator is certainly not 0...
    if(is_certain(den != 0)) {
      // ...then compute the intersection point...
      const typename FK::Point_3 result = fp - qp * ( (ap * abac) / den );
      //                       [Note: ap * abac is equal to det(ap, ab, ac).]

      // ...and check that its coordinates intervals are precise enough.
      const double prec = EPEC::FT::get_relative_precision_of_to_double();
      if(has_smaller_relative_precision(result.x(),prec) &&
         has_smaller_relative_precision(result.y(),prec) &&
         has_smaller_relative_precision(result.z(),prec))
      {
        // Then return the intersection point as a K::Point_3.
        CGAL::FPU_set_cw(CGAL_FE_TONEAREST);
        return typename K::Point_3(to_double(result.x()),
                                   to_double(result.y()),
                                   to_double(result.z()));
      }
    }
  }

  CGAL_BRANCH_PROFILER_BRANCH(tmp);

  // In case the filtering computation was not successful (not 'return' so
  // far), then recompute the same expression using an exact kernel...

  To_exact to_exact;
  Back_from_exact from_exact;

  const typename EK::Point_3 ep = to_exact(p), eq = to_exact(q);
  const typename EK::Point_3 ea = to_exact(a), eb = to_exact(b), ec = to_exact(c);
  const typename EK::Vector_3 ap = ep - ea;
  const typename EK::Vector_3 ab = eb - ea;
  const typename EK::Vector_3 ac = ec - ea;
  const typename EK::Vector_3 qp = eq - ep;
  const typename EK::Vector_3 abac = cross_product(ab, ac);

  // .... and return the result converted to K::Point_3.
  return from_exact(ep - qp * ( (ap * abac) / (qp * abac) ));
}

// And intersection function for Triangle_3 and Segment_3 that always
// returns a point or the empty Object. In case of degeneracy, the empty
// Object is returned as well.
template <class K>
typename cpp11::result_of<
  typename K::Intersect_3(typename K::Segment_3, typename K::Triangle_3)>::type
ts_intersection(const typename K::Triangle_3 &t,
                const typename K::Segment_3  &s,
                const K & k)
{
  typedef typename cpp11::result_of<
    typename K::Intersect_3(typename K::Segment_3, typename K::Triangle_3)
  >::type result_type;

  CGAL_MESH_3_BRANCH_PROFILER(std::string("coplanar/calls in   : ") +
                              std::string(CGAL_PRETTY_FUNCTION), tmp);
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t) ) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(s) ) ;

  typedef typename K::Point_3 Point_3;

  typename K::Construct_point_on_3 point_on =
    k.construct_point_on_3_object();

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Orientation_3 orientation =
    k.orientation_3_object();

  const Point_3 & a = vertex_on(t,0);
  const Point_3 & b = vertex_on(t,1);
  const Point_3 & c = vertex_on(t,2);
  const Point_3 & p = point_on(s,0);
  const Point_3 & q = point_on(s,1);

  const Orientation abcp = orientation(a,b,c,p);
  const Orientation abcq = orientation(a,b,c,q);

  switch ( abcp ) {
    case POSITIVE:
      switch ( abcq ) {
        case POSITIVE:
          // the segment lies in the positive open halfspaces defined by the
          // triangle's supporting plane
          return result_type();

        case NEGATIVE:
          // p sees the triangle in counterclockwise order
          if ( orientation(p,q,a,b) != POSITIVE
              && orientation(p,q,b,c) != POSITIVE
              && orientation(p,q,c,a) != POSITIVE )
          {
            // The intersection is a point
            return result_type( lp_intersection(p, q, a, b, c, k) );
          }
          else
            return result_type();

        default: // coplanar
          // q belongs to the triangle's supporting plane
          // p sees the triangle in counterclockwise order
          if(orientation(p,q,a,b) != POSITIVE
             && orientation(p,q,b,c) != POSITIVE
             && orientation(p,q,c,a) != POSITIVE)
          {
            return result_type(q);
          }
          else return result_type();
      }
    case NEGATIVE:
      switch ( abcq ) {
        case POSITIVE:
          // q sees the triangle in counterclockwise order
          if ( orientation(q,p,a,b) != POSITIVE
            && orientation(q,p,b,c) != POSITIVE
            && orientation(q,p,c,a) != POSITIVE )
          {
            // The intersection is a point
            return result_type( lp_intersection(p, q, a, b, c, k) );
          }
          else
            return result_type();

        case NEGATIVE:
          // the segment lies in the negative open halfspaces defined by the
          // triangle's supporting plane
          return result_type();

        default: // coplanar
          // q belongs to the triangle's supporting plane
          // p sees the triangle in clockwise order
          if(orientation(q,p,a,b) != POSITIVE
             && orientation(q,p,b,c) != POSITIVE
             && orientation(q,p,c,a) != POSITIVE)
          {
            return result_type(q);
          }
          else return result_type();
      }
    default: // coplanar
      switch ( abcq ) {
      case POSITIVE:
        // q sees the triangle in counterclockwise order
        if(orientation(q,p,a,b) != POSITIVE
           && orientation(q,p,b,c) != POSITIVE
           && orientation(q,p,c,a) != POSITIVE)
        {
          return result_type(p);
        } else
          return result_type();
      case NEGATIVE:
        // q sees the triangle in clockwise order
        if(orientation(p,q,a,b) != POSITIVE
           && orientation(p,q,b,c) != POSITIVE
           && orientation(p,q,c,a) != POSITIVE)
        {
          return result_type(p);
        } else
          return result_type();
      case COPLANAR:
        // the segment is coplanar with the triangle's supporting plane
        // we test whether the segment intersects the triangle in the common
        // supporting plane
        //
        // Laurent Rineau, 2016/10/10: this case is purposely ignored by
        // Mesh_3, because if the intersection is not a point, it is
        // ignored anyway.
        return result_type();
      default: // should not happen.
        CGAL_kernel_assertion(false);
        return result_type();
      }
  }
}



  ////////////////


template <class K>
typename cpp11::result_of<
  typename K::Intersect_3(typename K::Ray_3, typename K::Triangle_3)>::type
tr_intersection(const typename K::Triangle_3  &t,
                const typename K::Ray_3 &r,
                const K& k)
{
  CGAL_MESH_3_BRANCH_PROFILER(std::string("coplanar/calls in   : ") +
                              std::string(CGAL_PRETTY_FUNCTION), tmp);
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(t) ) ;
  CGAL_kernel_precondition( ! k.is_degenerate_3_object()(r) ) ;

  typedef typename cpp11::result_of<
    typename K::Intersect_3(typename K::Ray_3, typename K::Triangle_3)
  >::type result_type;

  typedef typename K::Point_3 Point_3;

  typename K::Construct_vertex_3 vertex_on =
    k.construct_vertex_3_object();

  typename K::Orientation_3 orientation =
    k.orientation_3_object();

  typename K::Construct_point_on_3 point_on =
    k.construct_point_on_3_object();

  typename Mesh_3::Vector_plane_orientation_3_static_filter<K> vector_plane_orient;


  // Here we use the fact that K::Triangle_3 internally stores three
  // points, and that K::Ray_3 internally stores two points.
  const Point_3& a = vertex_on(t,0);
  const Point_3& b = vertex_on(t,1);
  const Point_3& c = vertex_on(t,2);
  const Point_3& p = point_on(r,0);
  const Point_3& q = point_on(r,1);

  const Orientation ray_direction = vector_plane_orient(p, q, a, b, c);

  if(ray_direction == COPLANAR) return result_type();

  const Orientation abcp = orientation(a,b,c,p);

  if(abcp == COPLANAR) return result_type(); // p belongs to the triangle's
                                        // supporting plane

  if(ray_direction == abcp) return result_type();
  // The ray lies entirely in one of the two open halfspaces defined by the
  // triangle's supporting plane.

  // Here we know that the ray crosses the plane (abc)

  if ( orientation(p,q,a,b) != abcp
       && orientation(p,q,b,c) != abcp
       && orientation(p,q,c,a) != abcp )
    return result_type(lp_intersection(p, q, a, b, c, k));
  else
    return result_type();
}

  ////////////////

template < typename K_ >
class Robust_intersection_3_new
{
public:
  typedef typename K_::FT                             FT;

  typedef typename K_::Triangle_3                     Triangle_3;
  typedef typename K_::Segment_3                      Segment_3;
  typedef typename K_::Ray_3                          Ray_3;

  template <typename>
  struct result;

  template <typename F, typename A, typename B>
  struct result<F(A, B)> {
    typedef typename cpp11::result_of<typename K_::Intersect_3(A, B)>::type type;
  };

  typedef Exact_predicates_exact_constructions_kernel     EK;
  typedef Cartesian_converter<typename K_::Kernel, EK>    To_exact;
  typedef Cartesian_converter<EK, typename K_::Kernel>    Back_from_exact;

  template<class T1, class T2>
  typename cpp11::result_of<typename K_::Intersect_3(T1, T2)>::type
  operator() (const T1& t, const T2& s) const
  {
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Intersect_3 exact_intersection = EK().intersect_3_object();

    // Cartesian converters have an undocumented, optional< variant > operator
    return typename cpp11::result_of<typename K_::Intersect_3(T1, T2)>::type
      (back_from_exact(exact_intersection(to_exact(t), to_exact(s))));
  }

  typename cpp11::result_of<typename K_::Intersect_3(Segment_3, Triangle_3)>::type
  operator()(const Segment_3& s, const Triangle_3& t) const
  {
    return ts_intersection(t, s, K_());
  }

  typename cpp11::result_of<typename K_::Intersect_3(Segment_3, Triangle_3)>::type
  operator()(const Triangle_3& t, const Segment_3& s) const
  {
    return ts_intersection(t, s, K_());
  }

  typename cpp11::result_of<typename K_::Intersect_3(Ray_3, Triangle_3)>::type
  operator()(const Ray_3& r, const Triangle_3& t) const  {
    return tr_intersection(t, r, K_());
  }

  typename cpp11::result_of<typename K_::Intersect_3(Ray_3, Triangle_3)>::type
  operator()(const Triangle_3& t, const Ray_3& r) const
  {
    return tr_intersection(t, r, K_());
  }

};

template < typename K_ >
class Robust_intersection_3
{
public:
  typedef typename K_::FT                             FT;

  typedef typename K_::Triangle_3                     Triangle_3;
  typedef typename K_::Segment_3                      Segment_3;

  template <typename>
  struct result;

  template <typename F, typename A, typename B>
  struct result<F(A, B)> {
    typedef typename cpp11::result_of<typename K_::Intersect_3(A, B)>::type type;
  };

  typedef Exact_predicates_exact_constructions_kernel   EK;
  typedef Cartesian_converter<typename K_::Kernel, EK>    To_exact;
  typedef Cartesian_converter<EK, typename K_::Kernel>    Back_from_exact;

  template<class T1, class T2>
  typename cpp11::result_of<typename K_::Intersect_3(T1, T2)>::type
  operator() (const T1& t, const T2& s) const
  {
    // Switch to exact
    To_exact to_exact;
    Back_from_exact back_from_exact;
    EK::Intersect_3 exact_intersection = EK().intersect_3_object();

    // Cartesian converters have an undocumented, optional< variant > operator
    return typename cpp11::result_of<typename K_::Intersect_3(T1, T2)>::type
      (back_from_exact(exact_intersection(to_exact(t), to_exact(s))));
  }
};



/**
 * @struct Robust_intersection_traits_3
 */
template<class K_>
struct Robust_intersection_traits_3
: public K_
{
  typedef Robust_intersection_3<K_> Intersect_3;
  typedef Robust_intersection_traits_3<K_> Kernel;
  Intersect_3
  intersect_3_object() const
  {
    return Intersect_3();
  }

};

template<class K_>
struct Robust_intersection_traits_3_new
: public K_
{
  typedef Robust_intersection_3_new<K_> Intersect_3;
  typedef Robust_intersection_traits_3_new<K_> Kernel;
  Intersect_3
  intersect_3_object() const
  {
    return Intersect_3();
  }

};

} // end namespace Mesh_3

} //namespace CGAL

#endif // CGAL_MESH_3_ROBUST_INTERSECTION_TRAITS_3_H
