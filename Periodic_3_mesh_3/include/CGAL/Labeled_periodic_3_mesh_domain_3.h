// Copyright (c) 2014 INRIA Sophia-Antipolis (France).
// Copyright (c) 2017 GeometryFactory (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Stéphane Tayeb,
//                 Aymeric Pellé
//
#ifndef CGAL_PERIODIC_3_MESH_3_LABELED_PERIODIC_3_MESH_DOMAIN_3_H
#define CGAL_PERIODIC_3_MESH_3_LABELED_PERIODIC_3_MESH_DOMAIN_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/Periodic_3_mesh_triangulation_3.h>
#include <CGAL/internal/canonicalize_helper.h>

#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/intersections.h>
#include <CGAL/result_of.h>

#include <algorithm>

namespace CGAL
{

/**
 * \class Labeled_periodic_3_mesh_domain_3
 *
 * Function `f` must be defined over the fundamental domain and take his values into `N`.
 * Let `p` be a point.
 *  - `f(p)=0` means that p is outside domain.
 *  - `f(p)=a`, `a!=0` means that `p` is inside subdomain `a`.
 *
 *  Any boundary facet is labelled `<a,b>`, with `a<b`, where `a` and `b` are the
 *  tags of its incident subdomain.
 *  Thus, a boundary facet of the domain is labelled `<0,b>`, where `b!=0`.
 */
// The difference with Mesh_3's is very tiny if 'CGAL_PERIODIC_CANONICALIZE_DUAL_INTERSECTIONS'
// is not enabled: it resides in the call to label() which canonicalizes the point
// (that is, send them to the fundamental domain) before evaluation.
//
// This could be simplified with a virtual function label() (in Mesh_3)?
// Leaving this simplification to be done when other types of domains are added to P3M3.
template<class Function, class BGT, class Null_subdomain_index = Default>
class Labeled_periodic_3_mesh_domain_3
  : public Labeled_mesh_domain_3<Function, BGT,
                                 typename Default::Get<Null_subdomain_index,
                                                       CGAL::Null_subdomain_index>::type>
{
public:
  /// Null subdomain type
  typedef typename Default::Get<Null_subdomain_index,
                                CGAL::Null_subdomain_index>::type       Null;

  /// Base type
  typedef Labeled_mesh_domain_3<Function, BGT, Null>                    Base;

  /// Geometric object types
  typedef typename Base::Point_3            Point_3;
  typedef typename Base::Segment_3          Segment_3;
  typedef typename Base::Ray_3              Ray_3;
  typedef typename Base::Line_3             Line_3;
  typedef typename Base::Vector_3           Vector_3;
  typedef typename Base::Iso_cuboid_3       Iso_cuboid_3;
  typedef typename Base::Bbox_3             Bbox_3;
  typedef typename Base::Sphere_3           Sphere_3;

  /// Type of indexes for surface patch of the input complex
  typedef typename Base::Surface_patch       Surface_patch;
  typedef typename Base::Surface_patch_index Surface_patch_index;

  /// Type of indexes to characterize the lowest dimensional face of the input
  /// complex on which a vertex lie
  typedef typename Base::Index               Index;
  typedef typename Base::Intersection        Intersection;

  /// Type of indexes for cells of the input complex
  typedef typename Base::Subdomain_index     Subdomain_index;
  typedef typename Base::Subdomain           Subdomain;

  /// Public types
  typedef typename Base::FT                  FT;
  typedef BGT                                Geom_traits;

  /// Periodic traits used in the canonicalization of the points
  typedef typename details::Periodic_3_mesh_geom_traits_generator<BGT>::type Periodic_geom_traits;

  Labeled_periodic_3_mesh_domain_3(const Function& f,
                                   const Iso_cuboid_3& bbox,
                                   const FT& error_bound = FT(1e-6),
                                   Null null = Null(),
                                   CGAL::Random* p_rng = NULL)
    : Base(f, bbox, error_bound, null, p_rng),
      pgt(bbox)
  { }

  const Iso_cuboid_3& canonical_periodic_domain() const { return Base::bounding_box(); }
  const Periodic_geom_traits& periodic_geom_traits() const { return pgt; }

  Subdomain_index label(const Point_3& p) const {
    return Base::labeling_function()(P3T3::internal::robust_canonicalize_point(p, pgt));
  }
  Subdomain_index label(const Point_3& p, const bool b) const {
    return Base::labeling_function()(P3T3::internal::robust_canonicalize_point(p, pgt), b);
  }

  /**
   * Constructs  a set of \ccc{n} points on the surface, and output them to
   *  the output iterator \ccc{pts} whose value type is required to be
   *  \ccc{std::pair<Points_3, Index>}.
   */
  // This is a complete copy paste from the base class because Do_intersect
  // needs to be the nested class' and not the base class'.
  // To be simplified...
  struct Construct_initial_points
  {
    Construct_initial_points(const Labeled_periodic_3_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int nb_points = 12) const
    {
      // Create point_iterator on and in bounding_sphere
      typedef Random_points_on_sphere_3<Point_3> Random_points_on_sphere_3;
      typedef Random_points_in_sphere_3<Point_3> Random_points_in_sphere_3;

      const FT squared_radius = BGT().compute_squared_radius_3_object()(
          r_domain_.bounding_sphere(r_domain_.bounding_box()));

      const double radius = std::sqrt(CGAL::to_double(squared_radius));

      CGAL::Random& rng = *(r_domain_.random_number_generator());
      Random_points_on_sphere_3 random_point_on_sphere(radius, rng);
      Random_points_in_sphere_3 random_point_in_sphere(radius, rng);

      // Get some functors
      typename Periodic_geom_traits::Construct_segment_3 segment_3 =
        r_domain_.periodic_geom_traits().construct_segment_3_object();
      typename Periodic_geom_traits::Construct_vector_3 vector_3 =
        r_domain_.periodic_geom_traits().construct_vector_3_object();
      typename Periodic_geom_traits::Construct_translated_point_3 translate =
        r_domain_.periodic_geom_traits().construct_translated_point_3_object();
      typename Periodic_geom_traits::Construct_center_3 center =
        r_domain_.periodic_geom_traits().construct_center_3_object();

      // Get translation from origin to sphere center
      Point_3 center_pt = center(r_domain_.bounding_sphere(r_domain_.bounding_box()));
      const Vector_3 sphere_translation = vector_3(CGAL::ORIGIN, center_pt);

      // Create nb_point points
      int n = nb_points;
#ifdef CGAL_MESH_3_VERBOSE
      std::cerr << "construct initial points (nb_points: " << nb_points << ")\n";
#endif
      while ( 0 != n )
      {
        // Get a random segment
        const Point_3 random_point = translate(*random_point_on_sphere, sphere_translation);
        const Segment_3 random_seg = segment_3(center_pt, random_point);

        // Add the intersection to the output if it exists
        Surface_patch surface = r_domain_.do_intersect_surface_object()(random_seg);
        if ( surface )
        {
          const Point_3 intersect_pt = CGAL::cpp11::get<0>(
              r_domain_.construct_intersection_object()(random_seg));
          *pts++ = std::make_pair(intersect_pt,
                                  r_domain_.index_from_surface_patch_index(*surface));
          --n;

#ifdef CGAL_MESH_3_VERBOSE
          std::cerr << boost::format("\r             \r"
                                     "%1%/%2% initial point(s) found...")
                       % (nb_points - n)
                       % nb_points;
#endif
        }
        else
        {
          // Get a new random point into sphere as center of object
          // It may be necessary if the center of the domain is empty, e.g. torus
          // In general case, it is good for input point dispersion
          ++random_point_in_sphere;
          center_pt = translate(*random_point_in_sphere, sphere_translation);
        }
        ++random_point_on_sphere;
      }

#ifdef CGAL_MESH_3_VERBOSE
      std::cerr << "\n";
#endif
      return pts;
    }

  private:
    const Labeled_periodic_3_mesh_domain_3& r_domain_;
  };

  /// Returns Construct_initial_points object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  /**
   * Returns true if point~\ccc{p} is in the domain. If \ccc{p} is in the
   *  domain, the parameter index is set to the index of the subdomain
   *  including $p$. It is set to the default value otherwise.
   */
  struct Is_in_domain
  {
    Is_in_domain(const Labeled_periodic_3_mesh_domain_3& domain)
      : r_domain_(domain) {}

    Subdomain operator()(const Point_3& p) const
    {
      // null(f(p)) means p is outside the domain
      Subdomain_index index = r_domain_.label(p);
      if ( r_domain_.null_function()(index) )
        return Subdomain();
      else
        return Subdomain(index);
    }
  private:
    const Labeled_periodic_3_mesh_domain_3& r_domain_;
  };

  /// Returns Is_in_domain object
  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }

  /**
   * Returns true if the element \ccc{type} intersect properly any of the
   * surface patches describing the either the boundary of the domain
   * or some subdomain boundary.
   * \ccc{Type} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * Parameter index is set to the index of the intersected surface patch
   * if \ccc{true} is returned and to the default \ccc{Surface_patch_index}
   * value otherwise.
   */
  struct Do_intersect_surface
  {
    Do_intersect_surface(const Labeled_periodic_3_mesh_domain_3& domain)
      : r_domain_(domain)
    { }

    Surface_patch operator()(const Segment_3& s) const
    {
      return this->operator()(s.source(), s.target());
    }

    Surface_patch operator()(const Ray_3& r) const
    {
      return clip_to_segment(r);
    }

    Surface_patch operator()(const Line_3& l) const
    {
      return clip_to_segment(l);
    }

  private:
    /// Returns true if the points \c a & \c b do not belong to the same subdomain
    /// \c index is set to the surface index of subdomains f(a), f(b)
    Surface_patch operator()(const Point_3& a, const Point_3& b) const
    {
      // If f(a) != f(b), then [a,b] intersects some surface. Here we consider
      // [a,b] intersects surface_patch labelled <f(a),f(b)> (or <f(b),f(a)>).
      // It may be false, further rafinement will improve precision

      // This whole canonicalization of offsets process seems useless... Hiding it behind macros.
#ifdef CGAL_PERIODIC_CANONICALIZE_DUAL_INTERSECTIONS
      Iso_cuboid_3 pbb = r_domain_.canonical_periodic_domain();
      FT dimension [3] = { pbb.xmax()-pbb.xmin(),
                           pbb.ymax()-pbb.ymin(),
                           pbb.zmax()-pbb.zmin() };
      FT a_t [3] = { a.x(), a.y(), a.z() };
      FT b_t [3] = { b.x(), b.y(), b.z() };
      a_t[0] -= pbb.xmin();
      a_t[1] -= pbb.ymin();
      a_t[2] -= pbb.zmin();
      b_t[0] -= pbb.xmin();
      b_t[1] -= pbb.ymin();
      b_t[2] -= pbb.zmin();
      int o1 [3] = { static_cast<int>(a_t[0] / dimension[0]),
                     static_cast<int>(a_t[1] / dimension[1]),
                     static_cast<int>(a_t[2] / dimension[2]) };
      int o2 [3] = { static_cast<int>(b_t[0] / dimension[0]),
                     static_cast<int>(b_t[1] / dimension[1]),
                     static_cast<int>(b_t[2] / dimension[2]) };
      FT a_min [3] = { a.x(), a.y(), a.z() };
      FT b_min [3] = { b.x(), b.y(), b.z() };

      for (unsigned idx = 0; idx < 3; ++idx)
      {
        FT offset = dimension[idx] * static_cast<FT>((std::min)(o1[idx], o2[idx]));
        a_min[idx] -= offset;
        b_min[idx] -= offset;
      }

      const Point_3 pa(a_min[0], a_min[1], a_min[2]);
      const Point_3 pb(b_min[0], b_min[1], b_min[2]);
      Subdomain_index value_a = r_domain_.label(pa);
      Subdomain_index value_b = r_domain_.label(pb);
#else
      Subdomain_index value_a = r_domain_.label(a);
      Subdomain_index value_b = r_domain_.label(b);
#endif

      if ( value_a != value_b )
      {
        if( r_domain_.null_function()(value_a) && r_domain_.null_function()(value_b) )
          return Surface_patch();
        else
          return Surface_patch(r_domain_.make_surface_index(value_a, value_b));
      }

#ifdef CGAL_PERIODIC_CANONICALIZE_DUAL_INTERSECTIONS
      FT a_max [3] = { a.x(), a.y(), a.z() };
      FT b_max [3] = { b.x(), b.y(), b.z() };

      for (unsigned idx = 0; idx < 3; ++idx)
      {
        FT offset = dimension[idx] * static_cast<FT>((std::max)(o1[idx], o2[idx]));
        a_max[idx] -= offset;
        b_max[idx] -= offset;
      }

      pa = Point_3(a_max[0], a_max[1], a_max[2]);
      pb = Point_3(b_max[0], b_max[1], b_max[2]);
      value_a = r_domain_.label(pa);
      value_b = r_domain_.label(pb);
      if ( value_a != value_b )
      {
        if( r_domain_.null_function()(value_a) && r_domain_.null_function()(value_b) )
          return Surface_patch();
        else
          return Surface_patch(r_domain_.make_surface_index(value_a, value_b));
      }
#endif

      return Surface_patch();
    }

    /**
     * Clips \c query to a segment \c s, and call operator()(s)
     */
    template<typename Query>
    Surface_patch clip_to_segment(const Query& query) const
    {
      typename cpp11::result_of<typename BGT::Intersect_3(Query, Iso_cuboid_3)>::type
        clipped = CGAL::intersection(query, r_domain_.bounding_box());

      if(clipped)
#if CGAL_INTERSECTION_VERSION > 1
        if(const Segment_3* s = boost::get<Segment_3>(&*clipped))
          return this->operator()(*s);
#else
        if(const Segment_3* s = object_cast<Segment_3>(&clipped))
          return this->operator()(*s);
#endif

      return Surface_patch();
    }

  private:
    const Labeled_periodic_3_mesh_domain_3& r_domain_;
  };

  /// Returns Do_intersect_surface object
  Do_intersect_surface do_intersect_surface_object() const
  {
    return Do_intersect_surface(*this);
  }

  /**
   * Returns a point in the intersection of the primitive \ccc{type}
   * with some boundary surface.
   * \ccc{Type1} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * The integer \ccc{dimension} is set to the dimension of the lowest
   * dimensional face in the input complex containing the returned point, and
   * \ccc{index} is set to the index to be stored at a mesh vertex lying
   * on this face.
   */
  struct Construct_intersection
  {
    Construct_intersection(const Labeled_periodic_3_mesh_domain_3& domain)
      : r_domain_(domain)
    { }

    Intersection operator()(const Segment_3& s) const
    {
#ifndef CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      CGAL_precondition(r_domain_.do_intersect_surface_object()(s));
#endif // NOT CGAL_MESH_3_NO_LONGER_CALLS_DO_INTERSECT_3
      return this->operator()(s.source(),s.target());
    }

    Intersection operator()(const Ray_3& r) const
    {
      return clip_to_segment(r);
    }

    Intersection operator()(const Line_3& l) const
    {
      return clip_to_segment(l);
    }

  private:
    /**
     * Returns a point in the intersection of [a,b] with the surface
     * \c a must be the source point, and \c b the target point. It's important
     * because it drives bisection cuts.
     * Indeed, the returned point is the first intersection from \c a on \c [a,b]
     * with a subdomain surface.
     */
    Intersection operator()(const Point_3& a, const Point_3& b) const
    {
      // Functors
      typename Periodic_geom_traits::Compute_squared_distance_3 squared_distance =
        r_domain_.periodic_geom_traits().compute_squared_distance_3_object();
      typename Periodic_geom_traits::Construct_midpoint_3 midpoint =
        r_domain_.periodic_geom_traits().construct_midpoint_3_object();

#ifdef CGAL_PERIODIC_CANONICALIZE_DUAL_INTERSECTIONS
      Iso_cuboid_3 pbb = r_domain_.canonical_periodic_domain();
      FT dimension [3] = { pbb.xmax() - pbb.xmin(),
                           pbb.ymax() - pbb.ymin(),
                           pbb.zmax() - pbb.zmin() };
      FT a_t [3] = { a.x(), a.y(), a.z() };
      FT b_t [3] = { b.x(), b.y(), b.z() };
      a_t[0] -= pbb.xmin();
      a_t[1] -= pbb.ymin();
      a_t[2] -= pbb.zmin();
      b_t[0] -= pbb.xmin();
      b_t[1] -= pbb.ymin();
      b_t[2] -= pbb.zmin();

      int o1 [3] = { static_cast<int>(a_t[0] / dimension[0]),
                     static_cast<int>(a_t[1] / dimension[1]),
                     static_cast<int>(a_t[2] / dimension[2]) };
      int o2 [3] = { static_cast<int>(b_t[0] / dimension[0]),
                     static_cast<int>(b_t[1] / dimension[1]),
                     static_cast<int>(b_t[2] / dimension[2]) };

      FT a_min [3] = { a.x(), a.y(), a.z() };
      FT b_min [3] = { b.x(), b.y(), b.z() };
      for (unsigned idx = 0; idx < 3; ++idx)
      {
        FT offset = dimension[idx] * static_cast<FT>((std::min)(o1[idx], o2[idx]));
        a_min[idx] -= offset;
        b_min[idx] -= offset;
      }

      // Non const points
      Point_3 p1(a_min[0], a_min[1], a_min[2]);
      Point_3 p2(b_min[0], b_min[1], b_min[2]);
#else
      Point_3 p1 = a;
      Point_3 p2 = b;
#endif

      Point_3 mid = midpoint(p1, p2);

      // Cannot be const: those values are modified below.
      Subdomain_index value_at_p1 = r_domain_.label(p1);
      Subdomain_index value_at_p2 = r_domain_.label(p2);
      Subdomain_index value_at_mid = r_domain_.label(mid, true);

      // If both extremities are in the same subdomain, then there is no intersection.
      // This should not happen...
      if( value_at_p1 == value_at_p2 )
      {
#ifdef CGAL_PERIODIC_CANONICALIZE_DUAL_INTERSECTIONS
        FT a_max [3] = { a.x(), a.y(), a.z() };
        FT b_max [3] = { b.x(), b.y(), b.z() };
        for (unsigned idx = 0; idx < 3; ++idx)
        {
          FT offset = dimension[idx] * static_cast<FT>((std::max)(o1[idx], o2[idx]));
          a_max[idx] -= offset;
          b_max[idx] -= offset;
        }

        p1 = Point_3(a_max[0], a_max[1], a_max[2]);
        p2 = Point_3(b_max[0], b_max[1], b_max[2]);
        mid = midpoint(p1, p2);

        value_at_p1 = r_domain_.label(p1);
        value_at_p2 = r_domain_.label(p2);
        value_at_mid = r_domain_.label(mid, true);

        if( value_at_p1 == value_at_p2 )
#endif
          return Intersection();
      }

      if( r_domain_.null_function()(value_at_p1) && r_domain_.null_function()(value_at_p2) ) {
        return Intersection();
      }

      // Else lets find a point (by bisection)
      // Bisection ends when the point is closer to the surface than the error bound
      while(true)
      {
        // If the two points are sufficiently close, then we return the midpoint
        if ( squared_distance(p1, p2) < r_domain_.squared_error_bound_value() )
        {
          CGAL_assertion(value_at_p1 != value_at_p2 &&
                         ! ( r_domain_.null_function()(value_at_p1) && r_domain_.null_function()(value_at_p2) ) );

          const Surface_patch_index sp_index =
            r_domain_.make_surface_index(value_at_p1, value_at_p2);
          const Index index = r_domain_.index_from_surface_patch_index(sp_index);

          return Intersection(mid, index, 2);
        }

        // Otherwise, we must go on...
        // Here we consider that p1(a) is the source point. Thus, we keep p1 and
        // change p2 if f(p1)!=f(p2).
        // That allows us to find the first intersection from a of [a,b] with
        // a surface.
        if ( value_at_p1 != value_at_mid &&
             ! ( r_domain_.null_function()(value_at_p1) && r_domain_.null_function()(value_at_mid) ) )
        {
          p2 = mid;
          value_at_p2 = value_at_mid;
        }
        else
        {
          p1 = mid;
          value_at_p1 = value_at_mid;
        }

        mid = midpoint(p1, p2);
        value_at_mid = r_domain_.label(mid, true);
      }
    }

    /// Clips \c query to a segment \c s, and call operator()(s)
    template<typename Query>
    Intersection clip_to_segment(const Query& query) const
    {
      typename cpp11::result_of<typename BGT::Intersect_3(Query, Iso_cuboid_3)>::type
        clipped = CGAL::intersection(query, r_domain_.bounding_box());

      if(clipped)
#if CGAL_INTERSECTION_VERSION > 1
        if(const Segment_3* s = boost::get<Segment_3>(&*clipped))
          return this->operator()(*s);
#else
        if(const Segment_3* s = object_cast<Segment_3>(&clipped))
          return this->operator()(*s);
#endif

      return Intersection();
    }

  private:
    const Labeled_periodic_3_mesh_domain_3& r_domain_;
  };

  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

private:
  // The domain must know the periodic fundamental domain, which is in the traits.
  //
  // Note that the fundamental domain is also known through 'bounding_box()', but
  // a full geometric traits class is needed for 'construct_point()'
  // in canonicalization functions.
  Periodic_geom_traits pgt;

private:
  // Disabled copy constructor & assignment operator
  Labeled_periodic_3_mesh_domain_3(const Labeled_periodic_3_mesh_domain_3&);
  Labeled_periodic_3_mesh_domain_3& operator=(const Labeled_periodic_3_mesh_domain_3&);
};

} // namespace CGAL

#endif // CGAL_PERIODIC_3_MESH_3_LABELED_PERIODIC_3_MESH_DOMAIN_3_H
