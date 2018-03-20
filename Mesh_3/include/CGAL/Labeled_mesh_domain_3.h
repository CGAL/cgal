// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
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
//
// Author(s)     : Stéphane Tayeb, Aymeric PELLE
//
//******************************************************************************
// File Description :
// class Labeled_mesh_domain_3. See class description.
//******************************************************************************

#ifndef CGAL_LABELED_MESH_DOMAIN_3_H
#define CGAL_LABELED_MESH_DOMAIN_3_H

#include <CGAL/license/Mesh_3.h>

#include <CGAL/disable_warnings.h>

#include <CGAL/Mesh_3/config.h>

#include <CGAL/Bbox_3.h>

#include <CGAL/point_generators_3.h>

#include <boost/variant.hpp>
#include <boost/format.hpp>
#include <boost/optional.hpp>
#include <CGAL/tuple.h>
#include <CGAL/Origin.h>

#include <CGAL/Default.h>

#include <CGAL/internal/Mesh_3/Handle_IO_for_pair_of_int.h>

namespace CGAL {

struct Null_subdomain_index {
  template <typename T>
  bool operator()(const T& x) const { return 0 == x; }
};

/**
 * \class Labeled_mesh_domain_3
 *
 * Function f must take his values into N.
 * Let p be a Point.
 *  - f(p)=0 means that p is outside domain.
 *  - f(p)=a, a!=0 means that p is inside subdomain a.
 *
 *  Any boundary facet is labelled <a,b>, a<b, where a and b are the
 *  tags of it's incident subdomain.
 *  Thus, a boundary facet of the domain is labelled <0,b>, where b!=0.
 */
template<class Function, class BGT,
         class Null_subdomain_index = Default >
class Labeled_mesh_domain_3
{
  typedef typename Default::Get<Null_subdomain_index,
                                CGAL::Null_subdomain_index>::type Null;
public:
  /// Geometric object types
  typedef typename BGT::Point_3    Point_3;
  typedef typename BGT::Segment_3  Segment_3;
  typedef typename BGT::Ray_3      Ray_3;
  typedef typename BGT::Line_3     Line_3;
  typedef typename BGT::Vector_3   Vector_3;
  typedef typename BGT::Sphere_3   Sphere_3;
  typedef CGAL::Bbox_3             Bbox_3;

protected:
  typedef typename BGT::Iso_cuboid_3 Iso_cuboid_3;

public:
  // Kernel_traits compatibility
  typedef BGT R;
  // access Function type from inherited class
  typedef Function Fct;

  //-------------------------------------------------------
  // Index Types
  //-------------------------------------------------------
  /// Type of indexes for cells of the input complex
  typedef typename Function::return_type Subdomain_index;
  typedef boost::optional<Subdomain_index> Subdomain;
  /// Type of indexes for surface patch of the input complex
  typedef std::pair<Subdomain_index, Subdomain_index> Surface_patch_index;
  typedef boost::optional<Surface_patch_index> Surface_patch;
  /// Type of indexes to characterize the lowest dimensional face of the input
  /// complex on which a vertex lie
  typedef boost::variant<Subdomain_index, Surface_patch_index> Index;
  typedef CGAL::cpp11::tuple<Point_3,Index,int> Intersection;


  typedef typename BGT::FT FT;
  typedef BGT Geom_traits;

  /**
   * @brief Constructor
   */
  Labeled_mesh_domain_3(const Function& f,
                        const Sphere_3& bounding_sphere,
                        const FT& error_bound = FT(1e-3),
                        Null null = Null(),
                        CGAL::Random* p_rng = NULL);

  Labeled_mesh_domain_3(const Function& f,
                        const Bbox_3& bbox,
                        const FT& error_bound = FT(1e-3),
                        Null null = Null(),
                        CGAL::Random* p_rng = NULL);

  Labeled_mesh_domain_3(const Function& f,
                        const Iso_cuboid_3& bbox,
                        const FT& error_bound = FT(1e-3),
                        Null null = Null(),
                        CGAL::Random* p_rng = NULL);

  Labeled_mesh_domain_3(const Function& f,
                        const Sphere_3& bounding_sphere,
                        const FT& error_bound,
                        CGAL::Random* p_rng);

  Labeled_mesh_domain_3(const Function& f,
                        const Bbox_3& bbox,
                        const FT& error_bound,
                        CGAL::Random* p_rng);

  Labeled_mesh_domain_3(const Function& f,
                        const Iso_cuboid_3& bbox,
                        const FT& error_bound,
                        CGAL::Random* p_rng);
  /// Destructor
  virtual ~Labeled_mesh_domain_3()
  {
    if(delete_rng_)
      delete p_rng_;
  }


  /**
   * Constructs  a set of \ccc{n} points on the surface, and output them to
   *  the output iterator \ccc{pts} whose value type is required to be
   *  \ccc{std::pair<Points_3, Index>}.
   */
  struct Construct_initial_points
  {
    Construct_initial_points(const Labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

    template<class OutputIterator>
    OutputIterator operator()(OutputIterator pts, const int n = 12) const;

  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  /// Returns Construct_initial_points object
  Construct_initial_points construct_initial_points_object() const
  {
    return Construct_initial_points(*this);
  }

  /**
   * Returns a bounding box of the domain
   */
  Bbox_3 bbox() const {
    return bbox_.bbox();
  }

  /**
   * Returns true if point~\ccc{p} is in the domain. If \ccc{p} is in the
   *  domain, the parameter index is set to the index of the subdomain
   *  including $p$. It is set to the default value otherwise.
   */
  struct Is_in_domain
  {
    Is_in_domain(const Labeled_mesh_domain_3& domain) : r_domain_(domain) {}

    Subdomain operator()(const Point_3& p) const
    {
      // null(f(p)) means p is outside the domain
      Subdomain_index index = (r_domain_.function_)(p);
      if ( r_domain_.null(index) )
        return Subdomain();
      else
        return Subdomain(index);
    }
  private:
    const Labeled_mesh_domain_3& r_domain_;
  };

  /// Returns Is_in_domain object
  Is_in_domain is_in_domain_object() const { return Is_in_domain(*this); }

  /**
   * Returns true is the element \ccc{type} intersect properly any of the
   * surface patches describing the either the domain boundary or some
   * subdomain boundary.
   * \ccc{Type} is either \ccc{Segment_3}, \ccc{Ray_3} or \ccc{Line_3}.
   * Parameter index is set to the index of the intersected surface patch
   * if \ccc{true} is returned and to the default \ccc{Surface_patch_index}
   * value otherwise.
   */
  struct Do_intersect_surface
  {
    Do_intersect_surface(const Labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

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
    /// Returns true if points \c a & \c b do not belong to the same subdomain
    /// \c index is set to the surface index of subdomains f(a), f(b)
    Surface_patch operator()(const Point_3& a, const Point_3& b) const
    {
      // If f(a) != f(b), then [a,b] intersects some surface. Here we consider
      // [a,b] intersects surface_patch labelled <f(a),f(b)> (or <f(b),f(a)>).
      // It may be false, further rafinement will improve precision
      const Subdomain_index value_a = r_domain_.function_(a);
      const Subdomain_index value_b = r_domain_.function_(b);

      if ( value_a != value_b ) {
        if( r_domain_.null(value_a) && r_domain_.null(value_b) )
          return Surface_patch();
        else
          return Surface_patch(r_domain_.make_surface_index(value_a, value_b));
      }
      else
        return Surface_patch();
    }

    /**
     * Clips \c query to a segment \c s, and call operator()(s)
     */
    template<typename Query>
    Surface_patch clip_to_segment(const Query& query) const
    {
      typename cpp11::result_of<typename BGT::Intersect_3(Query, Iso_cuboid_3)>::type
        clipped = CGAL::intersection(query, r_domain_.bbox_);

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
    const Labeled_mesh_domain_3& r_domain_;
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
    Construct_intersection(const Labeled_mesh_domain_3& domain)
      : r_domain_(domain) {}

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
     * \c a must be the source point, and \c b the out point. It's important
     * because it drives bisection cuts.
     * Indeed, the returned point is the first intersection from \c [a,b]
     * with a subdomain surface.
     */
    Intersection operator()(const Point_3& a, const Point_3& b) const
    {
      // Functors
      typename BGT::Compute_squared_distance_3 squared_distance =
                                      BGT().compute_squared_distance_3_object();
      typename BGT::Construct_midpoint_3 midpoint =
                                      BGT().construct_midpoint_3_object();

      // Non const points
      Point_3 p1 = a;
      Point_3 p2 = b;
      Point_3 mid = midpoint(p1, p2);

      // Cannot be const: those values are modified below.
      Subdomain_index value_at_p1 = r_domain_.function_(p1);
      Subdomain_index value_at_p2 = r_domain_.function_(p2);
      Subdomain_index value_at_mid = r_domain_.function_(mid,true);

      // If both extremities are in the same subdomain,
      // there is no intersection.
      // This should not happen...
      if( value_at_p1 == value_at_p2 )
      {
        return Intersection();
      }
      if( r_domain_.null(value_at_p1) && r_domain_.null(value_at_p2) ) {
        return Intersection();
      }

      // Else lets find a point (by bisection)
      // Bisection ends when the point is near than error bound from surface
      while(true)
      {
        // If the two points are enough close, then we return midpoint
        if ( squared_distance(p1, p2) < r_domain_.squared_error_bound_ )
        {
          CGAL_assertion(value_at_p1 != value_at_p2 &&
             ! ( r_domain_.null(value_at_p1) && r_domain_.null(value_at_p2) ) );
          const Surface_patch_index sp_index =
            r_domain_.make_surface_index(value_at_p1, value_at_p2);
          const Index index = r_domain_.index_from_surface_patch_index(sp_index);
          return Intersection(mid, index, 2);
        }

        // Else we must go on
        // Here we consider that p1(a) is the source point. Thus, we keep p1 and
        // change p2 if f(p1)!=f(p2).
        // That allows us to find the first intersection from a of [a,b] with
        // a surface.
        if ( value_at_p1 != value_at_mid &&
             ! ( r_domain_.null(value_at_p1) && r_domain_.null(value_at_mid) ) )
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
        value_at_mid = r_domain_.function_(mid,true);
      }
    }

    /// Clips \c query to a segment \c s, and call operator()(s)
    template<typename Query>
    Intersection clip_to_segment(const Query& query) const
    {
      typename cpp11::result_of<typename BGT::Intersect_3(Query, Iso_cuboid_3)>::type
        clipped = CGAL::intersection(query, r_domain_.bbox_);

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
    const Labeled_mesh_domain_3& r_domain_;
  };

  /// Returns Construct_intersection object
  Construct_intersection construct_intersection_object() const
  {
    return Construct_intersection(*this);
  }

  /**
   * Returns the index to be stored in a vertex lying on the surface identified
   * by \c index.
   */
  Index index_from_surface_patch_index(const Surface_patch_index& index) const
  { return Index(index); }

  /**
   * Returns the index to be stored in a vertex lying in the subdomain
   * identified by \c index.
   */
  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return Index(index); }

  /**
   * Returns the \c Surface_patch_index of the surface patch
   * where lies a vertex with dimension 2 and index \c index.
   */
  Surface_patch_index surface_patch_index(const Index& index) const
  { return boost::get<Surface_patch_index>(index); }

  /**
   * Returns the index of the subdomain containing a vertex
   *  with dimension 3 and index \c index.
   */
  Subdomain_index subdomain_index(const Index& index) const
  { return boost::get<Subdomain_index>(index); }
  
  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;
  
  Index index_from_surface_index(const Surface_index& index) const
  { return index_from_surface_patch_index(index); }
  
  Surface_index surface_index(const Index& index) const
  { return surface_patch_index(index); }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------


private:
  /// Returns Surface_patch_index from \c i and \c j
  Surface_patch_index make_surface_index(const Subdomain_index i,
                                   const Subdomain_index j) const
  {
    if ( i < j ) return Surface_patch_index(i,j);
    else return Surface_patch_index(j,i);
  }

  /// Returns squared error bound from \c bbox and \c error
  FT squared_error_bound(const Iso_cuboid_3& bbox, const FT& error) const
  {
    typename BGT::Compute_squared_distance_3 squared_distance =
                                    BGT().compute_squared_distance_3_object();
    return squared_distance((bbox.min)(), (bbox.max)())*error*error/4;
  }

  /// Returns squared error bound from \c sphere and \c error
  FT squared_error_bound(const Sphere_3& sphere, const FT& error) const
  {
    typename BGT::Compute_squared_radius_3 squared_radius =
                                    BGT().compute_squared_radius_3_object();
    return squared_radius(sphere)*error*error;
  }

  /// Returns the bounding sphere of an Iso_cuboid_3
  Sphere_3 bounding_sphere(const Iso_cuboid_3& bbox) const
  {
    typename BGT::Construct_sphere_3 sphere = BGT().construct_sphere_3_object();
    return sphere((bbox.min)(), (bbox.max)());
  }

  /// Returns and Iso_cuboid_3 from a Bbox_3
  Iso_cuboid_3 iso_cuboid(const Bbox_3& bbox)
  {
    const Point_3 p_min(bbox.xmin(), bbox.ymin(), bbox.zmin());
    const Point_3 p_max(bbox.xmax(), bbox.ymax(), bbox.zmax());

    return Iso_cuboid_3(p_min,p_max);
  }

protected:
  /// Returns bounding box
  const Iso_cuboid_3& bounding_box() const { return bbox_; }

private:
  /// The function which answers subdomain queries
  const Function function_;
  /// The bounding box
  const Iso_cuboid_3 bbox_;
  /// The functor that decides which sub-domain indices correspond to the
  /// outside of the domain.
  Null null;
  /// The random number generator used by Construct_initial_points
  CGAL::Random* p_rng_;
  bool delete_rng_;
  /// Error bound relative to sphere radius
  FT squared_error_bound_;


private:
  // Disabled copy constructor & assignment operator
  typedef Labeled_mesh_domain_3<Function,BGT> Self;
  Labeled_mesh_domain_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Labeled_mesh_domain_3




//-------------------------------------------------------
// Method implementation
//-------------------------------------------------------
template<class F, class BGT, class Null>
Labeled_mesh_domain_3<F,BGT,Null>::Labeled_mesh_domain_3(
                       const F& f,
                       const Sphere_3& bounding_sphere,
                       const FT& error_bound,
                       Null null /* = Null() */,
                       CGAL::Random* p_rng)
: function_(f)
, bbox_(iso_cuboid(bounding_sphere.bbox()))
, null(null)
, p_rng_(p_rng == 0 ? new CGAL::Random(0) : p_rng)
, delete_rng_(p_rng == 0)
, squared_error_bound_(squared_error_bound(bounding_sphere,error_bound))
{
}

template<class F, class BGT, class Null>
Labeled_mesh_domain_3<F,BGT,Null>::Labeled_mesh_domain_3(
                       const F& f,
                       const Bbox_3& bbox,
                       const FT& error_bound,
                       Null null /* = Null() */,
                       CGAL::Random* p_rng)
: function_(f)
, bbox_(iso_cuboid(bbox))
, null(null)
, p_rng_(p_rng == 0 ? new CGAL::Random(0) : p_rng)
, delete_rng_(p_rng == 0)
, squared_error_bound_(squared_error_bound(bbox_,error_bound))
{
}

template<class F, class BGT, class Null>
Labeled_mesh_domain_3<F,BGT,Null>::Labeled_mesh_domain_3(
                       const F& f,
                       const Iso_cuboid_3& bbox,
                       const FT& error_bound,
                       Null null /* = Null() */,
                       CGAL::Random* p_rng)
: function_(f)
, bbox_(bbox)
, null(null)
, p_rng_(p_rng == 0 ? new CGAL::Random(0) : p_rng)
, delete_rng_(p_rng == 0)
, squared_error_bound_(squared_error_bound(bbox_,error_bound))
{
}




template<class F, class BGT, class Null>
Labeled_mesh_domain_3<F,BGT,Null>::Labeled_mesh_domain_3(
                       const F& f,
                       const Sphere_3& bounding_sphere,
                       const FT& error_bound,
                       CGAL::Random* p_rng)
: function_(f)
, bbox_(iso_cuboid(bounding_sphere.bbox()))
, null(Null())
, p_rng_(p_rng == 0 ? new CGAL::Random(0) : p_rng)
, delete_rng_(p_rng == 0)
, squared_error_bound_(squared_error_bound(bounding_sphere,error_bound))
{
}

template<class F, class BGT, class Null>
Labeled_mesh_domain_3<F,BGT,Null>::Labeled_mesh_domain_3(
                       const F& f,
                       const Bbox_3& bbox,
                       const FT& error_bound,
                       CGAL::Random* p_rng)
: function_(f)
, bbox_(iso_cuboid(bbox))
, null(Null())
, p_rng_(p_rng == 0 ? new CGAL::Random(0) : p_rng)
, delete_rng_(p_rng == 0)
, squared_error_bound_(squared_error_bound(bbox_,error_bound))
{
}

template<class F, class BGT, class Null>
Labeled_mesh_domain_3<F,BGT,Null>::Labeled_mesh_domain_3(
                       const F& f,
                       const Iso_cuboid_3& bbox,
                       const FT& error_bound,
                       CGAL::Random* p_rng)
: function_(f)
, bbox_(bbox)
, null(Null())
, p_rng_(p_rng == 0 ? new CGAL::Random(0) : p_rng)
, delete_rng_(p_rng == 0)
, squared_error_bound_(squared_error_bound(bbox_,error_bound))
{
}


template<class F, class BGT, class Null>
template<class OutputIterator>
OutputIterator
Labeled_mesh_domain_3<F,BGT,Null>::Construct_initial_points::operator()(
                                                    OutputIterator pts,
                                                    const int nb_points) const
{
  // Create point_iterator on and in bounding_sphere
  typedef Random_points_on_sphere_3<Point_3> Random_points_on_sphere_3;
  typedef Random_points_in_sphere_3<Point_3> Random_points_in_sphere_3;


  const FT squared_radius = BGT().compute_squared_radius_3_object()(
      r_domain_.bounding_sphere(r_domain_.bbox_));

  const double radius = std::sqrt(CGAL::to_double(squared_radius));

  CGAL::Random& rng = *(r_domain_.p_rng_);
  Random_points_on_sphere_3 random_point_on_sphere(radius, rng);
  Random_points_in_sphere_3 random_point_in_sphere(radius, rng);

  // Get some functors
  typename BGT::Construct_segment_3 segment_3 =
                              BGT().construct_segment_3_object();
  typename BGT::Construct_vector_3 vector_3 =
                              BGT().construct_vector_3_object();
  typename BGT::Construct_translated_point_3 translate =
                              BGT().construct_translated_point_3_object();
  typename BGT::Construct_center_3 center = BGT().construct_center_3_object();

  // Get translation from origin to sphere center
  Point_3 center_pt = center(r_domain_.bounding_sphere(r_domain_.bbox_));
  const Vector_3 sphere_translation = vector_3(CGAL::ORIGIN, center_pt);

  // Create nb_point points
  int n = nb_points;
#ifdef CGAL_MESH_3_VERBOSE
  std::cerr << "construct initial points:\n";
#endif
  while ( 0 != n )
  {
    // Get a random segment
    const Point_3 random_point = translate(*random_point_on_sphere,
                                           sphere_translation);
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


}  // end namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // LABELLED_MESH_TRAITS_3_H_
