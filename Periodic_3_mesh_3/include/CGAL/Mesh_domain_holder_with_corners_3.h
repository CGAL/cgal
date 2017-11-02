// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
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
// Author(s)     : Mikhail Bogdanov
//
#ifndef CGAL_MESH_DOMAIN_HOLDER_WITH_CORNERS_3_H
#define CGAL_MESH_DOMAIN_HOLDER_WITH_CORNERS_3_H

#include <CGAL/license/Periodic_3_mesh_3.h>

#include <CGAL/Periodic_3_mesh_3/config.h>

#include <CGAL/enum.h>
#include <CGAL/iterator.h>
#include <CGAL/number_utils.h>

#include <map>
#include <utility>

namespace CGAL {

/**
 * @class Mesh_domain_holder_with_corners_3
 *
 * Implements the `Periodic_3MeshDomain_3` concept through the object owned by this class.
 * Provides functionality of the concept `MeshDomainWithFeatures_3`
 * associated with the corners.
 *
 */
template < typename MeshDomain >
class Mesh_domain_holder_with_corners_3
{
  typedef MeshDomain Base;

public:
  // Index types
  typedef typename Base::Index     Index;
  typedef int                      Curve_index;
  typedef int                      Corner_index;

  typedef typename Base::R         Gt;
  typedef Gt                       R;
  typedef typename Base::Point_3   Point_3;
  typedef typename Base::FT        FT;

  typedef CGAL::Tag_true           Has_features;

  /// Constructors
  Mesh_domain_holder_with_corners_3(const Base& mesh_domain)
    : mesh_domain_(mesh_domain)
  { }

  /// Destructor
  ~Mesh_domain_holder_with_corners_3() { }

  // -----------------------------------
  // Delegation. The holder delegates all the functionality of the Periodic_3MeshDomain_3
  // concept to the member which implements a model of the concept.
  // -----------------------------------

  typedef typename Base::Subdomain        Subdomain;
  typedef typename Base::Surface_patch    Surface_patch;
  typedef typename Base::Intersection     Intersection;

  typedef typename Base::Iso_cuboid_3     Iso_cuboid_3;

  const Iso_cuboid_3& periodic_bounding_box() const
  {
    return mesh_domain_.periodic_bounding_box();
  }

  typedef typename Base::Construct_initial_points Construct_initial_points;

  Construct_initial_points construct_initial_points_object() const
  {
    return mesh_domain_.construct_initial_points_object();
  }

  typedef typename Base::Is_in_domain Is_in_domain;

  Is_in_domain is_in_domain_object() const
  {
    return mesh_domain_.is_in_domain_object();
  }

  typedef typename Base::Do_intersect_surface Do_intersect_surface;

  Do_intersect_surface do_intersect_surface_object() const
  {
    return mesh_domain_.do_intersect_surface_object();
  }

  typedef typename Base::Construct_intersection Construct_intersection;

  Construct_intersection construct_intersection_object() const
  {
    return mesh_domain_.construct_intersection_object();
  }

  typedef typename Base::Surface_patch_index Surface_patch_index;

  Index index_from_surface_patch_index(const Surface_patch_index& index) const
  {
    return mesh_domain_.index_from_surface_patch_index(index);
  }

  typedef typename Base::Subdomain_index Subdomain_index;

  Index index_from_subdomain_index(const Subdomain_index& index) const
  {
    return mesh_domain_.index_from_subdomain_index(index);
  }

  Surface_patch_index surface_patch_index(const Index& index) const
  {
    return mesh_domain_.surface_patch_index(index);
  }

  Subdomain_index subdomain_index(const Index& index) const
  {
    return mesh_domain_.subdomain_index(index);
  }

  // -----------------------------------
  // Backward Compatibility
  // -----------------------------------
#ifndef CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  typedef Surface_patch_index   Surface_index;

  Index index_from_surface_index(const Surface_index& index) const
  {
    return mesh_domain_.index_from_surface_index(index);
  }

  Surface_index surface_index(const Index& index) const
  {
    return mesh_domain_.surface_index(index);
  }
#endif // CGAL_MESH_3_NO_DEPRECATED_SURFACE_INDEX
  // -----------------------------------
  // End backward Compatibility
  // -----------------------------------

  // -----------------------------------
  // End Delegation
  // -----------------------------------

  /// OutputIterator value type is std::pair<Corner_index, Point_3>
  template <typename OutputIterator>
  OutputIterator get_corners(OutputIterator out) const;

  /// OutputIterator value type is CGAL::cpp0x::tuple<Curve_index,
  /// std::pair<Point_3,Index>, std::pair<Point_3,Index> >
  template <typename OutputIterator>
  OutputIterator get_curves(OutputIterator out) const;

  /// Returns the length of the curve segment, on the curve with index
  /// \c  curve_index, from \c p to \c q, in the orientation
  /// \c orientation
  ///
  /// If the curve containing \c p and \c q is a loop,
  /// the orientation gives identifies which portion of the loop
  /// corresponds to the arc, otherwise \c orientation must be compatible
  /// with the orientation of \c p and \c q on the curve segment.
  FT curve_segment_length(const Point_3& p, const Point_3 q,
                          const Curve_index& curve_index,
                          CGAL::Orientation orientation) const;

  /// Returns the length of the curve with index
  /// \c curve_index
  FT curve_length(const Curve_index& curve_index) const;

  /// Returns the geodesic distance between points p and q of curve
  /// \c curve_index
  FT signed_geodesic_distance(const Point_3& p, const Point_3& q,
                              const Curve_index& curve_index) const;

  /// Construct a point on curve \c curve_index at geodesic distance \c distance
  /// of \c starting_point
  Point_3
  construct_point_on_curve(const Point_3& starting_point,
                           const Curve_index& curve_index,
                           FT distance) const;

  /// Returns the sign of the orientation of p,q,r along curve segment
  /// of index \c index
  CGAL::Sign distance_sign_along_loop(const Point_3& p,
                                       const Point_3& q,
                                       const Point_3& r,
                                       const Curve_index& index) const;

  /// Returns an Index from a Curve_index
  Index index_from_Curve_index(const Curve_index& index) const
  { return Index(index); }

  /// Returns an Curve_index from an Index
  Curve_index curve_index(const Index& index) const
  { return boost::get<Curve_index>(index); }

  /// Returns an Index from a Corner_index
  Index index_from_corner_index(const Corner_index& index) const
  { return Index(index); }

  /// Returns an Corner_index from an Index
  Corner_index corner_index(const Index& index) const
  { return boost::get<Corner_index>(index); }

  template <typename InputIterator>
  void add_corners(InputIterator first, InputIterator last);

  void add_corner(const Point_3& point);

  /// Insert one edge into domain
  /// InputIterator value type is Point_3
  template <typename InputIterator>
  Curve_index insert_edge(InputIterator first, InputIterator last);

  /// Returns the sign of the geodesic distance between \c p and \c q
  /// Precondition: index is not a loop
  CGAL::Sign distance_sign(const Point_3& p, const Point_3& q,
                           const Curve_index& index) const;

  /// Returns `true` if curve \c curve_index is a loop
  bool is_loop(const Curve_index& index) const;

  /// Returns `true` if the portion of the curve segment of index \c index,
  /// from \c c1 to \c c2 in the orientation \c orientation, is covered by the spheres of
  /// centers \c c1 and \c c2 and squared radii \c sq_r1 and \c sq_r2
  /// respectively.
  bool is_curve_segment_covered(const Curve_index& index,
                                CGAL::Orientation orientation,
                                const Point_3& c1, const Point_3& c2,
                                const FT sq_r1, const FT sq_r2) const;

  Index index_from_curve_index(const Curve_index& index) const
  { return Index(index); }

private:
  /// Returns Index associated to p (p must be the coordinates of a corner
  /// point)
  Index point_corner_index(const Point_3& p) const;

private:
  typedef std::map<Point_3,Corner_index> Corners;

  Corners corners_;
  const Base& mesh_domain_;

private:
  // Disabled copy constructor & assignment operator
  typedef Mesh_domain_holder_with_corners_3 Self;
  Mesh_domain_holder_with_corners_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Mesh_domain_holder_with_corners_3

template <class MD_>
template <typename OutputIterator>
OutputIterator
Mesh_domain_holder_with_corners_3<MD_>::
get_corners(OutputIterator out) const
{
  for ( typename Corners::const_iterator cit = corners_.begin(),
        end = corners_.end();
        cit != end ; ++cit ) {
    *out++ = std::make_pair(cit->second,cit->first);
  }

  return out;
}

template <class MD_>
template <typename OutputIterator>
OutputIterator
Mesh_domain_holder_with_corners_3<MD_>::
get_curves(OutputIterator out) const
{
  return out;
}

template <class MD_>
typename Mesh_domain_holder_with_corners_3<MD_>::Index
Mesh_domain_holder_with_corners_3<MD_>::
point_corner_index(const Point_3& p) const
{
  typename Corners::const_iterator p_index_it = corners_.find(p);
  if ( p_index_it == corners_.end() )
  {
    CGAL_assertion(false);
    return Index();
  }

  return p_index_it->second;
}

template <class MD_>
typename Mesh_domain_holder_with_corners_3<MD_>::FT
Mesh_domain_holder_with_corners_3<MD_>::
curve_segment_length(const Point_3& /*p*/, const Point_3 /*q*/,
                     const Curve_index& /*curve_index*/,
                     CGAL::Orientation /*orientation*/) const
{
  assert(false);
  return FT();
}

template <class MD_>
typename Mesh_domain_holder_with_corners_3<MD_>::FT
Mesh_domain_holder_with_corners_3<MD_>::
curve_length(const Curve_index& /*curve_index*/) const
{
  assert(false);
  return FT();
}

template <class MD_>
typename Mesh_domain_holder_with_corners_3<MD_>::FT
Mesh_domain_holder_with_corners_3<MD_>::
signed_geodesic_distance(const Point_3& /* p */,
                         const Point_3& /* q */,
                         const Curve_index& /* curve_index */) const
{
  assert(false);
  return FT();
}

template <class MD_>
typename Mesh_domain_holder_with_corners_3<MD_>::Point_3
Mesh_domain_holder_with_corners_3<MD_>::
construct_point_on_curve(const Point_3& /* starting_point */,
                         const Curve_index& /* curve_index */,
                         FT /* distance */) const
{
  assert(false);

  return Point_3();
}

template <class MD_>
template <typename InputIterator>
void
Mesh_domain_holder_with_corners_3<MD_>::
add_corners(InputIterator first, InputIterator last)
{
  while ( first != last )
  {
    add_corner(*first);
    ++first;
  }
}

template <class MD_>
void
Mesh_domain_holder_with_corners_3<MD_>::
add_corner(const Point_3& point)
{
  //#warning modify corner_index!
  const Corner_index corner_index = 0;

  corners_.insert(std::make_pair(point, corner_index));
}

template <class MD_>
CGAL::Sign
Mesh_domain_holder_with_corners_3<MD_>::
distance_sign(const Point_3& /* p */, const Point_3& /* q */,
              const Curve_index& /* index */) const
{
  assert(false);
  return CGAL::Sign();
}

template <class MD_>
bool
Mesh_domain_holder_with_corners_3<MD_>::
is_loop(const Curve_index& /*index*/) const
{
  assert(false);
  return false;
}

template <class MD_>
bool
Mesh_domain_holder_with_corners_3<MD_>::
is_curve_segment_covered(const Curve_index& /*index*/,
                         CGAL::Orientation /*orientation*/,
                         const Point_3& /*c1*/, const Point_3& /*c2*/,
                         const FT /*sq_r1*/, const FT /*sq_r2*/) const
{
  assert(false);
  return false;
}

template <class MD_>
CGAL::Sign
Mesh_domain_holder_with_corners_3<MD_>::
distance_sign_along_loop(const Point_3& /* p */,
                          const Point_3& /* q */,
                          const Point_3& /* r */,
                          const Curve_index& /* index */) const
{
  assert(false);
  return CGAL::Sign();
}

} // namespace CGAL

#endif // CGAL_MESH_DOMAIN_HOLDER_WITH_CORNERS_3_H
