// Copyright (c) 2009-2010 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Stéphane Tayeb, Laurent Rineau
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_MESH_DOMAIN_WITH_POLYLINE_FEATURES_3_H
#define CGAL_MESH_DOMAIN_WITH_POLYLINE_FEATURES_3_H

#include <CGAL/license/Mesh_3.h>


#include <CGAL/iterator.h>
#include <CGAL/enum.h>
#include <CGAL/number_utils.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/is_streamable.h>
#include <CGAL/Real_timer.h>
#include <CGAL/property_map.h>
#include <CGAL/SMDS_3/internal/indices_management.h>

#include <vector>
#include <set>
#include <map>
#include <algorithm>
#include <type_traits>

#include <variant>
#include <memory>
#include <deque>

namespace CGAL {

/// @cond CGAL_DOCUMENT_INTERNALS
namespace Mesh_3 {
namespace internal {

template <typename Kernel>
class Polyline
{
  typedef typename Kernel::Point_3  Point_3;
  typedef typename Kernel::Segment_3 Segment_3;
  typedef typename Kernel::FT       FT;

  typedef std::vector<Point_3>      Data;

public:
  typedef typename Data::const_iterator const_iterator;

  Polyline() {}
  ~Polyline() {}

  /// adds a point at the end of the polyline
  void add_point(const Point_3& p)
  {
    if( points_.empty() || p != end_point() ) {
      points_.push_back(p);
    }
  }

  /// returns the starting point of the polyline
  const Point_3& start_point() const
  {
    CGAL_assertion( ! points_.empty() );
    return points_.front();
  }

  /// returns the ending point of the polyline
  const Point_3& end_point() const
  {
    CGAL_assertion( ! points_.empty() );
    return points_.back();
  }

  /// returns `true` if the polyline is not degenerate
  bool is_valid() const
  {
    return points_.size() > 1;
  }

  /// returns `true` if polyline is a loop
  bool is_loop() const
  {
    return start_point() == end_point();
  }

  const_iterator next(const_iterator it, Orientation orientation) const {
    if(orientation == POSITIVE) {
      CGAL_assertion(it != (points_.end() - 1));
      if(it == (points_.end() - 2)) {
        CGAL_assertion(is_loop());
        it = points_.begin();
      } else {
        ++it;
      }
    } else {
      CGAL_assertion(orientation == NEGATIVE);
      CGAL_assertion(it != points_.begin());
      if(it == (points_.begin() + 1)) {
        CGAL_assertion(is_loop());
        it = points_.end() - 1;
      } else {
        --it;
      }
    }
    return it;
  }

  bool is_curve_segment_covered(CGAL::Orientation orientation,
                                const Point_3& c1, const Point_3& c2,
                                const FT sq_r1, const FT sq_r2) const
  {
    CGAL_assertion(orientation != CGAL::ZERO);
    typename Kernel::Has_on_bounded_side_3 cover_pred =
      Kernel().has_on_bounded_side_3_object();

    typedef typename Kernel::Sphere_3 Sphere_3;
    const Sphere_3 s1(c1, sq_r1);
    const Sphere_3 s2(c2, sq_r2);

    const_iterator c1_it = locate(c1);
    const_iterator c2_it = locate(c2);

    if(orientation == CGAL::NEGATIVE) {
      ++c1_it;
      ++c2_it;
      CGAL_assertion(c1_it != points_.end());
      CGAL_assertion(c2_it != points_.end());
    }

    if(c1_it == c2_it) return cover_pred(s1, s2, c1, c2);
    const_iterator next_it = this->next(c1_it, orientation);

    if(!cover_pred(s1, s2, c1, *next_it)) return false;

    for(const_iterator it = next_it; it != c2_it; /* in body */) {
      next_it = this->next(it, orientation);
      if(!cover_pred(s1, s2, *it, *next_it)) return false;
      it = next_it;
    } // end loop ]c1_it, c2_it[

    return cover_pred(s1, s2, *c2_it, c2);
  }

  FT curve_segment_length(const Point_3& p, const Point_3 q,
                          CGAL::Orientation orientation) const
  {
    CGAL_assertion(orientation != CGAL::ZERO);
    const_iterator p_it = locate(p);
    const_iterator q_it = locate(q);
    return curve_segment_length(p, q, orientation, p_it, q_it);
  }

  FT curve_segment_length(const Point_3& p, const Point_3 q,
                          CGAL::Orientation orientation,
                          const_iterator p_it,
                          const_iterator q_it) const
  {
    CGAL_assertion(orientation != CGAL::ZERO);

    if(p_it == q_it) {
      const CGAL::Comparison_result cmp = compare_distance(*p_it,p,q);
      if( (cmp != LARGER  && orientation == POSITIVE) ||
          (cmp != SMALLER && orientation == NEGATIVE) )
      {
        // If the orientation of `p` and `q` on the segment is compatible
        // with `orientation`, then return the distance between the two
        // points.
        return distance(p, q);
      }
    }

    if(orientation == CGAL::NEGATIVE) {
      ++p_it;
      ++q_it;
      CGAL_assertion(p_it != points_.end());
      CGAL_assertion(q_it != points_.end());
    }

    const_iterator next_it = this->next(p_it, orientation);
    FT result = distance(p, *next_it);
    for(const_iterator it = next_it; it != q_it; /* in body */) {
      next_it = this->next(it, orientation);
      result += distance(*it, *next_it);
      it = next_it;
    } // end loop ]p_it, q_it[
    result += distance(*q_it, q);
    return result;
  }


  /// returns the angle at the first point.
  /// \pre The polyline must be a loop.
  Angle angle_at_first_point() const {
    CGAL_precondition(is_loop());
    const Point_3& first = points_.front();
    const Point_3& next_p = points_[1];
    const Point_3& prev = points_[points_.size() - 2];
    return angle(prev, first, next_p);
  }

  /// returns the length of the polyline
  FT length() const
  {
    //TODO: cache result
    FT result (0);
    const_iterator it = points_.begin();
    const_iterator previous = it++;

    for ( const_iterator end = points_.end() ; it != end ; ++it, ++previous )
    {
      result += distance(*previous,*it);
    }

    return result;
  }

  /// returns the signed geodesic distance between `p` and `q`.
  FT signed_geodesic_distance(const Point_3& p, const Point_3& q) const
  {
    // Locate p & q on polyline
    const_iterator pit = locate(p);
    const_iterator qit = locate(q,false);

    // If p and q are in the same segment of the polyline
    if ( pit == qit )
    {
      FT result = distance(p,q);

      // Find the closest point to *pit
      if ( compare_distance(*pit,p,q) != CGAL::LARGER )
      { return result; }
      else
      { return -result; }
    }
    if(is_loop()) {
      const FT positive_distance = curve_segment_length(p, q, CGAL::POSITIVE, pit, qit);
      const FT negative_distance = curve_segment_length(p, q, CGAL::NEGATIVE, pit, qit);
      return (positive_distance < negative_distance)
        ?    positive_distance
        : (- negative_distance);
    } else {
      return (pit <= qit)
        ?     curve_segment_length(p, q, CGAL::POSITIVE, pit, qit)
        : ( - curve_segment_length(p, q, CGAL::NEGATIVE, pit, qit) );
    }
  }


  /// returns a point at geodesic distance `distance` from p along the
  /// polyline. The polyline is oriented from starting point to end point.
  /// The distance could be negative.
  Point_3 point_at(const Point_3& p, FT distance) const
  {
    // use first point of the polyline instead of p
    distance += curve_segment_length(start_point(),p,CGAL::POSITIVE);

    // If polyline is a loop, ensure that distance is given from start_point()
    if ( is_loop() )
    {
      if ( distance < FT(0) ) { distance += length(); }
      else if ( distance > length() ) { distance -= length(); }
    }

    CGAL_assertion( distance >= FT(0) );
    CGAL_assertion( distance <= length() );

    // Initialize iterators
    const_iterator pit = points_.begin();
    const_iterator previous = pit++;

    // Iterate to find which segment contains the point we want to construct
    FT segment_length = this->distance(*previous,*pit);
    while ( distance > segment_length )
    {
      distance -= segment_length;

      // Increment iterators and update length
      ++previous;
      ++pit;

      if (pit == points_.end())
        return *previous;

      segment_length = this->distance(*previous,*pit);
    }

    // return point at distance from current segment source
    typedef typename Kernel::Vector_3 Vector_3;
    Vector_3 v (*previous, *pit);

    return (*previous) + (distance / CGAL::sqrt(v.squared_length())) * v;
  }

  bool are_ordered_along(const Point_3& p, const Point_3& q) const
  {
    CGAL_precondition(!is_loop());

    // Locate p & q on polyline
    const_iterator pit = locate(p);
    const_iterator qit = locate(q,true);

    // Points are not located on the same segment
    if ( pit != qit ) { return (pit <= qit); }

    // pit == qit, then we have to sort p&q along (pit,pit+1)
    return ( compare_distance(*pit,p,q) != CGAL::LARGER );
  }

private:
  const_iterator first_segment_source() const
  {
    CGAL_precondition(is_valid());
    return points_.begin();
  }

  const_iterator last_segment_source() const
  {
    CGAL_precondition(is_valid());
    return (points_.end() - 2);
  }

  /// returns an iterator on the starting point of the segment of the
  /// polyline which contains p
  /// if end_point_first is true, then --end is returned instead of begin
  /// if p is the starting point of a loop.
  const_iterator locate(const Point_3& p, bool end_point_first=false) const
  {
    CGAL_precondition(is_valid());

    // First look if p is one of the points of the polyline
    const_iterator result = std::find(points_.begin(), points_.end(), p);
    if ( result != points_.end() )
    {
      if ( result != points_.begin() )
      { return --result; }
      else
      {
        // Treat loops
        if ( end_point_first && p == end_point() )
        { return last_segment_source(); }
        else
        { return result; }
      }
    }

    CGAL_assertion(result == points_.end());

    // Get result by projecting p on the polyline
    const_iterator it = points_.begin();
    const_iterator previous = it;
    Segment_3 nearest_segment;
    const_iterator nearest_vertex = it;
    result = nearest_vertex;
    bool nearest_is_a_segment = false;

    while ( ++it != points_.end() )
    {
      Segment_3 seg (*previous, *it);

      if(nearest_is_a_segment)
      {
        if(compare_distance(p, *it, nearest_segment) == CGAL::SMALLER)
        {
          nearest_vertex = it;
          nearest_is_a_segment = false;
          result = it;
          if (possibly(angle(*previous, *it, p) == CGAL::ACUTE) &&
              compare_distance(p, seg, *nearest_vertex) == CGAL::SMALLER)
          {
            nearest_segment = seg;
            nearest_is_a_segment = true;
            result = previous;
          }
        }
        else if(compare_distance(p, seg, nearest_segment) == CGAL::SMALLER)
        {
          nearest_segment = seg;
          result = previous;
        }
      }
      else {
        if(compare_distance(p, *it, *nearest_vertex) == CGAL::SMALLER)
        {
          nearest_vertex = it;
          result = it;
        }
        if ((nearest_vertex != it ||
             possibly(angle(*previous, *it, p) == CGAL::ACUTE)) &&
            compare_distance(p, seg, *nearest_vertex) == CGAL::SMALLER)
        {
          nearest_segment = seg;
          nearest_is_a_segment = true;
          result = previous;
        }
      }
      previous = it;
    } // end the while loop on the vertices of the polyline


    if(result == points_.begin()) {
      return (end_point_first && !nearest_is_a_segment) ? last_segment_source() : points_.begin();
    } else {
      return result;
    }
  }

  // FT squared_distance(const Point_3& p, const Point_3& q) const
  // {
  //   typename Kernel::Compute_squared_distance_3 sq_distance =
  //     Kernel().compute_squared_distance_3_object();
  //   return sq_distance(p,q);
  // }

  FT distance(const Point_3& p, const Point_3& q) const
  {
    return CGAL::sqrt(squared_distance(p, q));
  }

  Angle angle(const Point_3& p,
              const Point_3& angle_vertex_point,
              const Point_3& q) const
  {
    typename Kernel::Angle_3 compute_angle =  Kernel().angle_3_object();
    return compute_angle(p,angle_vertex_point,q);
  }

  template <typename T1, typename T2>
  CGAL::Sign compare_distance(const Point_3& p,
                              const T1& obj1,
                              const T2& obj2) const
  {
    typename Kernel::Compare_distance_3 compare_distance =
      Kernel().compare_distance_3_object();
    return compare_distance(p,obj1,obj2);
  }

public:
  Data points_;
}; // end class Polyline


template <typename GT, typename MapIterator>
struct Mesh_domain_segment_of_curve_primitive{
  typedef typename std::iterator_traits<MapIterator>::value_type Polyline;
  typedef typename Polyline::const_iterator PointIterator;

  typedef std::pair<Curve_index, PointIterator> Id;

  typedef typename std::iterator_traits<PointIterator>::value_type Point;

  typedef typename GT::Segment_3 Datum;

  Id id_;

  Mesh_domain_segment_of_curve_primitive(Id id) : id_(id) {}

  const Id& id() const { return id_; }

  const Point& reference_point() const {
    return *(id_.second);
  }

  Datum datum() const {
    return Datum(*id_.second, *(id_.second+1));
  }
}; // end Mesh_domain_segment_of_curve_primitive

template <typename MDwPF, bool patch_id_is_streamable>
struct Display_incidences_to_patches_aux {
  template <typename Container, typename Point>
  void operator()(std::ostream& os, Point p, typename MDwPF::Curve_index id,
                  const Container&) const;
};

template <typename MDwPF> //specialization when patch_id_is_streamable == false
struct Display_incidences_to_patches_aux<MDwPF, false> {
  template <typename Container, typename Point>
  void operator()(std::ostream& os, Point p,
                  typename MDwPF::Curve_index id,
                  const Container&) const;
};

template <typename MDwPF, bool curve_id_is_streamable>
struct Display_incidences_to_curves_aux {
  template <typename Container, typename Point>
  void operator()(std::ostream& os, Point p, typename MDwPF::Curve_index id,
                  const Container&) const;
};

template <typename MDwPF> //specialization when curve_id_is_streamable == false
struct Display_incidences_to_curves_aux<MDwPF, false> {
  template <typename Container, typename Point>
  void operator()(std::ostream& os, Point p,  typename MDwPF::Curve_index id,
                  const Container&) const;
};

} // end of namespace CGAL::Mesh_3::internal
} // end of namespace CGAL::Mesh_3
/// @endcond

/*!
\ingroup PkgMesh3Domains

The class `Mesh_domain_with_polyline_features_3` enables the user
to add some 0- and 1-dimensional
features into any model of the `MeshDomain_3` concept.
The 1-dimensional features are described as polylines
whose endpoints are the added corners.

\tparam MD is the type of the domain which is extended. It has to be a model of the `MeshDomain_3` concept.

\cgalModels{MeshDomainWithFeatures_3}

\sa `MeshPolyline_3`
\sa `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT>`
*/
template < typename MD >
class Mesh_domain_with_polyline_features_3
  : public MD
{
  typedef Mesh_domain_with_polyline_features_3<MD>   Self;

public:
  /// \name Types
  /// @{

  typedef typename MD::Surface_patch_index           Surface_patch_index;
  typedef typename MD::Subdomain_index               Subdomain_index;
  typedef int                                        Curve_index;
  typedef int                                        Corner_index;

#ifdef DOXYGEN_RUNNING
  typedef unspecified_type                           Index;
#else
  typedef typename Mesh_3::internal::Index_generator_with_features<
    typename MD::Subdomain_index,
    Surface_patch_index,
    Curve_index,
    Corner_index>::type                              Index;
#endif

  typedef CGAL::Tag_true                             Has_features;
  typedef typename MD::R::FT                         FT;

  /// @}

#ifndef CGAL_NO_DEPRECATED_CODE
  typedef Curve_index                                Curve_segment_index;
#endif

  typedef typename MD::R                             GT;
  typedef GT                                         R;
  typedef typename MD::Point_3                       Point_3;

  /// \name Creation
  /// @{

  // forwards the arguments to the constructor of the base class.
  template <typename ... T>
  Mesh_domain_with_polyline_features_3(const T& ...o)
    : MD(o...)
    , current_corner_index_(1)
    , current_curve_index_(1)
    , curves_aabb_tree_is_built(false) {}

  Mesh_domain_with_polyline_features_3(const Mesh_domain_with_polyline_features_3&) = default;

  /// @}

  /// \name Operations
  /// @{

  /// @cond CGAL_DOCUMENT_INTERNALS

  /// adds a 0-dimensional feature in the domain.
  Corner_index add_corner(const Point_3& p);

  /// Overload where the last parameter `out` is not `CGAL::Emptyset_iterator()`.
  template <typename InputIterator, typename IndicesOutputIterator>
  IndicesOutputIterator
  add_corners(InputIterator first, InputIterator end,
              IndicesOutputIterator out  /*= CGAL::Emptyset_iterator()*/);

  /*!
    adds 0-dimensional features in the domain.

    The value type of `InputIterator` must be `Point_3`.
  */
  template <typename InputIterator>
  void
  add_corners(InputIterator first, InputIterator end)
  { add_corners(first, end, CGAL::Emptyset_iterator()); }

  Corner_index register_corner(const Point_3& p, const Curve_index& index);
  Corner_index add_corner_with_context(const Point_3& p, const Surface_patch_index& index);

  /// Overload where the last parameter `out` is not
  /// `CGAL::Emptyset_iterator()`.
  template <typename InputIterator, typename IndicesOutputIterator>
  IndicesOutputIterator
  add_features(InputIterator first, InputIterator end,
               IndicesOutputIterator out /*= CGAL::Emptyset_iterator()*/);

  template <typename InputIterator,
            typename PolylinePMap,
            typename IncidentPatchesIndicesPMap,
            typename IndicesOutputIterator>
  IndicesOutputIterator
  add_features_and_incidences
  (InputIterator first, InputIterator end,
   PolylinePMap polyline_pmap,
   IncidentPatchesIndicesPMap incident_patches_indices_pmap,
   IndicesOutputIterator out /* = CGAL::Emptyset_iterator() */);

  template <typename InputIterator, typename IndicesOutputIterator>
  IndicesOutputIterator
  add_features_with_context(InputIterator first, InputIterator end,
                            IndicesOutputIterator out /*=
                                                        CGAL::Emptyset_iterator()*/);

  /// @endcond

  /*!
    adds 1-dimensional features in the domain.

    The value type of `InputIterator` must be a model of the concept `MeshPolyline_3`.
  */
  template <typename InputIterator>
  void
  add_features(InputIterator first, InputIterator end)
  { add_features(first, end, CGAL::Emptyset_iterator()); }

  /// @cond CGAL_DOCUMENT_INTERNALS

  /// Undocumented function, kept for backward-compatibility with existing code
  template <typename InputIterator>
  void
  add_features_with_context(InputIterator first, InputIterator end)
  { add_features_with_context(first, end, CGAL::Emptyset_iterator()); }

  /// @endcond

  /*!
    adds 1-dimensional features (curves) from the range `[first, end)` in the domain with their incidences
    with 2-dimensional features (patches) of the domain.

    \tparam InputIterator input iterator over curves
    \tparam PolylinePMap is a model of `ReadablePropertyMap` with key type
      `std::iterator_traits<InputIterator>::%reference` and a value type
      that is a model of `MeshPolyline_3`.
    \tparam IncidentPatchesIndicesPMap is a model of `ReadablePropertyMap`
      with key type `std::iterator_traits<InputIterator>::%reference` and a
      value type that is a range of `Surface_patch_index`.

    \param first iterator to the first curve of the sequence
    \param end past-the-end iterator of the sequence of curves
    \param polyline_pmap the property map that provides access to the
      polyline, model of `MeshPolyline_3`, from the `%reference` type of
      the iterator
    \param incident_patches_indices_pmap the property map that provides
      access to the set of indices of the surface patches that are incident to
      a given 1D-feature (curve)
  */
  template <typename InputIterator,
            typename PolylinePMap,
            typename IncidentPatchesIndicesPMap>
  void
  add_features_and_incidences
  (InputIterator first, InputIterator end,
   PolylinePMap polyline_pmap,
   IncidentPatchesIndicesPMap incident_patches_indices_pmap)
  {
    add_features_and_incidences(first, end, polyline_pmap,
                                incident_patches_indices_pmap,
                                CGAL::Emptyset_iterator());
  }

  /// @}

  /// \name Implementation of the concept MeshDomainWithFeatures_3
  /// The following methods implement the requirements of the concept
  /// `MeshDomainWithFeatures_3`.
  /// @{

  /// implements `MeshDomainWithFeatures_3::get_corners()`.
  /// OutputIterator is `std::pair<Corner_index, Point_3>`
  template <typename OutputIterator>
  OutputIterator get_corners(OutputIterator out) const;

  /// implements `MeshDomainWithFeatures_3::get_curves()`.
  /// OutputIterator value type is std::tuple<Curve_index,
  /// std::pair<Point_3,Index>, std::pair<Point_3,Index> >
  template <typename OutputIterator>
  OutputIterator get_curves(OutputIterator out) const;

  /// implements `MeshDomainWithFeatures_3::curve_segment_length()`.
  FT curve_segment_length(const Point_3& p, const Point_3 q,
                          const Curve_index& curve_index,
                          CGAL::Orientation orientation) const;

  /// implements `MeshDomainWithFeatures_3::curve_length()`.
  FT curve_length(const Curve_index& curve_index) const;

  /// implements `MeshDomainWithFeatures_3::construct_point_on_curve()`.
  Point_3
  construct_point_on_curve(const Point_3& starting_point,
                           const Curve_index& curve_index,
                           FT distance) const;
  /// implements `MeshDomainWithFeatures_3::distance_sign_along_loop()`.
  CGAL::Sign distance_sign_along_loop(const Point_3& p,
                                      const Point_3& q,
                                      const Point_3& r,
                                      const Curve_index& index) const;

  /// implements `MeshDomainWithFeatures_3::distance_sign()`.
  CGAL::Sign distance_sign(const Point_3& p, const Point_3& q,
                           const Curve_index& index) const;

  /// implements `MeshDomainWithFeatures_3::is_loop()`.
  bool is_loop(const Curve_index& index) const;

  /// implements `MeshDomainWithFeatures_3::is_curve_segment_covered()`.
  bool is_curve_segment_covered(const Curve_index& index,
                                CGAL::Orientation orientation,
                                const Point_3& c1, const Point_3& c2,
                                const FT sq_r1, const FT sq_r2) const;

  /**
   * Returns the index to be stored in a vertex lying on the surface identified
   * by `index`.
   */
  Index index_from_surface_patch_index(const Surface_patch_index& index) const
  { return Index(index); }

  /**
   * Returns the index to be stored in a vertex lying in the subdomain
   * identified by `index`.
   */
  Index index_from_subdomain_index(const Subdomain_index& index) const
  { return Index(index); }

  /// returns an `Index` from a `Curve_index`
  Index index_from_curve_index(const Curve_index& index) const
  { return Index(index); }

  /// returns an `Index` from a `Corner_index`
  Index index_from_corner_index(const Corner_index& index) const
  { return Index(index); }

  /**
   * Returns the `Surface_patch_index` of the surface patch
   * where lies a vertex with dimension 2 and index `index`.
   */
  Surface_patch_index surface_patch_index(const Index& index) const
  { return Mesh_3::internal::get_index<Surface_patch_index>(index); }

  /**
   * Returns the index of the subdomain containing a vertex
   *  with dimension 3 and index `index`.
   */
  Subdomain_index subdomain_index(const Index& index) const
  { return Mesh_3::internal::get_index<Subdomain_index>(index); }

  /// returns a `Curve_index` from an `Index`
  Curve_index curve_index(const Index& index) const
  { return Mesh_3::internal::get_index<Curve_index>(index); }

  /// returns a `Corner_index` from an `Index`
  Corner_index corner_index(const Index& index) const
  { return Mesh_3::internal::get_index<Corner_index>(index); }

  /// @cond CGAL_DOCUMENT_INTERNALS
#ifndef CGAL_NO_DEPRECATED_CODE
  CGAL_DEPRECATED_MSG("deprecated: use curve_index() instead")
  Curve_index curve_segment_index(const Index& index) const {
    return curve_index(index);
  }
#endif // CGAL_NO_DEPRECATED_CODE

  FT signed_geodesic_distance(const Point_3& p, const Point_3& q,
                              const Curve_index& curve_index) const;

  template <typename Surf_p_index, typename IncidenceMap>
  void reindex_patches(const std::vector<Surf_p_index>& map, IncidenceMap& incidence_map);

  template <typename Surf_p_index>
  void reindex_patches(const std::vector<Surf_p_index>& map);

  template <typename IndicesOutputIterator>
  IndicesOutputIterator
  get_incidences(Curve_index id, IndicesOutputIterator out) const;

  template <typename IndicesOutputIterator>
  IndicesOutputIterator
  get_corner_incidences(Corner_index id, IndicesOutputIterator out) const;

  template <typename IndicesOutputIterator>
  IndicesOutputIterator
  get_corner_incident_curves(Corner_index id, IndicesOutputIterator out) const;

  typedef std::set<Surface_patch_index> Surface_patch_index_set;

  const Surface_patch_index_set&
  get_incidences(Curve_index id) const;

  void display_corner_incidences(std::ostream& os, Point_3, Corner_index id);

  /// Insert one edge into domain
  /// InputIterator value type is Point_3
  template <typename InputIterator>
  Curve_index insert_edge(InputIterator first, InputIterator end);
  /// @endcond

  /// @}

private:
  void compute_corners_incidences();

  /// returns Index associated to p (p must be the coordinates of a corner
  /// point)
  Index point_corner_index(const Point_3& p) const;

private:
  typedef std::map<Point_3,Corner_index> Corners;

  typedef Mesh_3::internal::Polyline<GT> Polyline;   // defined before curve accessors
  typedef std::deque<Polyline> Edges;
  typedef std::deque<Surface_patch_index_set> Edges_incidences;
  typedef std::deque<std::set<Curve_index>> Corners_tmp_incidences;
  typedef std::map<Corner_index, Surface_patch_index_set> Corners_incidences;

  typedef Mesh_3::internal::Mesh_domain_segment_of_curve_primitive<
    GT, typename Edges::const_iterator> Curves_primitives;
  typedef CGAL::AABB_traits_3<GT, Curves_primitives> AABB_curves_traits;

  Corners corners_;
  Corners_tmp_incidences corners_tmp_incidences_;
  Corner_index current_corner_index_;
  Corners_incidences corners_incidences_;

  Edges edges_;
  Curve_index current_curve_index_;
  Edges_incidences edges_incidences_;

  // --- curve accessors (now placed after all definitions) ---
  Polyline& curve(Curve_index index) {
    CGAL_precondition(index > 0 && std::size_t(index) <= edges_.size());
    return edges_[index - 1];
  }
  const Polyline& curve(Curve_index index) const {
    CGAL_precondition(index > 0 && std::size_t(index) <= edges_.size());
    return edges_[index - 1];
  }
  // ----------------------------------------------------------

public:
  /// @cond CGAL_DOCUMENT_INTERNALS
  typedef CGAL::AABB_tree<AABB_curves_traits> Curves_AABB_tree;

private:
  mutable std::shared_ptr<Curves_AABB_tree> curves_aabb_tree_ptr_;
  mutable bool curves_aabb_tree_is_built;

public:
  const Corners_incidences& corners_incidences_map() const
  { return corners_incidences_; }

  const Curves_AABB_tree& curves_aabb_tree() const {
    if(!curves_aabb_tree_is_built) build_curves_aabb_tree();
    return *curves_aabb_tree_ptr_;
  }
  Curve_index maximal_curve_index() const {
    if(edges_incidences_.empty()) return Curve_index();
    return static_cast<Curve_index>(edges_incidences_.size());
  }

  void build_curves_aabb_tree() const {
#ifdef CGAL_MESH_3_VERBOSE
    std::cerr << "Building curves AABB tree...";
    CGAL::Real_timer timer;
    timer.start();
#endif
    if(curves_aabb_tree_ptr_) {
      curves_aabb_tree_ptr_->clear();
    } else {
      curves_aabb_tree_ptr_ = std::make_shared<Curves_AABB_tree>();
    }
    for(typename Edges::const_iterator
          edges_it = edges_.begin(),
          edges_end = edges_.end();
        edges_it != edges_end; ++edges_it)
    {
      const Polyline& polyline = *edges_it;
      Curve_index curve_id = 1 + static_cast<Curve_index>(edges_it - edges_.begin());
      for(typename Polyline::const_iterator
            pit = polyline.points_.begin(),
            end = polyline.points_.end() - 1;
          pit != end; ++pit)
      {
        curves_aabb_tree_ptr_->insert(std::make_pair(curve_id, pit));
      }
    }
    curves_aabb_tree_ptr_->build();
    curves_aabb_tree_is_built = true;
#ifdef CGAL_MESH_3_VERBOSE
    timer.stop();
    std::cerr << " done (" << timer.time() * 1000 << " ms)" << std::endl;
#endif
  }

  /// @endcond

}; // class Mesh_domain_with_polyline_features_3



template <class MD_>
template <typename OutputIterator>
OutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
get_corners(OutputIterator out) const
{
  for ( typename Corners::const_iterator
       cit = corners_.begin(), end = corners_.end() ; cit != end ; ++cit )
  {
    *out++ = std::make_pair(cit->second,cit->first);
  }
  return out;
}

template <class MD_>
template <typename OutputIterator>
OutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
get_curves(OutputIterator out) const
{
  Curve_index idx = 1;
  for ( typename Edges::const_iterator
       eit = edges_.begin(), end = edges_.end() ; eit != end ; ++eit, ++idx )
  {
    CGAL_assertion( eit->is_valid() );

    const Point_3& p = eit->start_point();
    const Point_3& q = eit->end_point();

    Index p_index, q_index;
    if ( ! eit->is_loop() )
    {
      p_index = point_corner_index(p);
      q_index = point_corner_index(q);
    }
    else
    {
      p_index = index_from_curve_index(idx);
      q_index = p_index;
    }

    *out++ = {idx,
              std::make_pair(p,p_index),
              std::make_pair(q,q_index)};
  }
  return out;
}


template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Index
Mesh_domain_with_polyline_features_3<MD_>::
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
typename Mesh_domain_with_polyline_features_3<MD_>::FT
Mesh_domain_with_polyline_features_3<MD_>::
curve_segment_length(const Point_3& p, const Point_3 q,
                     const Curve_index& curve_index,
                     CGAL::Orientation orientation) const
{
  CGAL_assertion(curve_index > 0 &&
                 std::size_t(curve_index) <= edges_.size());
  const Polyline& poly = curve(curve_index);
  return poly.curve_segment_length(p, q, orientation);
}


template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::FT
Mesh_domain_with_polyline_features_3<MD_>::
curve_length(const Curve_index& curve_index) const
{
  CGAL_assertion(curve_index > 0 &&
                 std::size_t(curve_index) <= edges_.size());
  return curve(curve_index).length();
}


template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Point_3
Mesh_domain_with_polyline_features_3<MD_>::
construct_point_on_curve(const Point_3& starting_point,
                         const Curve_index& curve_index,
                         FT distance) const
{
  CGAL_assertion(curve_index > 0 &&
                 std::size_t(curve_index) <= edges_.size());
  const Polyline& poly = curve(curve_index);
  return poly.point_at(starting_point, distance);
}


/// @cond CGAL_DOCUMENT_INTERNALS
template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Corner_index
Mesh_domain_with_polyline_features_3<MD_>::
add_corner(const Point_3& p)
{
  typename Corners::iterator cit = corners_.lower_bound(p);

  if(cit != corners_.end() && !(corners_.key_comp()(p, cit->first)))
    return cit->second;

  const Corner_index index = current_corner_index_++;
  corners_.insert(cit, std::make_pair(p, index));

  corners_tmp_incidences_.emplace_back();
  corners_incidences_.insert(std::make_pair(index, Surface_patch_index_set()));

  return index;
}


template <class MD_>
template <typename InputIterator, typename IndicesOutputIterator>
IndicesOutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
add_corners(InputIterator first, InputIterator end,
            IndicesOutputIterator indices_out)
{
  while ( first != end )
    *indices_out++ = add_corner(*first++);
  return indices_out;
}

template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Corner_index
Mesh_domain_with_polyline_features_3<MD_>::
register_corner(const Point_3& p, const Curve_index& curve_index)
{
  Corner_index index = add_corner(p);
  corners_tmp_incidences_[index - 1].insert(curve_index);
  return index;
}


template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Corner_index
Mesh_domain_with_polyline_features_3<MD_>::
add_corner_with_context(const Point_3& p, const Surface_patch_index& surface_patch_index)
{
  Corner_index index = add_corner(p);
  Surface_patch_index_set& incidences = corners_incidences_[index];
  incidences.insert(surface_patch_index);
  return index;
}
/// @endcond


template <class MD_>
template <typename InputIterator, typename IndicesOutputIterator>
IndicesOutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
add_features(InputIterator first, InputIterator end,
             IndicesOutputIterator indices_out)
{
  while ( first != end )
  {
    *indices_out++ = insert_edge(first->begin(), first->end());
    ++first;
  }
  compute_corners_incidences();
  return indices_out;
}

/// @cond CGAL_DOCUMENT_INTERNALS
namespace details {

template <typename PolylineWithContext>
struct Get_content_from_polyline_with_context
{
  typedef Get_content_from_polyline_with_context Self;
  typedef PolylineWithContext key_type;
  typedef typename PolylineWithContext::Bare_polyline value_type;
  typedef const value_type& reference;
  typedef boost::readable_property_map_tag category;

  friend reference get(const Self&, const key_type& polyline) {
    return polyline.polyline_content;
  }
};

template <typename PolylineWithContext>
struct Get_patches_id_from_polyline_with_context
{
  typedef Get_patches_id_from_polyline_with_context Self;
  typedef PolylineWithContext key_type;
  typedef typename PolylineWithContext::Context::Patches_ids value_type;
  typedef const value_type& reference;
  typedef boost::readable_property_map_tag category;

  friend reference get(const Self&, const key_type& polyline) {
    return polyline.context.adjacent_patches_ids;
  }
};

} // end namespace details
/// @endcond

template <class MD_>
template <typename InputIterator,
          typename PolylinePMap,
          typename IncidentPatchesIndicesPMap,
          typename IndicesOutputIterator>
IndicesOutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
add_features_and_incidences(InputIterator first, InputIterator end,
                            PolylinePMap polyline_pmap,
                            IncidentPatchesIndicesPMap inc_patches_ind_pmap,
                            IndicesOutputIterator indices_out)
{
  for( ; first != end ; ++first )
  {
    const typename boost::property_traits<PolylinePMap>::reference
      polyline = get(polyline_pmap, *first);
    const typename boost::property_traits<IncidentPatchesIndicesPMap>::reference
      patches_ids = get(inc_patches_ind_pmap, *first);

    Curve_index curve_id = insert_edge(polyline.begin(), polyline.end());
    edges_incidences_[curve_id - 1].insert(patches_ids.begin(), patches_ids.end());
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    std::cerr << "Curve #" << curve_id << " is incident to the following patches: {";
    for(auto id: patches_ids) {
      std::cerr << " " << id;
    }
    std::cerr << "}\n";
#endif
    *indices_out++ = curve_id;
  }

  compute_corners_incidences();
  return indices_out;
}

/// @cond CGAL_DOCUMENT_INTERNALS
template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::FT
Mesh_domain_with_polyline_features_3<MD_>::
signed_geodesic_distance(const Point_3& p, const Point_3& q,
                         const Curve_index& curve_index) const
{
  CGAL_assertion(curve_index > 0 &&
                 std::size_t(curve_index) <= edges_.size());
  const Polyline& poly = curve(curve_index);
  return poly.signed_geodesic_distance(p,q);
}


template <class MD_>
template <typename InputIterator, typename IndicesOutputIterator>
IndicesOutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
add_features_with_context(InputIterator first, InputIterator end,
                          IndicesOutputIterator indices_out)
{
  typedef typename std::iterator_traits<InputIterator>::value_type Pwc;
  return add_features_and_incidences
    (first, end,
     details::Get_content_from_polyline_with_context<Pwc>(),
     details::Get_patches_id_from_polyline_with_context<Pwc>(),
     indices_out);
}

template <class MD_>
template <typename Surf_p_index, typename IncidenceMap>
void
Mesh_domain_with_polyline_features_3<MD_>::
reindex_patches(const std::vector<Surf_p_index>& map,
                IncidenceMap& incidence_map)
{
  for(auto& patch_index_set : incidence_map)
  {
    Surface_patch_index_set new_index_set;
    for(Surface_patch_index idx : patch_index_set)
    {
      CGAL_assertion(std::size_t(idx) < map.size());
      new_index_set.insert(map[idx]);
    }
    patch_index_set = std::move(new_index_set);
  }
}

template <class MD_>
template <typename Surf_p_index>
void
Mesh_domain_with_polyline_features_3<MD_>::
reindex_patches(const std::vector<Surf_p_index>& map,
                std::map<Corner_index, Surface_patch_index_set>& incidence_map)
{
  for(auto& pair : incidence_map)
  {
    Surface_patch_index_set& patch_index_set = pair.second;
    Surface_patch_index_set new_index_set;
    for(Surface_patch_index idx : patch_index_set)
    {
      CGAL_assertion(std::size_t(idx) < map.size());
      new_index_set.insert(map[idx]);
    }
    patch_index_set = std::move(new_index_set);
  }
}

template <class MD_>
template <typename Surf_p_index>
void
Mesh_domain_with_polyline_features_3<MD_>::
reindex_patches(const std::vector<Surf_p_index>& map)
{
  reindex_patches(map, edges_incidences_);
  reindex_patches(map, corners_incidences_);
}

template <class MD_>
template <typename IndicesOutputIterator>
IndicesOutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
get_incidences(Curve_index id,
               IndicesOutputIterator indices_out) const
{
  if(id <= 0 || std::size_t(id) > edges_incidences_.size())
    return indices_out;
  const Surface_patch_index_set& incidences = edges_incidences_[id - 1];
  return std::copy(incidences.begin(), incidences.end(), indices_out);
}

template <class MD_>
template <typename IndicesOutputIterator>
IndicesOutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
get_corner_incidences(Corner_index id,
                      IndicesOutputIterator indices_out) const
{
  typename Corners_incidences::const_iterator it = corners_incidences_.find(id);
  if(it == corners_incidences_.end())
    return indices_out;
  const Surface_patch_index_set& incidences = it->second;
  return std::copy(incidences.begin(), incidences.end(), indices_out);
}

template <class MD_>
template <typename IndicesOutputIterator>
IndicesOutputIterator
Mesh_domain_with_polyline_features_3<MD_>::
get_corner_incident_curves(Corner_index id,
                           IndicesOutputIterator indices_out) const
{
  if(id <= 0 || std::size_t(id) > corners_tmp_incidences_.size())
    return indices_out;
  const std::set<Curve_index>& incidences = corners_tmp_incidences_[id - 1];
  return std::copy(incidences.begin(), incidences.end(), indices_out);
}
/// @endcond

/// @cond CGAL_DOCUMENT_INTERNALS
namespace Mesh_3 {
namespace internal {

template <typename MDwPF_, bool curve_id_is_streamable>
template <typename Container2, typename Point>
void
Display_incidences_to_curves_aux<MDwPF_,curve_id_is_streamable>::
operator()(std::ostream& os, Point p, typename MDwPF_::Curve_index id,
           const Container2& corners_tmp_incidences_of_id) const
{
  os << "Corner #" << id << " (" << p
     << ") is incident to the following curves: {";
  for(typename MDwPF_::Curve_index curve_index :
                corners_tmp_incidences_of_id)
  {
    os << " " << curve_index;
  }
  os << " }\n";
}

template <class MDwPF_>
template <typename Container2, typename Point>
void
Display_incidences_to_curves_aux<MDwPF_,false>::
operator()(std::ostream& os, Point p, typename MDwPF_::Curve_index id,
           const Container2& corners_tmp_incidences_of_id) const
{
  os << "Corner #" << id << " (" << p
     << ") is incident to "
     << corners_tmp_incidences_of_id .size()
     << " curve(s).\n";
}

template <typename MDwPF_, bool patch_id_is_streamable>
template <typename Container, typename Point>
void
Display_incidences_to_patches_aux<MDwPF_,patch_id_is_streamable>::
operator()(std::ostream& os, Point p, typename MDwPF_::Curve_index id,
           const Container& corners_incidences_of_id) const
{
  os << "Corner #" << id << " (" << p
     << ") is incident to the following patches: {";
  for(typename MDwPF_::Surface_patch_index i :
                corners_incidences_of_id)
  {
    os << " " << i;
  }
  os << " }\n";
}

template <class MDwPF_>
template <typename Container, typename Point>
void
Display_incidences_to_patches_aux<MDwPF_,false>::
operator()(std::ostream& os, Point p, typename MDwPF_::Curve_index id,
           const Container& corners_incidences_id) const
{
  os << "Corner #" << id << " (" << p << ") is incident to "
     << corners_incidences_id.size()
     << " surface patch(es).\n";
}

} // end namespace Mesh_3::internal
} // end namespace Mesh_3
/// @endcond

/// @cond CGAL_DOCUMENT_INTERNALS
template <class MD_>
void
Mesh_domain_with_polyline_features_3<MD_>::
display_corner_incidences(std::ostream& os, Point_3 p, Corner_index id)
{
  typedef Mesh_domain_with_polyline_features_3<MD_> Mdwpf;
  typedef is_streamable<Surface_patch_index> i_s_spi;
  typedef is_streamable<Curve_index> i_s_csi;

  typedef Mesh_3::internal::Display_incidences_to_curves_aux<Mdwpf,i_s_csi::value> D_i_t_c;
  typedef Mesh_3::internal::Display_incidences_to_patches_aux<Mdwpf,i_s_spi::value> D_i_t_p;
  D_i_t_c()(os, p, id, corners_tmp_incidences_[id - 1]);
  D_i_t_p()(os, p, id, corners_incidences_[id]);
}
/// @endcond

template <class MD_>
void
Mesh_domain_with_polyline_features_3<MD_>::
compute_corners_incidences()
{
  for(typename Corners::iterator
        cit = corners_.begin(), end = corners_.end();
      cit != end; )
  {
    const Corner_index id = cit->second;
    const std::set<Curve_index>&
      corner_tmp_incidences = corners_tmp_incidences_[id - 1];

    if(corner_tmp_incidences.size() == 1 &&
       is_loop(*corner_tmp_incidences.begin()))
    {
      const Curve_index curve_id = *corner_tmp_incidences.begin();
      const Polyline& polyline = curve(curve_id);
      if(polyline.angle_at_first_point() == OBTUSE) {
        typename Corners::iterator to_erase = cit;
        ++cit;
        corners_.erase(to_erase);
        continue;
      }
    }

    Surface_patch_index_set& incidences = corners_incidences_[id];
    for(Curve_index curve_index : corner_tmp_incidences)
    {
      get_incidences(curve_index,
                     std::inserter(incidences,
                                   incidences.begin()));
    }
#if CGAL_MESH_3_PROTECTION_DEBUG & 1
    display_corner_incidences(std::cerr, cit->first, id);
#endif
    ++cit;
  }
}

/// @cond CGAL_DOCUMENT_INTERNALS
template <class MD_>
const typename Mesh_domain_with_polyline_features_3<MD_>::Surface_patch_index_set&
Mesh_domain_with_polyline_features_3<MD_>::
get_incidences(Curve_index id) const
{
  CGAL_precondition(id > 0 && std::size_t(id) <= edges_incidences_.size());
  return edges_incidences_[id - 1];
}

template <class MD_>
template <typename InputIterator>
typename Mesh_domain_with_polyline_features_3<MD_>::Curve_index
Mesh_domain_with_polyline_features_3<MD_>::
insert_edge(InputIterator first, InputIterator end)
{
  CGAL_assertion(std::distance(first,end) > 1);

  const Curve_index curve_index = current_curve_index_++;

  register_corner(*first, curve_index);
  if ( *first != *std::prev(end) )
  {
    register_corner(*std::prev(end), curve_index);
  }

  edges_.emplace_back();
  CGAL_assertion(curve_index == static_cast<Curve_index>(edges_.size()));
  Polyline& polyline = edges_.back();

  while ( first != end )
  {
    polyline.add_point(*first++);
  }

  edges_incidences_.emplace_back();
  return curve_index;
}
/// @endcond

template <class MD_>
CGAL::Sign
Mesh_domain_with_polyline_features_3<MD_>::
distance_sign(const Point_3& p, const Point_3& q,
              const Curve_index& index) const
{
  CGAL_assertion(index > 0 && std::size_t(index) <= edges_.size());
  const Polyline& poly = curve(index);
  CGAL_precondition( ! poly.is_loop() );

  if ( p == q )
    return CGAL::ZERO;
  else if ( poly.are_ordered_along(p,q) )
    return CGAL::POSITIVE;
  else
    return CGAL::NEGATIVE;
}


template <class MD_>
CGAL::Sign
Mesh_domain_with_polyline_features_3<MD_>::
distance_sign_along_loop(const Point_3& p,
                         const Point_3& q,
                         const Point_3& r,
                         const Curve_index& index) const
{
  CGAL_assertion(p != q);
  CGAL_assertion(p != r);
  CGAL_assertion(r != q);
  CGAL_assertion(index > 0 && std::size_t(index) <= edges_.size());
  const Polyline& poly = curve(index);
  CGAL_assertion(poly.is_loop());

  FT pq = poly.curve_segment_length(p,q,CGAL::POSITIVE);
  FT pr = poly.curve_segment_length(p,r,CGAL::POSITIVE);

  if ( pq <= pr ) { return CGAL::POSITIVE; }
  else { return CGAL::NEGATIVE; }
}

template <class MD_>
bool
Mesh_domain_with_polyline_features_3<MD_>::
is_loop(const Curve_index& index) const
{
  CGAL_assertion(index > 0 && std::size_t(index) <= edges_.size());
  return curve(index).is_loop();
}

template <class MD_>
bool
Mesh_domain_with_polyline_features_3<MD_>::
is_curve_segment_covered(const Curve_index& index,
                         CGAL::Orientation orientation,
                         const Point_3& c1, const Point_3& c2,
                         const FT sq_r1, const FT sq_r2) const
{
  CGAL_assertion(index > 0 && std::size_t(index) <= edges_.size());
  const Polyline& poly = curve(index);
  return poly.is_curve_segment_covered(orientation,
                                       c1, c2, sq_r1, sq_r2);
}

} //namespace CGAL

#endif // CGAL_MESH_DOMAIN_WITH_POLYLINE_FEATURES_3_H