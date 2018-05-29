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
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
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
#include <CGAL/AABB_traits.h>
#include <CGAL/is_streamable.h>
#include <CGAL/Real_timer.h>
#include <CGAL/property_map.h>

#include <vector>
#include <set>
#include <map>
#include <algorithm>

#include <boost/next_prior.hpp> // for boost::prior and boost::next
#include <boost/variant.hpp>
#include <boost/foreach.hpp>

namespace CGAL {

/// @cond DEVELOPERS
namespace internal {
namespace Mesh_3 {

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

  /// Add a point at the end of the polyline
  void add_point(const Point_3& p)
  {
    if( points_.empty() || p != end_point() ) {
      points_.push_back(p);
    }
  }

  /// Returns the starting point of the polyline
  const Point_3& start_point() const
  {
    CGAL_assertion( ! points_.empty() );
    return points_.front();
  }

  /// Returns the ending point of the polyline
  const Point_3& end_point() const
  {
    CGAL_assertion( ! points_.empty() );
    return points_.back();
  }

  /// Returns `true` if the polyline is not degenerated
  bool is_valid() const
  {
    return points_.size() > 1;
  }

  /// Returns `true` if polyline is a loop
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


  /// Returns the angle at the first point.
  /// \pre The polyline must be a loop.
  Angle angle_at_first_point() const {
    CGAL_precondition(is_loop());
    const Point_3& first = points_.front();
    const Point_3& next_p = points_[1];
    const Point_3& prev = points_[points_.size() - 2];
    return angle(prev, first, next_p);
  }

  /// Returns the length of the polyline
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

  /// Returns signed geodesic distance between \c p and \c q
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
        ?     curve_segment_length(p, q, CGAL::POSITIVE)
        : ( - curve_segment_length(p, q, CGAL::NEGATIVE) );
    }
  }


  /// Returns a point at geodesic distance \c distance from p along the
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

  /// Returns an iterator on the starting point of the segment of the
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
        if(compare_distance(p, seg, nearest_segment) == CGAL::SMALLER)
        {
          nearest_segment = seg;
          result = previous;
        }
        if(compare_distance(p, *it, nearest_segment) == CGAL::SMALLER)
        {
          nearest_vertex = it;
          nearest_is_a_segment = false;
          result = it;
        }
      }
      else {
        if(compare_distance(p, *it, *nearest_vertex) == CGAL::SMALLER)
        {
          nearest_vertex = it;
          result = it;
        }
        if(compare_distance(p, seg, *nearest_vertex) == CGAL::SMALLER)
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


template <typename Gt, typename MapIterator>
struct Mesh_domain_segment_of_curve_primitive{
  typedef typename std::iterator_traits<MapIterator>::value_type Map_value_type;
  typedef typename Map_value_type::first_type Curve_id;
  typedef typename Map_value_type::second_type Polyline;

  typedef std::pair<MapIterator,
                    typename Polyline::const_iterator> Id;

  typedef typename std::iterator_traits<
    typename Polyline::const_iterator>::value_type Point;

  typedef typename Gt::Segment_3 Datum;

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

} // end of namespace CGAL::internal::Mesh_3
} // end of namespace CGAL::internal
/// @endcond

/*!
\ingroup PkgMesh_3Domains

The class `Mesh_domain_with_polyline_features_3` is designed to allow the user
to add some 0- and 1-dimensional
features into any model of the `MeshDomain_3` concept.
The 1-dimensional features are described as polylines
whose endpoints are the added corners.

\tparam MeshDomain_3 is the type
of the domain which should be extended.
It has to be a model of the `MeshDomain_3` concept.

\cgalModels `MeshDomainWithFeatures_3`

\sa `MeshDomain_3`
\sa `MeshPolyline_3`
\sa `CGAL::Implicit_mesh_domain_3<Function,BGT>`
\sa `CGAL::Polyhedral_mesh_domain_3<Polyhedron,IGT,TriangleAccessor>`
\sa `CGAL::Labeled_image_mesh_domain_3<Image,BGT>`

*/
template < typename MeshDomain_3 >
class Mesh_domain_with_polyline_features_3
  : public MeshDomain_3
{
public:
/// \name Types
/// @{
  typedef typename MeshDomain_3::Index    Index;

  typedef typename MeshDomain_3::Surface_patch_index
                                  Surface_patch_index;

  typedef int                     Curve_index;
  typedef int                     Corner_index;
  typedef CGAL::Tag_true           Has_features;
  typedef typename MeshDomain_3::R::FT          FT;
/// @}

#ifndef DOXYGEN_RUNNING

#ifndef CGAL_NO_DEPRECATED_CODE
  typedef Curve_index Curve_segment_index;
#endif

  typedef typename MeshDomain_3::R         Gt;
  typedef Gt                       R;
  typedef typename MeshDomain_3::Point_3   Point_3;
#endif // DOXYGEN_RUNNING

#ifndef CGAL_CFG_NO_CPP0X_VARIADIC_TEMPLATES
/// \name Creation
/// @{

/*!
Constructor. Forwards the arguments to the constructor
of the base class.
*/
  template <typename ... T>
  Mesh_domain_with_polyline_features_3(const T& ...t)
    : MeshDomain_3(t...)
    , current_corner_index_(1)
    , current_curve_index_(1)
    , curves_aabb_tree_is_built(false) {}
/// @}
#else
  /// Constructors
  /// Call the base class constructor
  Mesh_domain_with_polyline_features_3()
    : MeshDomain_3()
    , current_corner_index_(1)
    , current_curve_index_(1)
    , curves_aabb_tree_is_built(false) {}

  template <typename T1>
  Mesh_domain_with_polyline_features_3(const T1& o1)
    : MeshDomain_3(o1)
    , current_corner_index_(1)
    , current_curve_index_(1)
    , curves_aabb_tree_is_built(false) {}

  template <typename T1, typename T2>
  Mesh_domain_with_polyline_features_3(const T1& o1, const T2& o2)
    : MeshDomain_3(o1, o2)
    , current_corner_index_(1)
    , current_curve_index_(1)
    , curves_aabb_tree_is_built(false) {}

  template <typename T1, typename T2, typename T3>
  Mesh_domain_with_polyline_features_3(const T1& o1, const T2& o2,
                                       const T3& o3)
    : MeshDomain_3(o1, o2, o3)
    , current_corner_index_(1)
    , current_curve_index_(1)
    , curves_aabb_tree_is_built(false) {}

  template <typename T1, typename T2, typename T3, typename T4>
  Mesh_domain_with_polyline_features_3(const T1& o1, const T2& o2,
                                       const T3& o3, const T4& o4)
    : MeshDomain_3(o1, o2, o3, o4)
    , current_corner_index_(1)
    , current_curve_index_(1)
    , curves_aabb_tree_is_built(false) {}
#endif

/// \name Operations
/// @{

  /// @cond DEVELOPERS
  /// @{

  /// Add a 0-dimensional feature in the domain.
  Corner_index add_corner(const Point_3& p);

  /// Overloads where the last parameter \c out is not `CGAL::Emptyset_iterator()`.
  template <typename InputIterator, typename IndicesOutputIterator>
  IndicesOutputIterator
  add_corners(InputIterator first, InputIterator end,
              IndicesOutputIterator out  /*= CGAL::Emptyset_iterator()*/);

  /*!
    Add 0-dimensional features in the domain. The value type of `InputIterator` must
    be `Point_3`.
  */
  template <typename InputIterator>
  void
  add_corners(InputIterator first, InputIterator end)
  { add_corners(first, end, CGAL::Emptyset_iterator()); }

  Corner_index register_corner(const Point_3& p, const Curve_index& index);
  Corner_index add_corner_with_context(const Point_3& p, const Surface_patch_index& index);

  /// Overloads where the last parameter \c out is not
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
   IncidentPatchesIndicesPMap incident_paches_indices_pmap,
   IndicesOutputIterator out /* = CGAL::Emptyset_iterator() */);
  
  template <typename InputIterator, typename IndicesOutputIterator>
  IndicesOutputIterator
  add_features_with_context(InputIterator first, InputIterator end,
                            IndicesOutputIterator out /*=
                                                        CGAL::Emptyset_iterator()*/);
  /// @}
  /// \endcond
  /*!
    Add 1-dimensional features in the domain. `InputIterator` value type must 
    be a model of the concept `MeshPolyline_3`.
  */
  template <typename InputIterator>
  void
  add_features(InputIterator first, InputIterator end)
  { add_features(first, end, CGAL::Emptyset_iterator()); }

  /// @cond DEVELOPERS
  /// Undocumented function, kept for backward-compatibility with existing
  /// code
  template <typename InputIterator>
  void
  add_features_with_context(InputIterator first, InputIterator end)
  { add_features_with_context(first, end, CGAL::Emptyset_iterator()); }
  /// @endcond

  /*!
    Add 1-dimensional features (curves) from the range `[first, end)` in the domain with their incidences
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
/// The following methods implement the requirement of the concept
/// `MeshDomainWithFeatures_3`.
/// @{

  /// Implements `MeshDomainWithFeatures_3::get_corners()`.
  /// OutputIterator value type is std::pair<Corner_index, Point_3>
  template <typename OutputIterator>
  OutputIterator get_corners(OutputIterator out) const;

  /// Implements `MeshDomainWithFeatures_3::get_curves()`.
  /// OutputIterator value type is CGAL::cpp11::tuple<Curve_index,
  /// std::pair<Point_3,Index>, std::pair<Point_3,Index> >
  template <typename OutputIterator>
  OutputIterator get_curves(OutputIterator out) const;

  /// Implements `MeshDomainWithFeatures_3::curve_segment_length()`.
  FT curve_segment_length(const Point_3& p, const Point_3 q,
                          const Curve_index& curve_index,
                          CGAL::Orientation orientation) const;

  /// Implements `MeshDomainWithFeatures_3::curve_length()`.
  FT curve_length(const Curve_index& curve_index) const;

  /// Implements `MeshDomainWithFeatures_3::construct_point_on_curve()`.
  Point_3
  construct_point_on_curve(const Point_3& starting_point,
                           const Curve_index& curve_index,
                           FT distance) const;
  /// Implements `MeshDomainWithFeatures_3::distance_sign_along_loop()`.
  CGAL::Sign distance_sign_along_loop(const Point_3& p,
                                      const Point_3& q,
                                      const Point_3& r,
                                      const Curve_index& index) const;

  /// Implements `MeshDomainWithFeatures_3::distance_sign()`.
  CGAL::Sign distance_sign(const Point_3& p, const Point_3& q,
                           const Curve_index& index) const;

  /// Implements `MeshDomainWithFeatures_3::is_loop()`.
  bool is_loop(const Curve_index& index) const;

  /// Implements `MeshDomainWithFeatures_3::is_curve_segment_covered()`.
  bool is_curve_segment_covered(const Curve_index& index,
                                CGAL::Orientation orientation,
                                const Point_3& c1, const Point_3& c2,
                                const FT sq_r1, const FT sq_r2) const;

  /// Returns an `Index` from a `Curve_index`
  Index index_from_curve_index(const Curve_index& index) const
  { return Index(index); }

  /// Returns a `Curve_index` from an `Index`
  Curve_index curve_index(const Index& index) const
  { return boost::get<Curve_index>(index); }

  /// Returns an `Index` from a `Corner_index`
  Index index_from_corner_index(const Corner_index& index) const
  { return Index(index); }

  /// Returns a `Corner_index` from an `Index`
  Corner_index corner_index(const Index& index) const
  { return boost::get<Corner_index>(index); }

  /// @cond DEVELOPERS
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
private:
  void compute_corners_incidences();

  /// Returns Index associated to p (p must be the coordinates of a corner
  /// point)
  Index point_corner_index(const Point_3& p) const;

private:
  typedef std::map<Point_3,Corner_index> Corners;

  typedef internal::Mesh_3::Polyline<Gt> Polyline;
  typedef std::map<Curve_index, Polyline> Edges;
  typedef std::map<Curve_index, Surface_patch_index_set > Edges_incidences;
  typedef std::map<Corner_index, std::set<Curve_index> > Corners_tmp_incidences;
  typedef std::map<Corner_index, Surface_patch_index_set > Corners_incidences;

  typedef internal::Mesh_3::Mesh_domain_segment_of_curve_primitive<
    Gt,
    typename Edges::const_iterator> Curves_primitives;

  typedef CGAL::AABB_traits<Gt,
                            Curves_primitives> AABB_curves_traits;

  Corners corners_;
  Corners_tmp_incidences corners_tmp_incidences_;
  Corner_index current_corner_index_;
  Corners_incidences corners_incidences_;

  Edges edges_;
  Curve_index current_curve_index_;
  Edges_incidences edges_incidences_;

public:
  /// @cond DEVELOPERS
  typedef CGAL::AABB_tree<AABB_curves_traits> Curves_AABB_tree;

private:
  mutable Curves_AABB_tree curves_aabb_tree_;
  mutable bool curves_aabb_tree_is_built;

public:
  const Corners_incidences& corners_incidences_map() const
  { return corners_incidences_; }

  const Curves_AABB_tree& curves_aabb_tree() const {
    if(!curves_aabb_tree_is_built) build_curves_aabb_tree();
    return curves_aabb_tree_;
  }
  Curve_index maximal_curve_index() const {
    if(edges_incidences_.empty()) return Curve_index();
    return boost::prior(edges_incidences_.end())->first;
  }

  void build_curves_aabb_tree() const {
#if CGAL_MESH_3_VERBOSE
    std::cerr << "Building curves AABB tree...";
    CGAL::Real_timer timer;
    timer.start();
#endif
    curves_aabb_tree_.clear();
    for(typename Edges::const_iterator
          edges_it = edges_.begin(),
          edges_end = edges_.end();
        edges_it != edges_end; ++edges_it)
    {
      const Polyline& polyline = edges_it->second;
      for(typename Polyline::const_iterator
            pit = polyline.points_.begin(),
            end = polyline.points_.end() - 1;
          pit != end; ++pit)
      {
        curves_aabb_tree_.insert(std::make_pair(edges_it, pit));
      }
    }
    curves_aabb_tree_.build();
    curves_aabb_tree_is_built = true;
#if CGAL_MESH_3_VERBOSE
    timer.stop();
    std::cerr << " done (" << timer.time() * 1000 << " ms)" << std::endl;
#endif
  } // end build_curves_aabb_tree()
  /// @endcond
private:
  // Disabled copy constructor & assignment operator
  typedef Mesh_domain_with_polyline_features_3 Self;
  Mesh_domain_with_polyline_features_3(const Self& src);
  Self& operator=(const Self& src);

};  // end class Mesh_domain_with_polyline_features_3



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
  for ( typename Edges::const_iterator
       eit = edges_.begin(), end = edges_.end() ; eit != end ; ++eit )
  {
    CGAL_assertion( eit->second.is_valid() );

    const Point_3& p = eit->second.start_point();
    const Point_3& q = eit->second.end_point();

    Index p_index, q_index;
    if ( ! eit->second.is_loop() )
    {
      p_index = point_corner_index(p);
      q_index = point_corner_index(q);
    }
    else
    {
      p_index = index_from_curve_index(eit->first);
      q_index = p_index;
    }

    *out++ = CGAL::cpp11::make_tuple(eit->first,
                                     std::make_pair(p,p_index),
                                     std::make_pair(q,q_index));
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
  // Get corresponding polyline
  typename Edges::const_iterator eit = edges_.find(curve_index);
  CGAL_assertion(eit != edges_.end());

  return eit->second.curve_segment_length(p, q, orientation);
}


template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::FT
Mesh_domain_with_polyline_features_3<MD_>::
curve_length(const Curve_index& curve_index) const
{
  // Get corresponding polyline
  typename Edges::const_iterator eit = edges_.find(curve_index);
  CGAL_assertion(eit != edges_.end());

  return eit->second.length();
}


template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Point_3
Mesh_domain_with_polyline_features_3<MD_>::
construct_point_on_curve(const Point_3& starting_point,
                         const Curve_index& curve_index,
                         FT distance) const
{
  // Get corresponding polyline
  typename Edges::const_iterator eit = edges_.find(curve_index);
  CGAL_assertion(eit != edges_.end());

  // Return point at geodesic_distance distance from starting_point
  return eit->second.point_at(starting_point,distance);
}


/// @cond DEVELOPERS
template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::Corner_index
Mesh_domain_with_polyline_features_3<MD_>::
add_corner(const Point_3& p)
{
  typename Corners::iterator cit = corners_.lower_bound(p);

  // If the corner already exists, return its assigned Corner_index...
  if(cit != corners_.end() && !(corners_.key_comp()(p, cit->first)))
    return cit->second;

  // ... otherwise, insert it!
  const Corner_index index = current_corner_index_++;
  corners_.insert(cit, std::make_pair(p, index));

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
  // 'add_corner' will itself seek if 'p' is already a corner, and, in that case,
  // return the Corner_index that has been assigned to this position.
  Corner_index index = add_corner(p);
  corners_tmp_incidences_[index].insert(curve_index);

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
  // Insert one edge for each element
  while ( first != end )
  {
    *indices_out++ = insert_edge(first->begin(), first->end());
    ++first;
  }
  compute_corners_incidences();
  return indices_out;
}

/// @cond DEVELOPERS
namespace details {

template <typename PolylineWithContext>
struct Get_content_from_polyline_with_context {
  typedef Get_content_from_polyline_with_context Self;
  typedef const PolylineWithContext& key_type;
  typedef const typename PolylineWithContext::Bare_polyline& value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  friend value_type get(const Self, key_type polyline) {
    return polyline.polyline_content;
  }
}; // end Get_content_from_polyline_with_context<PolylineWithContext>

template <typename PolylineWithContext>
struct Get_patches_id_from_polyline_with_context {
  typedef Get_patches_id_from_polyline_with_context Self;
  typedef const PolylineWithContext& key_type;
  typedef const typename PolylineWithContext::Context::Patches_ids& value_type;
  typedef value_type reference;
  typedef boost::readable_property_map_tag category;
  friend value_type get(const Self, key_type polyline) {
    return polyline.context.adjacent_patches_ids;
  }
}; // end Get_patches_id_from_polyline_with_context<PolylineWithContext>

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
  // Insert one edge for each element
  for( ; first != end ; ++first )
  {
    const typename boost::property_traits<PolylinePMap>::reference
      polyline = get(polyline_pmap, *first);
    const typename boost::property_traits<IncidentPatchesIndicesPMap>::reference
      patches_ids = get(inc_patches_ind_pmap, *first);

    Curve_index curve_id = insert_edge(polyline.begin(), polyline.end());
    edges_incidences_[curve_id].insert(patches_ids.begin(), patches_ids.end());
    *indices_out++ = curve_id;
  }

  compute_corners_incidences();
  return indices_out;
}

/// @cond DEVELOPERS
template <class MD_>
typename Mesh_domain_with_polyline_features_3<MD_>::FT
Mesh_domain_with_polyline_features_3<MD_>::
signed_geodesic_distance(const Point_3& p, const Point_3& q,
                         const Curve_index& curve_index) const
{
  // Get corresponding polyline
  typename Edges::const_iterator eit = edges_.find(curve_index);
  CGAL_assertion(eit != edges_.end());

  // Compute geodesic_distance
  return eit->second.signed_geodesic_distance(p,q);
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
  BOOST_FOREACH(typename IncidenceMap::value_type& pair,
                incidence_map)
  {
    Surface_patch_index_set& patch_index_set = pair.second;
    Surface_patch_index_set new_index_set;
    for(typename Surface_patch_index_set::const_iterator
        it = patch_index_set.begin(), end = patch_index_set.end();
        it != end; ++it)
    {
      CGAL_assertion(std::size_t(*it) < map.size());
      new_index_set.insert(map[*it]);
    }
    pair.second = new_index_set;
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
  typename Edges_incidences::const_iterator it = edges_incidences_.find(id);

  if(it == edges_incidences_.end())
    return indices_out;

  const Surface_patch_index_set& incidences = it->second;

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
  typename Corners_tmp_incidences::const_iterator it = corners_tmp_incidences_.find(id);
  if(it == corners_tmp_incidences_.end())
    return indices_out;

  const std::set<Curve_index>& incidences = it->second;
  return std::copy(incidences.begin(), incidences.end(), indices_out);
}
/// @endcond

/// @cond DEVELOPERS
namespace internal { namespace Mesh_3 {

template <typename MDwPF_, bool curve_id_is_streamable>
// here 'curve_id_is_streamable' is true
template <typename Container2, typename Point>
void
Display_incidences_to_curves_aux<MDwPF_,curve_id_is_streamable>::
operator()(std::ostream& os, Point p, typename MDwPF_::Curve_index id,
           const Container2& corners_tmp_incidences_of_id) const
{
  os << "Corner #" << id << " (" << p
     << ") is incident to the following curves: {";
  BOOST_FOREACH(typename MDwPF_::Curve_index curve_index,
                corners_tmp_incidences_of_id)
  {
    os << " " << curve_index;
  }
  os << " }\n";
}

template <class MDwPF_>
// here 'curve_id_is_streamable' is false
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
// here 'patch_id_is_streamable' is true
template <typename Container, typename Point>
void
Display_incidences_to_patches_aux<MDwPF_,patch_id_is_streamable>::
operator()(std::ostream& os, Point p, typename MDwPF_::Curve_index id,
           const Container& corners_incidences_of_id) const
{
  os << "Corner #" << id << " (" << p
     << ") is incident to the following patches: {";
  BOOST_FOREACH(typename MDwPF_::Surface_patch_index i,
                corners_incidences_of_id)
  {
    os << " " << i;
  }
  os << " }\n";
}

template <class MDwPF_>
// here 'patch_id_is_streamable' is false
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

}} // end namespaces internal::Mesh_3:: and internal::
/// @endcond

/// @cond DEVELOPERS
template <class MD_>
void
Mesh_domain_with_polyline_features_3<MD_>::
display_corner_incidences(std::ostream& os, Point_3 p, Corner_index id)
{
  typedef Mesh_domain_with_polyline_features_3<MD_> Mdwpf;
  typedef is_streamable<Surface_patch_index> i_s_spi;
  typedef is_streamable<Curve_index> i_s_csi;

  using namespace internal::Mesh_3;
  typedef Display_incidences_to_curves_aux<Mdwpf,i_s_csi::value> D_i_t_c;
  typedef Display_incidences_to_patches_aux<Mdwpf,i_s_spi::value> D_i_t_p;
  D_i_t_c()(os, p, id, corners_tmp_incidences_[id]);
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
      cit != end; /* the loop variable is incremented in the  body */)
  {
    const Corner_index id = cit->second;

    const typename Corners_tmp_incidences::mapped_type&
      corner_tmp_incidences = corners_tmp_incidences_[id];

    // If the corner is incident to only one curve, and that curve is a
    // loop, then remove the corner from the set, only if the angle is not
    // acute. If the angle is acute, the corner must remain as a corner,
    // to deal correctly with the angle.
    if(corner_tmp_incidences.size() == 1 &&
       is_loop(*corner_tmp_incidences.begin()))
    {
      const Curve_index curve_id = *corner_tmp_incidences.begin();
      const Polyline& polyline = edges_[curve_id];
      if(polyline.angle_at_first_point() == OBTUSE) {
        typename Corners::iterator to_erase = cit;
        ++cit;
        corners_.erase(to_erase);
        continue;
      }
    }

    Surface_patch_index_set& incidences = corners_incidences_[id];

    BOOST_FOREACH(Curve_index curve_index, corner_tmp_incidences)
    {
      get_incidences(curve_index,
                     std::inserter(incidences,
                                   incidences.begin()));
    }
#ifdef CGAL_MESH_3_PROTECTION_DEBUG
    display_corner_incidences(std::cerr, cit->first, id);
#endif // CGAL_MESH_3_PROTECTION_DEBUG

    // increment the loop variable
    ++cit;
  }
}

/// @cond DEVELOPERS
template <class MD_>
const typename Mesh_domain_with_polyline_features_3<MD_>::Surface_patch_index_set&
Mesh_domain_with_polyline_features_3<MD_>::
get_incidences(Curve_index id) const
{
  typename Edges_incidences::const_iterator it = edges_incidences_.find(id);
  CGAL_assertion(it != edges_incidences_.end());

  return it->second;
}

template <class MD_>
template <typename InputIterator>
typename Mesh_domain_with_polyline_features_3<MD_>::Curve_index
Mesh_domain_with_polyline_features_3<MD_>::
insert_edge(InputIterator first, InputIterator end)
{
  CGAL_assertion(std::distance(first,end) > 1);

  const Curve_index curve_index = current_curve_index_++;

  // Fill corners
  //
  // For a loop, the "first" point of the loop is registered as a
  // corner. If at the end, during the call to
  // 'compute_corners_incidences()', that corner is incident only to a
  // loop, then it will be removed from the set of corners.
  register_corner(*first, curve_index);
  if ( *first != *boost::prior(end) )
  {
    register_corner(*boost::prior(end), curve_index);
  }

  // Create a new polyline
  std::pair<typename Edges::iterator,bool> insertion =
    edges_.insert(std::make_pair(curve_index,Polyline()));

  // Fill polyline with data
  while ( first != end )
  {
    insertion.first->second.add_point(*first++);
  }
  return curve_index;
}
/// @endcond

template <class MD_>
CGAL::Sign
Mesh_domain_with_polyline_features_3<MD_>::
distance_sign(const Point_3& p, const Point_3& q,
              const Curve_index& index) const
{
  typename Edges::const_iterator eit = edges_.find(index);
  CGAL_assertion(eit != edges_.end());
  CGAL_precondition( ! eit->second.is_loop() );

  if ( p == q )
    return CGAL::ZERO;
  else if ( eit->second.are_ordered_along(p,q) )
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

  // Find edge
  typename Edges::const_iterator eit = edges_.find(index);
  CGAL_assertion(eit != edges_.end());
  CGAL_assertion(eit->second.is_loop());

  FT pq = eit->second.curve_segment_length(p,q,CGAL::POSITIVE);
  FT pr = eit->second.curve_segment_length(p,r,CGAL::POSITIVE);

  // Compare pq and pr
  if ( pq <= pr ) { return CGAL::POSITIVE; }
  else { return CGAL::NEGATIVE; }
}

template <class MD_>
bool
Mesh_domain_with_polyline_features_3<MD_>::
is_loop(const Curve_index& index) const
{
  // Find edge
  typename Edges::const_iterator eit = edges_.find(index);
  CGAL_assertion(eit != edges_.end());

  return eit->second.is_loop();
}

template <class MD_>
bool
Mesh_domain_with_polyline_features_3<MD_>::
is_curve_segment_covered(const Curve_index& index,
                         CGAL::Orientation orientation,
                         const Point_3& c1, const Point_3& c2,
                         const FT sq_r1, const FT sq_r2) const
{
  typename Edges::const_iterator eit = edges_.find(index);
  CGAL_assertion(eit != edges_.end());

  return eit->second.is_curve_segment_covered(orientation,
                                              c1, c2, sq_r1, sq_r2);
}



} //namespace CGAL


#endif // CGAL_MESH_DOMAIN_WITH_POLYLINE_FEATURES_3_H
