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
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez, Camille Wormser
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef CGAL_AABB_TRAITS_H_
#define CGAL_AABB_TRAITS_H_

#include <CGAL/Bbox_3.h>
#include <CGAL/AABB_intersections.h>

#include <boost/optional.hpp>

/// \file AABB_traits.h

namespace CGAL {

/// \addtogroup PkgAABB_tree
/// @{

/// The class AABB_traits is a model of the concept \ref
/// AABBTraits. This traits class handles any type of 3D geometric
/// primitives provided that the proper intersection tests and
/// constructions are implemented. It handles points, rays, lines and
/// segments as query types for intersection detection and
/// computations, and it handles points as query type for distance
/// queries. 
/// \tparam GeomTraits must  be a model of the concept \ref AABBGeomTraits,
/// snd provide the geometric types as well as the intersection tests and computations.
/// \tparam Primitive must be a model of the concept \ref AABBPrimitive and provide the
/// type of primitives stored in the AABB_tree.
///
/// \sa `AABBTraits`
/// \sa `AABB_tree`
/// \sa `AABBPrimitive`
template<typename GeomTraits, typename AABB_primitive>
class AABB_traits
{
  typedef typename CGAL::Object Object;
public:
  typedef AABB_traits<GeomTraits, AABB_primitive> AT;
  // AABBTraits concept types
  typedef typename GeomTraits::FT FT;
  typedef AABB_primitive Primitive;

  typedef typename std::pair<Object,typename Primitive::Id> Object_and_primitive_id;
  typedef typename std::pair<typename GeomTraits::Point_3, typename Primitive::Id> Point_and_primitive_id;

  // types for search tree
  /// \name Types
  /// @{

  /// Point query type.
  typedef typename GeomTraits::Point_3 Point_3;

  /// additionnal types for the search tree, required by the RangeSearchTraits concept
  /// \bug This is not documented for now in the AABBTraits concept.
  typedef typename GeomTraits::Iso_cuboid_3 Iso_cuboid_3;

  /// 
  typedef typename CGAL::Bbox_3 Bounding_box;

  /// @}

  typedef typename GeomTraits::Sphere_3 Sphere_3;
  typedef typename GeomTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_3; 
  typedef typename GeomTraits::Construct_cartesian_const_iterator_3 Construct_cartesian_const_iterator_3;
  typedef typename GeomTraits::Construct_center_3 Construct_center_3;
  typedef typename GeomTraits::Compute_squared_radius_3 Compute_squared_radius_3;
  typedef typename GeomTraits::Construct_min_vertex_3 Construct_min_vertex_3;
  typedef typename GeomTraits::Construct_max_vertex_3 Construct_max_vertex_3;  
  typedef typename GeomTraits::Construct_iso_cuboid_3 Construct_iso_cuboid_3;
  

  /// Default constructor.
  AABB_traits() { };


  typedef typename GeomTraits::Compute_squared_distance_3 Squared_distance;
  Squared_distance squared_distance_object() const { return GeomTraits().compute_squared_distance_3_object(); }

  /**
   * @internal
   * @brief Sorts [first,beyond[
   * @param first iterator on first element
   * @param beyond iterator on beyond element
   * @param bbox the bounding box of [first,beyond[
   *
   * Sorts the range defined by [first,beyond[. Sort is achieved on bbox longuest
   * axis, using the comparison function `<dim>_less_than` (dim in {x,y,z})
   */
class Sort_primitives
{
public:
template<typename PrimitiveIterator>
void operator()(PrimitiveIterator first,
                PrimitiveIterator beyond,
                const typename AT::Bounding_box& bbox) const
  {
    PrimitiveIterator middle = first + (beyond - first)/2;
    switch(longest_axis(bbox))
    {
    case AT::CGAL_AXIS_X: // sort along x
      std::nth_element(first, middle, beyond, less_x);
      break;
    case AT::CGAL_AXIS_Y: // sort along y
      std::nth_element(first, middle, beyond, less_y);
      break;
    case AT::CGAL_AXIS_Z: // sort along z
      std::nth_element(first, middle, beyond, less_z);
      break;
    default:
      CGAL_error();
    }
  }
};

Sort_primitives sort_primitives_object() {return Sort_primitives();}


  /*
   * Computes the bounding box of a set of primitives
   * @param first an iterator on the first primitive
   * @param beyond an iterator on the past-the-end primitive
   * @return the bounding box of the primitives of the iterator range
   */
   class Compute_bbox {
public:
template<typename ConstPrimitiveIterator>
typename AT::Bounding_box operator()(ConstPrimitiveIterator first,
                                     ConstPrimitiveIterator beyond) const
  {
    typename AT::Bounding_box bbox = compute_bbox(*first);
    for(++first; first != beyond; ++first)
    {
      bbox = bbox + compute_bbox(*first);
    }
    return bbox;
  }
};

Compute_bbox compute_bbox_object() {return Compute_bbox();}


class Do_intersect {
public:
  template<typename Query>
  bool operator()(const Query& q, const Bounding_box& bbox) const
  {
    return CGAL::do_intersect(q, bbox);
  }

  template<typename Query>
  bool operator()(const Query& q, const Primitive& pr) const
  {
    return GeomTraits().do_intersect_3_object()(q, pr.datum());
  }
};

Do_intersect do_intersect_object() {return Do_intersect();}

class Intersection {
public:
template<typename Query>
boost::optional<typename AT::Object_and_primitive_id>
operator()(const Query& query, const typename AT::Primitive& primitive) const
{
  typedef boost::optional<Object_and_primitive_id> Intersection;

  CGAL::Object object = GeomTraits().intersect_3_object()(primitive.datum(),query);
  if ( object.empty() )
    return Intersection();
  else
    return Intersection(Object_and_primitive_id(object,primitive.id()));
}
};

Intersection intersection_object() {return Intersection();}


  // This should go down to the GeomTraits, i.e. the kernel
  class Closest_point {
      typedef typename AT::Point_3 Point;
      typedef typename AT::Primitive Primitive;
  public:
      Point operator()(const Point& p, const Primitive& pr, const Point& bound) const
      {
          return CGAL::nearest_point_3(p, pr.datum(), bound);
      }
  };

  // This should go down to the GeomTraits, i.e. the kernel
  // and the internal implementation should change its name from
  // do_intersect to something like does_contain (this is what we compute,
  // this is not the same do_intersect as the spherical kernel)
  class Compare_distance {
      typedef typename AT::Point_3 Point;
      typedef typename AT::FT FT;
      typedef typename AT::Primitive Primitive;
  public:
      template <class Solid>
      CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const Point& bound) const
      {
          return GeomTraits().do_intersect_3_object()
          (GeomTraits().construct_sphere_3_object()
          (p, GeomTraits().compute_squared_distance_3_object()(p, bound)), pr)?
          CGAL::SMALLER : CGAL::LARGER;
      }

      template <class Solid>
      CGAL::Comparison_result operator()(const Point& p, const Solid& pr, const FT& sq_distance) const
      {
        return GeomTraits().do_intersect_3_object()
          (GeomTraits().construct_sphere_3_object()(p, sq_distance),
           pr) ?
          CGAL::SMALLER : 
          CGAL::LARGER;
      }
  };

  Closest_point closest_point_object() {return Closest_point();}
  Compare_distance compare_distance_object() {return Compare_distance();}


private:
  /**
   * @brief Computes bounding box of one primitive
   * @param pr the primitive
   * @return the bounding box of the primitive \c pr
   */
  static Bounding_box compute_bbox(const Primitive& pr)
  {
    return pr.datum().bbox();
  }

  typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1,
                 CGAL_AXIS_Z = 2} Axis;

  static Axis longest_axis(const Bounding_box& bbox);

  /// Comparison functions
  static bool less_x(const Primitive& pr1, const Primitive& pr2)
  { return pr1.reference_point().x() < pr2.reference_point().x(); }
  static bool less_y(const Primitive& pr1, const Primitive& pr2)
  { return pr1.reference_point().y() < pr2.reference_point().y(); }
  static bool less_z(const Primitive& pr1, const Primitive& pr2)
  { return pr1.reference_point().z() < pr2.reference_point().z(); }

};  // end class AABB_traits


//-------------------------------------------------------
// Private methods
//-------------------------------------------------------
template<typename GT, typename P>
typename AABB_traits<GT,P>::Axis
AABB_traits<GT,P>::longest_axis(const Bounding_box& bbox)
{
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();
  const double dz = bbox.zmax() - bbox.zmin();

  if(dx>=dy)
  {
    if(dx>=dz)
    {
      return CGAL_AXIS_X;
    }
    else // dz>dx and dx>=dy
    {
      return CGAL_AXIS_Z;
    }
  }
  else // dy>dx
  {
    if(dy>=dz)
    {
      return CGAL_AXIS_Y;
    }
    else  // dz>dy and dy>dx
    {
      return CGAL_AXIS_Z;
    }
  }
}

/// @}

}  // end namespace CGAL

#endif // CGAL_AABB_TRAITS_H_
