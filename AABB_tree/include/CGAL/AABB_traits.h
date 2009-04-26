// Copyright (c) 2009 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: $
// $Id: $
//
//
// Author(s) : Stéphane Tayeb, Pierre Alliez
//
//******************************************************************************
// File Description :
//
//******************************************************************************

#ifndef AABB_TRAITS_H_
#define AABB_TRAITS_H_

#include <CGAL/Bbox_3.h>
#include <CGAL/AABB_intersections.h>

namespace CGAL {

/**
 * @class AABB_traits
 *
 *
 */
template<typename GeomTraits, typename AABB_primitive>
class AABB_traits
{
public:
  /// Ray query type
  typedef typename GeomTraits::Ray_3 Ray_3;
  /// Line query type
  typedef typename GeomTraits::Line_3 Line_3;
  /// Segment query type
  typedef typename GeomTraits::Segment_3 Segment_3;

  // TODO: delete once "inside..." disappears
  typedef typename GeomTraits::Triangle_3 Triangle_3;

  /// AABBTraits concept types
  typedef typename CGAL::Bbox_3 Bounding_box;

  typedef typename AABB_primitive Primitive;
  typedef typename AABB_primitive::Datum Datum;

  typedef typename GeomTraits::Sphere_3 Sphere;
  typedef typename GeomTraits::Point_3 Projection;
  typedef typename GeomTraits::Point_3 Intersection;
  typedef typename GeomTraits::Point_3 Projection_query;

  // types for search tree
  typedef typename GeomTraits::Cartesian_const_iterator_3 Cartesian_const_iterator_3;
  typedef typename GeomTraits::Construct_cartesian_const_iterator_3
      Construct_cartesian_const_iterator_3;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Iso_cuboid_3 Iso_cuboid_3;
  typedef typename GeomTraits::Sphere_3 Sphere_3;
  typedef typename GeomTraits::Construct_iso_cuboid_3 Construct_iso_cuboid_3;
  typedef typename GeomTraits::Construct_min_vertex_3 Construct_min_vertex_3;
  typedef typename GeomTraits::Construct_max_vertex_3 Construct_max_vertex_3;
  typedef typename GeomTraits::Construct_center_3 Construct_center_3;
  typedef typename GeomTraits::Compute_squared_radius_3 Compute_squared_radius_3;
  typedef typename GeomTraits::FT FT;

  /// Constructor
  AABB_traits() { };

  /// Non-virtual Destructor
  ~AABB_traits() { };

  /// Comparison functions
  static bool x_less_than(const Primitive& pr1, const Primitive& pr2)
  { return pr1.xref() < pr2.xref(); }
  static bool y_less_than(const Primitive& pr1, const Primitive& pr2)
  { return pr1.yref() < pr2.yref(); }
  static bool z_less_than(const Primitive& pr1, const Primitive& pr2)
  { return pr1.zref() < pr2.zref(); }

  /// UNDOCUMENTED FEATURE
  /// TODO: see what to do
  /**
   * @brief Sorts [first,beyond[
   * @param first iterator on first element
   * @param beyond iterator on beyond element
   * @param bbox the bounding box of [first,beyond[
   *
   * Sorts the range defined by [first,beyond[. Sort is achieved on bbox longuest
   * axis, using the comparison function <dim>_less_than (dim in {x,y,z})
   */
  template<typename PrimitiveIterator>
  void sort_primitives(PrimitiveIterator first,
                       PrimitiveIterator beyond,
                       const Bounding_box& bbox) const;

  /**
   * Computes the bounding box of a set of primitives
   * @param first an iterator on the first primitive
   * @param beyond an iterator on the past-the-end primitive
   * @return the bounding box of the primitives of the iterator range
   */
  template<typename ConstPrimitiveIterator>
  Bounding_box compute_bbox(ConstPrimitiveIterator first,
                            ConstPrimitiveIterator beyond) const;

  template<typename Query>
  bool do_intersect(const Query& q, const Bounding_box& bbox) const;

  template<typename Query>
  bool do_intersect(const Query& q, const Primitive& pr) const;

  template<typename Query>
  bool intersection(const Query& q,
                    const Primitive& pr,
                    Intersection& intersection) const;

  Sphere sphere(const Projection_query& center,
                const Projection& hint) const
  {
    return Sphere(center, GeomTraits().compute_squared_distance_3_object()
                                                          (center, hint));
  }

  bool intersection(const Sphere& sphere,
                    const Primitive& pr,
                    Projection& projection_return) const;

  bool is_smaller(const Sphere& a, const Sphere& b) const;

private:
  /// Private types
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;

private:
  /**
   * @brief Computes bounding box of one primitive
   * @param pr the primitive
   * @return the bounding box of the primitive \c pr
   */
  Bounding_box compute_bbox(const Primitive& pr) const
  {
    return pr.datum().bbox();
  }

  typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1,
                 CGAL_AXIS_Z = 2} Axis;

  Axis longest_axis(const Bounding_box& bbox) const;

  bool is_inside_triangle_3(Point_3& p, const Triangle_3& t) const;

private:
  // Disabled copy constructor & assignment operator
  typedef AABB_traits<GeomTraits, Primitive> Self;
  AABB_traits(const Self& src);
  Self& operator=(const Self& src);

};  // end class AABB_traits

template<typename GT, typename P>
template<typename PrimitiveIterator>
void
AABB_traits<GT,P>::sort_primitives(PrimitiveIterator first,
                                    PrimitiveIterator beyond,
                                    const Bounding_box& bbox) const
{
  PrimitiveIterator middle = first + (beyond - first)/2;
  switch(longest_axis(bbox))
  {
  case CGAL_AXIS_X: // sort along x
    std::nth_element(first, middle, beyond, x_less_than);
    break;
  case CGAL_AXIS_Y: // sort along y
    std::nth_element(first, middle, beyond, y_less_than);
    break;
  case CGAL_AXIS_Z: // sort along z
    std::nth_element(first, middle, beyond, z_less_than);
    break;
  default:
    CGAL_error();
  }
}

template<typename GT, typename P>
template<typename ConstPrimitiveIterator>
typename AABB_traits<GT,P>::Bounding_box
AABB_traits<GT,P>::compute_bbox(ConstPrimitiveIterator first,
                                 ConstPrimitiveIterator beyond) const
{
  Bounding_box bbox = compute_bbox(*first);
  for(++first; first != beyond; ++first)
  {
    bbox = bbox + compute_bbox(*first);
  }
  return bbox;
}


template<typename GT, typename P>
template<typename Query>
bool
AABB_traits<GT,P>::do_intersect(const Query& q,
                                 const Bounding_box& bbox) const
{
  // AABB tree package call
  // TODO: extend kernel
  return CGAL::do_intersect(q, bbox);
}


template<typename GT, typename P>
template<typename Query>
bool
AABB_traits<GT,P>::do_intersect(const Query& q,
                                 const P& pr) const
{
  return GT().do_intersect_3_object()(q, pr.datum());
}


template<typename GT, typename P>
template<typename Query>
bool
AABB_traits<GT,P>::intersection(const Query& q,
                                 const P& pr,
                                 Intersection& result) const
{
  // TODO: implement a real intersection construction method
  // do_intersect is needed here because we construct intersection between
  // pr.datum().supporting_plane() and q
  if ( ! do_intersect(q,pr) )
  {
    return false;
  }

  // AABB tree package call
  // TODO: extend kernel
  Datum datum = pr.datum();
  CGAL::Object intersection_obj = CGAL::intersection(datum, q);

  return CGAL::assign(result, intersection_obj);
}


// PA: CAREFUL: the ad-hoc code here must be removed.

template<typename GT, typename P>
bool
AABB_traits<GT,P>::intersection(const Sphere& sphere,
                                 const P& pr,
                                 Projection& projected) const
{
  typedef typename P::Datum Triangle_3;

  const Triangle_3 triangle = pr.datum();
  projected = triangle.supporting_plane().projection(sphere.center());

  // If point is projected outside sphere, return false
  if ( sphere.bounded_side(projected) == CGAL::ON_UNBOUNDED_SIDE )
  {
    return false;
  }

  if( is_inside_triangle_3(projected, triangle) )  // projected is modified
  {
    return true;
  }

  if(sphere.bounded_side(projected) == CGAL::ON_UNBOUNDED_SIDE)
  {
    return false;
  }

  return true;
}


template<typename GT, typename P>
bool
AABB_traits<GT,P>::is_smaller(const Sphere& a, const Sphere& b) const
{
  CGAL_precondition(a.center() == b.center());

  return ( GT().compute_squared_radius_3_object()(a)
           < GT().compute_squared_radius_3_object()(b) );
}



//-------------------------------------------------------
// Private methods
//-------------------------------------------------------
template<typename GT, typename P>
typename AABB_traits<GT,P>::Axis
AABB_traits<GT,P>::longest_axis(const Bounding_box& bbox) const
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

// PA: ad-hoc code to be removed
template<typename GT, typename P>
bool
AABB_traits<GT,P>::is_inside_triangle_3(Point_3& p,
                                         const Triangle_3& t) const
{
  typedef typename GT::Vector_3 Vector;
  typedef typename GT::Line_3 Line;

  Vector w = CGAL::cross_product(t.vertex(1) - t.vertex(0),
                                 t.vertex(2) - t.vertex(0));
  bool out = false;
  for(int i = 0; i < 3; ++i)
  {
    Vector v = CGAL::cross_product(t.vertex(i+1) - t.vertex(i),
                                   p - t.vertex(i));
    if(v*w < 0)
    { // p is outside, on the side of (i, i+1)
      out = true;
      if(   (p - t.vertex(i))*(t.vertex(i+1) - t.vertex(i)) >= 0
         && (p - t.vertex(i+1))*(t.vertex(i) - t.vertex(i+1)) >= 0 )
      {
        p = Line(t.vertex(i), t.vertex(i+1)).projection(p);
        return false;
      }
    }
  }

  if(out)
  {
    typename GT::Compare_distance_3 c = GT().compare_distance_3_object();
    if(c(p, t.vertex(1), t.vertex(0)) == CGAL::LARGER &&
       c(p, t.vertex(2), t.vertex(0)) == CGAL::LARGER)
    {
       p = t.vertex(0);
       return false;
    }
    if(c(p, t.vertex(0), t.vertex(1)) == CGAL::LARGER &&
       c(p, t.vertex(2), t.vertex(1)) == CGAL::LARGER)
    {
       p = t.vertex(1);
       return false;
    }
    p = t.vertex(2);
    return false;
  }

  return true;
}


}  // end namespace CGAL

#endif // AABB_TRAITS_H_
