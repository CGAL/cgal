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

#include <CGAL/AABB_tree/Ray_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Bbox_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Segment_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Plane_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Triangle_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Line_3_Bbox_3_do_intersect.h>
#include <CGAL/AABB_tree/Sphere_3_Bbox_do_intersect.h>
#include <CGAL/AABB_tree/Triangle_3_segment_3_intersection.h>
#include <CGAL/AABB_tree/Triangle_3_ray_3_intersection.h>


namespace CGAL {

/**
 * @class AABB_traits
 *
 *
 */
template<typename GeomTraits, typename Primitive>
class AABB_traits
{
public:
  /// Ray query type
  typedef typename GeomTraits::Ray_3 Ray_3;
  /// Line query type
  typedef typename GeomTraits::Line_3 Line_3;
  /// Segment query type
  typedef typename GeomTraits::Segment_3 Segment_3;

  /// AABBTraits concept types
  typedef typename CGAL::Bbox_3 Bounding_box;

  typedef typename Primitive::Data Data;

  typedef typename GeomTraits::Sphere_3 Sphere;
  typedef typename GeomTraits::Point_3 Projection;
  typedef typename GeomTraits::Point_3 Intersection;
  typedef typename GeomTraits::Point_3 Projection_query;


  /// Constructor
  AABB_traits() { };

  /// Non-virtual Destructor
  ~AABB_traits() { };

  /// Comparison functions
  static bool x_less_than(const Primitive& pr1, const Primitive& pr2);
  static bool y_less_than(const Primitive& pr1, const Primitive& pr2);
  static bool z_less_than(const Primitive& pr1, const Primitive& pr2);

  /// UNDOCUMENTED FEATURE
  /// TODO: see what to do
  /**
   * @brief Sorts [first,last[
   * @param first iterator on first element
   * @param last iterator on last element
   * @param bbox the bounding box of [first,last[
   *
   * Sorts the range defined by [first,last[. Sort is achieved on bbox longuest
   * axis, using the comparison function <dim>_less_than (dim in {x,y,z})
   */
  template<typename PrimitiveIterator>
  void sort_primitives(PrimitiveIterator first,
                       PrimitiveIterator last,
                       const Bounding_box& bbox) const;

  /**
   * Computes the bounding box of a set of primitives
   * @param first an iterator on the first primitive
   * @param last an iterator on the last primitive
   * @return the bounding box of the primitives of the iterator range
   */
  template<typename ConstPrimitiveIterator>
  Bounding_box compute_bbox(ConstPrimitiveIterator first,
                            ConstPrimitiveIterator last) const;

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
  typedef typename TrianglePrimitive::Triangle_3 Triangle_3;

private:
  /**
   * @brief Computes bounding box of one primitive
   * @param pr the primitive
   * @return the bounding box of the primitive \c pr
   */
  Bounding_box compute_bbox(const Primitive& pr) const
  {
    return pr.data().bbox();
  }

  typedef enum { CGAL_AXIS_X = 0,
                 CGAL_AXIS_Y = 1,
                 CGAL_AXIS_Z = 2} Axis;

  Axis longest_axis(const Bounding_box& bbox) const;

  bool is_inside_triangle_3(Point_3& p, const Triangle_3& t) const;

private:
  // Disabled copy constructor & assignment operator
  typedef AABB_traits<GeomTraits, TrianglePrimitive> Self;
  AABB_traits(const Self& src);
  Self& operator=(const Self& src);

};  // end class AABB_traits



template<typename GT, typename TP>
bool
AABB_traits<GT,TP>::x_less_than(const Primitive& pr1,
                                const Primitive& pr2)
{
  const FT& ax1 = pr1.data().vertex(0).x();
  const FT& bx1 = pr1.data().vertex(1).x();
  const FT& cx1 = pr1.data().vertex(2).x();

  const FT& ax2 = pr2.data().vertex(0).x();
  const FT& bx2 = pr2.data().vertex(1).x();
  const FT& cx2 = pr2.data().vertex(2).x();

  return (ax1+bx1+cx1) < (ax2+bx2+cx2);
}

template<typename GT, typename TP>
bool
AABB_traits<GT,TP>::y_less_than(const Primitive& pr1,
                                const Primitive& pr2)
{
  const FT& ay1 = pr1.data().vertex(0).y();
  const FT& by1 = pr1.data().vertex(1).y();
  const FT& cy1 = pr1.data().vertex(2).y();

  const FT& ay2 = pr2.data().vertex(0).y();
  const FT& by2 = pr2.data().vertex(1).y();
  const FT& cy2 = pr2.data().vertex(2).y();

  return (ay1+by1+cy1) < (ay2+by2+cy2);
}

template<typename GT, typename TP>
bool
AABB_traits<GT,TP>::z_less_than(const Primitive& pr1,
                                const Primitive& pr2)
{
  const FT& az1 = pr1.data().vertex(0).z();
  const FT& bz1 = pr1.data().vertex(1).z();
  const FT& cz1 = pr1.data().vertex(2).z();

  const FT& az2 = pr2.data().vertex(0).z();
  const FT& bz2 = pr2.data().vertex(1).z();
  const FT& cz2 = pr2.data().vertex(2).z();

  return (az1+bz1+cz1) < (az2+bz2+cz2);
}


template<typename GT, typename TP>
template<typename PrimitiveIterator>
void
AABB_traits<GT,TP>::sort_primitives(PrimitiveIterator first,
                                    PrimitiveIterator last,
                                    const Bounding_box& bbox) const
{
  PrimitiveIterator middle = first + (last - first)/2;
  switch(longest_axis(bbox))
  {
  case CGAL_AXIS_X: // sort along x
    std::nth_element(first, middle, last, x_less_than);
    break;
  case CGAL_AXIS_Y: // sort along y
    std::nth_element(first, middle, last, y_less_than);
    break;
  case CGAL_AXIS_Z: // sort along z
    std::nth_element(first, middle, last, z_less_than);
    break;
  default:
    CGAL_error();
  }
}

template<typename GT, typename TP>
template<typename ConstPrimitiveIterator>
typename AABB_traits<GT,TP>::Bounding_box
AABB_traits<GT,TP>::compute_bbox(ConstPrimitiveIterator first,
                                 ConstPrimitiveIterator last) const
{
  Bounding_box bbox = compute_bbox(*first);
  for(++first; first != last; ++first)
  {
    bbox = bbox + compute_bbox(*first);
  }
  return bbox;
}


template<typename GT, typename TP>
template<typename Query>
bool
AABB_traits<GT,TP>::do_intersect(const Query& q,
                                 const Bounding_box& bbox) const
{
  // AABB tree package call
  // TODO: extend kernel
  return CGAL::do_intersect(q, bbox);
}


template<typename GT, typename TP>
template<typename Query>
bool
AABB_traits<GT,TP>::do_intersect(const Query& q,
                                 const Primitive& pr) const
{
  return GT().do_intersect_3_object()(q, pr.data());
}


template<typename GT, typename TP>
template<typename Query>
bool
AABB_traits<GT,TP>::intersection(const Query& q,
                                 const Primitive& pr,
                                 Intersection& intersection) const
{
  // TODO: implement a real intersection construction method
  // do_intersect is needed here because we construct intersection between
  // pr.data().supporting_plane() and q
  if ( ! do_intersect(q,pr) )
  {
    return false;
  }

  // AABB tree package call
  // TODO: extend kernel
  Object intersection_obj = CGAL::intersection(pr.data(), q);

  return CGAL::assign(intersection, intersection_obj);
}



template<typename GT, typename TP>
bool
AABB_traits<GT,TP>::intersection(const Sphere& sphere,
                                 const Primitive& pr,
                                 Projection& projected) const
{
  typedef typename TP::Data Triangle_3;

  const Triangle_3 triangle = pr.data();
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


template<typename GT, typename TP>
bool
AABB_traits<GT,TP>::is_smaller(const Sphere& a, const Sphere& b) const
{
  CGAL_precondition(a.center() == b.center());

  return ( GT().compute_squared_radius_3_object()(a)
           < GT().compute_squared_radius_3_object()(b) );
}



//-------------------------------------------------------
// Private methods
//-------------------------------------------------------
template<typename GT, typename TP>
typename AABB_traits<GT,TP>::Axis
AABB_traits<GT,TP>::longest_axis(const Bounding_box& bbox) const
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


template<typename GT, typename TP>
bool
AABB_traits<GT,TP>::is_inside_triangle_3(Point_3& p,
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
