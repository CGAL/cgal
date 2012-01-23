// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// $URL$
// $Id$
// 
//
// Author(s)     : Ron Wein   <wein@post.tau.ac.il>
//                 (based on old version by Eyal Flato)
//                 Efi Fogel  <efif@post.tau.ac.il>

#ifndef CGAL_ARR_SIMPLE_POINT_LOCATION_H
#define CGAL_ARR_SIMPLE_POINT_LOCATION_H

/*! \file
 * Definition of the Arr_simple_point_location<Arrangement> template.
 */

#include <CGAL/Arr_point_location/Arr_point_location.h>
#include <CGAL/Arrangement_2/Arr_traits_adaptor_2.h>
#include <CGAL/Object.h>

#include <boost/variant.hpp>
#include <boost/optional.hpp>

namespace CGAL {

/*! \class
 * A class that answers point-location and vertical ray-shooting queries
 * on a 2D arrangement using a simple algorithm.
 * The Arrangement parameter corresponds to an arrangement instantiation
 * of type Arrangement_on_surface_2<GeomTraits, TopTraits>.
 */
template <typename Arrangement_>
class Arr_simple_point_location {
public:
  typedef Arrangement_                                   Arrangement_2;
  typedef typename Arrangement_2::Geometry_traits_2      Geometry_traits_2;
  typedef typename Arrangement_2::Topology_traits        Topology_traits;

  typedef typename Arrangement_2::Vertex_const_handle    Vertex_const_handle;
  typedef typename Arrangement_2::Halfedge_const_handle  Halfedge_const_handle;
  typedef typename Arrangement_2::Face_const_handle      Face_const_handle;

  typedef typename Geometry_traits_2::Point_2            Point_2;
  typedef typename Geometry_traits_2::X_monotone_curve_2 X_monotone_curve_2;

#if CGAL_POINT_LOCATION_VERSION < 2
  typedef CGAL::Object                                   result_type;
#else
  typedef typename boost::variant<Vertex_const_handle,
                                  Halfedge_const_handle,
                                  Face_const_handle>     variant_type;
  typedef typename boost::optional<variant_type>         result_type;
#endif

protected:
  typedef typename Topology_traits::Dcel                 Dcel;
  typedef Arr_traits_basic_adaptor_2<Geometry_traits_2>  Traits_adaptor_2;

  // This function returns either make_object() or a result_type constructor
  // to generate return values. The Object version takes a dummy template
  // argument, which is needed for the return of the other option, e.g.,
  // boost::optional<boost::variant> >.
  // In theory a one parameter variant could be returned, but this _could_
  // lead to conversion overhead, and so we rather go for the real type.
  // Overloads for empty returns are also provided.
#if CGAL_POINT_LOCATION_VERSION < 2
  template<typename T>
  inline CGAL::Object result_return(T t) const { return CGAL::make_object(t); }

  inline CGAL::Object result_return() const { return CGAL::Object(); }
#else
  template<typename T>
  inline result_type result_return(T t) const { return result_type(t); }

  inline result_type result_return() const { return result_type(); }
#endif // CGAL_POINT_LOCATION_VERSION < 2
  
  // Data members:
  const Arrangement_2*    p_arr;        // The associated arrangement.  
  const Traits_adaptor_2* geom_traits;  // Its associated geometry traits.
  const Topology_traits*  top_traits;   // Its associated topology traits.

public:

  /*! Default constructor. */
  Arr_simple_point_location() : 
    p_arr(NULL),
    geom_traits(NULL),
    top_traits(NULL)
  {}
        
  /*! Constructor given an arrangement. */
  Arr_simple_point_location(const Arrangement_2& arr) :
    p_arr(&arr)
  {
    geom_traits =
      static_cast<const Traits_adaptor_2*>(p_arr->geometry_traits());
    top_traits = p_arr->topology_traits();
  }

  /*! Attach an arrangement object. */
  void attach(const Arrangement_2& arr) 
  {
    p_arr = &arr;
    geom_traits =
      static_cast<const Traits_adaptor_2*>(p_arr->geometry_traits());
    top_traits = p_arr->topology_traits();
  }

  /*! Detach from the current arrangement object. */
  void detach()
  {
    p_arr = NULL;
    geom_traits = NULL;
    top_traits = NULL;
  }
 
  /*!
   * Locate the arrangement feature containing the given point.
   * \param p The query point.
   * \return An object representing the arrangement feature containing the
   *         query point. This object is either a Face_const_handle or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type locate(const Point_2& p) const;

  /*!
   * Locate the arrangement feature which a upward vertical ray emanating from
   * the given point hits.
   * \param p The query point.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either an empty object or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type ray_shoot_up(const Point_2& p) const
  {
    return (_vertical_ray_shoot(p, true));
  }

  /*!
   * Locate the arrangement feature which a downward vertical ray emanating
   * from the given point hits.
   * \param p The query point.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either an empty object or a
   *         Halfedge_const_handle or a Vertex_const_handle.
   */
  result_type ray_shoot_down(const Point_2& p) const
  {
    return (_vertical_ray_shoot(p, false));
  }

protected:

  /*!
   * Locate the arrangement feature which a vertical ray emanating from the
   * given point hits (not inculding isolated vertices).
   * \param p The query point.
   * \param shoot_up Indicates whether the ray is directed upward or downward.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either a Halfedge_const_handle,
   *         a Vertex_const_handle or an empty object.
   */
  result_type _base_vertical_ray_shoot(const Point_2& p, bool shoot_up) const;

  /*!
   * Locate the arrangement feature which a vertical ray emanating from the
   * given point hits, considering isolated vertices.
   * \param p The query point.
   * \param shoot_up Indicates whether the ray is directed upward or downward.
   * \return An object representing the arrangement feature the ray hits.
   *         This object is either a Halfedge_const_handle,
   *         a Vertex_const_handle or an empty object.
   */
  result_type _vertical_ray_shoot(const Point_2& p, bool shoot_up) const;

  /*!
   * Find the first halfedge with a given source vertex, when going clockwise
   * from "6 o'clock" around this vertex.
   * \param v The given vertex.
   * \return The first halfedge.
   */
  Halfedge_const_handle _first_around_vertex(Vertex_const_handle v) const;
};

} //namespace CGAL

// The member-function definitions can be found under:
#include <CGAL/Arr_point_location/Arr_simple_point_location_impl.h>

#endif
