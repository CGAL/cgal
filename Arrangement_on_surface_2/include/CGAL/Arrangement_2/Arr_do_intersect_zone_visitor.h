// Copyright (c) 2005,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ophir Setter      <ophirset@post.tau.ac.il>
//
#ifndef CGAL_ARR_DO_INTERSECT_ZONE_VISITOR_H
#define CGAL_ARR_DO_INTERSECT_ZONE_VISITOR_H

/*! \file
 * Definition of the Arr_do_intersect_zone_visitor_2 class.
 */

namespace CGAL {

/*! \class
 * A visitor class for Arrangement_zone_2, which check whether
 * a given x-monotone curve intersects the arrangment.
 * The class shouldbe templated by an Arrangement_2 class.
 */
template <class Arrangement_>
class Arr_do_intersect_zone_visitor
{
public:

  typedef Arrangement_                                Arrangement_2;

  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_2::Face_handle         Face_handle;

  typedef typename Arrangement_2::Point_2             Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2  X_monotone_curve_2;

  typedef std::pair<Halfedge_handle, bool>            Result;

private:

  const Halfedge_handle         invalid_he;    // Invalid halfedge.
  const Vertex_handle           invalid_v;     // Invalid vertex.

  bool                          m_intersect;   // Boolean to hold the answer.

public:

  /*! Constructor. */
  Arr_do_intersect_zone_visitor () :
    invalid_he (),
    invalid_v (),
    m_intersect (false)
  {}

  /*! Initialize the visitor with an arrangement object. */
  void init (Arrangement_2 *)
  {
    m_intersect = false;
  }

  /*!
   * Handle the a subcurve located in the interior of a given face.
   * \param cv The subcurve.
   * \param face The face containing cv's interior.
   * \param left_v The vertex that corresponds to the left endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param left_he The halfedge that contains the left endpoint of cv
   *               (or an invalid handle if no such halfedge exists).
   * \param right_v The vertex that corresponds to the right endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param right_he The halfedge that contains the right endpoint of cv
   *                 (or an invalid handle if no such halfedge exists).
   * \return A handle to the halfedge obtained from the insertion of the
   *         subcurve into the arrangement.
   */
  Result found_subcurve (const X_monotone_curve_2&,
                         Face_handle,
                         Vertex_handle left_v, Halfedge_handle left_he,
                         Vertex_handle right_v, Halfedge_handle right_he)
  { 
    if ((left_v == invalid_v) && (right_v == invalid_v) &&
        (left_he == invalid_he) && (right_he == invalid_he))
    {
      // The current subcurve just lies inside the given face, and its
      // endpoints are not incident to any valid vertex or edge, so it does
      // not intersect the arrangement.
      return (Result (invalid_he, false));
    }

    // We found an intersection. Note we return a result indicating that the
    // zone-computation can stop here.
    m_intersect = true;
    return (Result (invalid_he, true));
  }

  /*!
   * Handle the a subcurve that overlaps a given edge.
   * \param cv The overlapping subcurve.
   * \param he The overlapped halfedge (directed from left to right).
   * \param left_v The vertex that corresponds to the left endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \param right_v The vertex that corresponds to the right endpoint of cv
   *               (or an invalid handle if no such arrangement vertex exists).
   * \return A handle to the halfedge obtained from the insertion of the
   *         overlapping subcurve into the arrangement.
   */
  Result found_overlap (const X_monotone_curve_2&,
                        Halfedge_handle,
                        Vertex_handle, Vertex_handle)
  {
    // We found an overlap (hence an intersection). Note we return a result
    // indicating that the zone-computation can stop here.
    m_intersect = true;
    return (Result (invalid_he, true));
  }

  bool do_intersect () const
  {
    return (m_intersect);
  }
};

} //namespace CGAL

#endif
