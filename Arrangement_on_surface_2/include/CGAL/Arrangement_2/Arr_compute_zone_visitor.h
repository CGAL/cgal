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
#ifndef CGAL_ARR_COMPUTE_ZONE_VISITOR_H
#define CGAL_ARR_COMPUTE_ZONE_VISITOR_H

#include <CGAL/license/Arrangement_on_surface_2.h>


/*! \file
 * Definition of the Arr_compute_zone_visitor class.
 */

namespace CGAL {

/*! \class
 * A visitor class for Arrangement_zone_2, which outputs the
 * zone of an x-monotone curve. Meaning, it output the arrangment's
 * vertices, edges and faces that the x-monotone curve intersects.
 * The class should be templated by an Arrangement_2 class, and by an
 * output iterator of CGAL Objects, where we store all arrangement
 * features the x-monotone curve intersects.
 */
template <class Arrangement_, class OutputIterator_>
class Arr_compute_zone_visitor
{
public:
  
  typedef OutputIterator_                             OutputIterator;
  typedef Arrangement_                                Arrangement_2;

  typedef typename Arrangement_2::Vertex_handle       Vertex_handle;
  typedef typename Arrangement_2::Halfedge_handle     Halfedge_handle;
  typedef typename Arrangement_2::Face_handle         Face_handle;

  typedef typename Arrangement_2::Point_2             Point_2;
  typedef typename Arrangement_2::X_monotone_curve_2  X_monotone_curve_2;

  typedef std::pair<Halfedge_handle, bool>            Result;

private:

  const Halfedge_handle      invalid_he;   // Invalid halfedge.
  const Vertex_handle        invalid_v;    // Invalid vertex.

  OutputIterator&            out_iter;     // for outputing the zone objects.
                                           // Its value type is CGAL::Object.
  bool                       output_left;  // Determines wheter we should
                                           // output the left end point of a
                                           // subcurve (to avoid outputing
                                           // the same feature twice).

public:

  /*! Constructor. */
  Arr_compute_zone_visitor (OutputIterator& oi) :
    invalid_he(),
    invalid_v(),
    out_iter (oi),
    output_left (true)
  {}

  /*! Initialize the visitor. */
  void init (Arrangement_2 *)
  {
    output_left = true;
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
                         Face_handle face,
                         Vertex_handle left_v, Halfedge_handle left_he,
                         Vertex_handle right_v, Halfedge_handle right_he)
  {
    if (output_left)
    {
      // Only the first subcurve should output the arrangement feature incident
      // to its left endpoint. This way we avoid reporting the same feature
      // twice.
      if (left_v != invalid_v)
      {
        *out_iter = CGAL::make_object (left_v);
        ++out_iter;
      }
      else if (left_he != invalid_he)
      {
        *out_iter = CGAL::make_object (left_he);
        ++out_iter;
      }

      output_left = false;
    }

    // Report the face that contains the interior of the subcurve.
    *out_iter = CGAL::make_object (face);
    ++out_iter;

    // If the right endpoint of the subcurve is incident to an arrangement
    // vertex or an arrangement edge, report this feature.
    if (right_v != invalid_v)
    {
      *out_iter = CGAL::make_object(right_v);
      ++out_iter;
    }
    else if (right_he != invalid_he)
    {
      *out_iter = CGAL::make_object(right_he);
      ++out_iter;
    }

    // We did not modify the arrangement, so we return an invalid handle
    // and a flag indicating that the zone-computation process should continue.
    return (Result (invalid_he, false));
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
                        Halfedge_handle he,
                        Vertex_handle left_v, Vertex_handle right_v)
  {
    if (output_left)
    {
      // Only the first subcurve should output the arrangement feature incident
      // to its left endpoint. This way we avoid reporting the same feature
      // twice.
      if (left_v != invalid_v)
      {
        *out_iter = CGAL::make_object (left_v);
        ++out_iter;
      }

      output_left = false;
    }

    // Report the arrangement edge the curve currently overlaps.
    *out_iter = CGAL::make_object(he);
    ++out_iter;

    // If the right endpoint of the overlapping subcurve is incident to an
    // arrangement vertex, report this vertex as well.
    if (right_v != invalid_v)
    {
      *out_iter = CGAL::make_object (right_v);
      ++out_iter;
    }

    // We did not modify the arrangement, so we return an invalid handle
    // and a flag indicating that the zone-computation process should continue.
    return (Result (invalid_he, false));
  }
};


} //namespace CGAL

#endif
