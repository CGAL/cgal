// Copyright (c) 2006  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Ron Wein   <wein_r@yahoo.com>

#ifndef CGAL_UNION_OF_SEGMENT_CYCLES_2_H
#define CGAL_UNION_OF_SEGMENT_CYCLES_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Minkowski_sum_2/Union_of_cycles_2.h>

namespace CGAL {

/*! \class
 * An auxiliary class for computing the union of the interiors of segment
 * cycles (polygon boundaries or convolution cycles).
 */
template <class Traits_, class Polygon_>
class Union_of_segment_cycles_2 : private Union_of_cycles_2<Traits_>
{
public:

  typedef Traits_                                 Traits_2;
  typedef Polygon_                                Polygon_2;

private:

  // Base-class types:
  typedef Union_of_cycles_2<Traits_2>             Base;
  typedef typename Base::Point_2                  Point_2;
  typedef typename Base::X_monotone_curve_2       X_monotone_curve_2;

  typedef typename Base::Arrangement_2            Arrangement_2;
  typedef typename Base::Vertex_handle            Vertex_handle;
  typedef typename Base::Halfedge_handle          Halfedge_handle;
  typedef typename Base::Face_handle              Face_handle;
  typedef typename Base::Vertex_iterator          Vertex_iterator;
  typedef typename Base::Edge_iterator            Edge_iterator;
  typedef typename Base::Halfedge_iterator        Halfedge_iterator;
  typedef typename Base::Face_iterator            Face_iterator;
  typedef typename Base::Inner_ccb_iterator       Inner_ccb_iterator;
  typedef typename Base::Halfedge_around_vertex_circulator
                                             Halfedge_around_vertex_circulator;
  typedef typename Base::Ccb_halfedge_circulator  Ccb_halfedge_circulator;

public:

  /*! Default constructor. */
  Union_of_segment_cycles_2 () :
    Base()
  {}

  /*!
   * Compute the union of the interiors of the segment cycles.
   * \param begin An iterator for the first segment in the range.
   * \param end A past-the-end iterator for the segment range.
   * \param out_bound Output: A polygon representing the union boundary.
   * \param holes Output: An output iterator of the holes in the union.
   * \return A past-the-end iterator for the holes.
   */
  template <class InputIterator, class OutputIterator>
  OutputIterator operator() (InputIterator begin, InputIterator end,
                             Polygon_2& out_bound,
                             OutputIterator holes) const
  {
    // Construct the arrangement of all segments.
    Arrangement_2                    arr;

    this->_construct_arrangement (begin, end, arr);

    // Produce the result. First set the outer boundary of the union, given
    // as the inner boundary of the single hole in the unbounded face.
    Face_iterator                    fit;
    const Face_handle                uf = arr.unbounded_face();
    Inner_ccb_iterator               iccb_it = uf->inner_ccbs_begin();
    Ccb_halfedge_circulator          first, circ;
    Halfedge_handle                  he;

    out_bound.erase (out_bound.vertices_begin(), out_bound.vertices_end());

    circ = first = *iccb_it;
    do
    {
      out_bound.push_back (circ->source()->point());
      --circ;

    } while (circ != first);
    ++iccb_it;

    // Locate the holes in the union: Go over all arrangement faces.
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    {
      CGAL_assertion (fit->data() != this->UNVISITED);

      // If a bounded face has an inside count that equals 0, it forms a hole
      // in the union.
      if (! fit->is_unbounded() && fit->data() == 0)
      {
        Polygon_2   pgn_hole;

        circ = first = fit->outer_ccb();
        do
        {
          pgn_hole.push_back (circ->source()->point());
          --circ;

        } while (circ != first);

        // Insert it to the containers of holes in the Minkowski sum.
        *holes = pgn_hole;
        ++holes;
      }
    }

    return (holes);
  }

};

} //namespace CGAL

#endif
