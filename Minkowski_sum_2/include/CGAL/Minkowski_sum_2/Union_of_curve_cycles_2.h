// Copyright (c) 2006  Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein   <wein_r@yahoo.com>

#ifndef CGAL_UNION_OF_CURVE_CYCLES_2_H
#define CGAL_UNION_OF_CURVE_CYCLES_2_H

#include <CGAL/license/Minkowski_sum_2.h>


#include <CGAL/Minkowski_sum_2/Union_of_cycles_2.h>

namespace CGAL {

/*! \class
 * An auxiliary class for computing the union of the interiors of curve
 * cycles (convolution cycles).
 */
template <class Traits_, class GeneralPolygon_>
class Union_of_curve_cycles_2 : private Union_of_cycles_2<Traits_>
{
public:

  typedef Traits_                                 Traits_2;
  typedef GeneralPolygon_                         General_polygon_2;

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

  /*! \class
   * A circulator adapter that iterates over all curves in a CCB.
   */
  class Ccb_curve_iterator
  {
  public:

    typedef Ccb_curve_iterator                      Self;
    typedef typename Union_of_cycles_2<Traits_>::Ccb_halfedge_circulator
                                                    Ccb_halfedge_circulator;
    typedef typename Union_of_cycles_2<Traits_>::X_monotone_curve_2
                                                    value_type;
    typedef std::forward_iterator_tag               iterator_category;
    typedef const value_type&                       reference;
    typedef const value_type*                       pointer;
    typedef int                                     difference_type;

  private:

    Ccb_halfedge_circulator    _first;    // The first halfedge.
    Ccb_halfedge_circulator    _circ;     // The current circulator.
    bool                       _done;     // Indicates whether we completed
                                          // a full traversal of the CCB.

  public:

    /*! Default constructor. */
    Ccb_curve_iterator () :
      _done (true)
    {}

    /*!
     * Constructor from a circulator.
     * \param circ A circulator for the first halfedge in the CCB.
     * \param done (true) in order to create a past-the-end iterator.
     */
    Ccb_curve_iterator (Ccb_halfedge_circulator circ,
                        bool done = false) :
      _first (circ),
      _circ (circ),
      _done (done)
    {}

    /*! Dereference operators. */
    reference operator* () const
    {
      return (_circ->curve());
    }

    pointer operator-> () const
    {
      return (&(_circ->curve()));
    }

    /*! Equality operators.*/
    bool operator== (const Self& it) const
    {
      return (_done == it._done && _circ == it._circ);
    }

    bool operator!= (const Self& it) const
    {
      return (_done != it._done || _circ != it._circ);
    }

    /*! Increment operators. */
    Self& operator++ ()
    {
      if (! _done)
      {
        --_circ;

        if (_circ == _first)
          _done = true;
      }

      return (*this);
    }

    Self operator++ (int )
    {
      Self   temp = *this;

      if (! _done)
      {
        --_circ;

        if (_circ == _first)
          _done = true;
      }

      return (temp);
    }

  };

  /*! Default constructor. */
  Union_of_curve_cycles_2 () :
    Base()
  {}

  /*!
   * Compute the union of the interiors of the curve cycles.
   * \param begin An iterator for the first curve in the range.
   * \param end A past-the-end iterator for the curve range.
   * \param out_bound Output: A generalized polygon representing the outer
   *                          boundary of the union.
   * \param holes Output: An output iterator of the holes in the union.
   * \return A past-the-end iterator for the holes.
   */
  template <class InputIterator, class OutputIterator>
  OutputIterator operator() (InputIterator begin, InputIterator end,
                             General_polygon_2& out_bound,
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
    Ccb_halfedge_circulator          circ;

    circ = *iccb_it;
    out_bound = General_polygon_2 (Ccb_curve_iterator (circ),
                                   Ccb_curve_iterator (circ, true));

    // Locate the holes in the union: Go over all arrangement faces.
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    {
      CGAL_assertion (fit->data() != this->UNVISITED);

      // If a bounded face has an inside count that equals 0, it forms a hole
      // in the union.
      if (! fit->is_unbounded() && fit->data() == 0)
      {
        circ = fit->outer_ccb();
        *holes = General_polygon_2 (Ccb_curve_iterator (circ),
                                    Ccb_curve_iterator (circ, true));
        ++holes;
      }
    }

    return (holes);
  }


  /*!
   * Compute the inverse of the union of the interiors of the curve cycles,
   * and return the result as a sequence of polygons.
   * \param begin An iterator for the first curve in the range.
   * \param end A past-the-end iterator for the curve range.
   * \param oi Output: An output iterator of the polygons.
   * \return A past-the-end iterator for the polygons.
   */
  template <class InputIterator, class OutputIterator>
  OutputIterator inverse (InputIterator begin, InputIterator end,
                          OutputIterator oi) const
  {
    // Construct the arrangement of all segments.
    Arrangement_2                    arr;

    this->_construct_arrangement (begin, end, arr);

    // Go over all arrangement faces, and output each face whose data is
    // negative.
    Face_iterator                    fit;
    Ccb_halfedge_circulator          circ;

    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit)
    {
      CGAL_assertion (fit->data() != this->UNVISITED);

      if (fit->data() < 0)
      {
        circ = fit->outer_ccb();
        *oi = General_polygon_2 (Ccb_curve_iterator (circ),
                                 Ccb_curve_iterator (circ, true));
        ++oi;
      }
    }

    return (oi);
  }

};

} //namespace CGAL

#endif
