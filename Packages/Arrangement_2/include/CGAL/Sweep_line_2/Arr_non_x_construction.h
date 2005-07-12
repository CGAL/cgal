// Copyright (c) 1997  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>

/*! \file
 * Definition of the Arr_non_x_construction<Arrangement> class.
 */

#ifndef ARR_NON_X_CONSTRUCTION_H
#define ARR_NON_X_CONSTRUCTION_H

#include <CGAL/Basic_sweep_line_2.h>
#include <CGAL/Sweep_line_2/Arr_construction_event.h>
#include <CGAL/Sweep_line_2/Arr_construction_curve.h>
#include <CGAL/Sweep_line_2/Arr_construction_visitor.h>
#include <CGAL/assertions.h>
#include <list>
#include <algorithm>

CGAL_BEGIN_NAMESPACE

/*! \class
 * An auxiliray class for performing aggragated insertion of a range of 
 * non-intersecting x-monotone curves into an arrangement using the sweep-line
 * algorithm.
 */
template <class Arrangement_>
class Arr_non_x_construction 
{
  typedef Arrangement_                              Arrangement_2;
  typedef typename Arrangement_2::Halfedge_handle   Halfedge_handle;
  typedef typename Arrangement_2::Edge_iterator     Edge_iterator;
  typedef typename Arrangement_2::Vertex_handle     Vertex_handle;
  typedef typename Arrangement_2::Vertex_iterator   Vertex_iterator;
  typedef typename Arrangement_2::Traits_2          Traits_2;
  typedef Arr_construction_curve<Traits_2>          Subcurve; 
  typedef Arr_construction_event<Traits_2,
                                 Subcurve,
                                 Halfedge_handle>   Event;
  typedef typename Traits_2::X_monotone_curve_2     X_monotone_curve_2;
  typedef typename Traits_2::Point_2                Point_2;
 
  typedef Arr_construction_visitor<Traits_2,
                                   Arrangement_2,
                                   Event,
                                   Subcurve>        Visitor;

  
 
  typedef Basic_sweep_line_2<Traits_2,
                             Visitor,
                             Subcurve,
                             Event>                 Basic_sweep_line_2;
 
public:

  /*! Constructor. */
  Arr_non_x_construction (Arrangement_2& arr) :
    m_arr (&arr),
    m_traits (arr.get_traits()),
    m_visitor (&arr),
    m_sweep_line (m_traits, &m_visitor)
  {}

  /*!
   * Insert a range of x-monotone curves into the arrangement.
   * \param begin An iterator for the first x-monotone curve in the range.
   * \param end A past-the-end iterator for the range.
   * \pre The value-type of the iterators should be X_monotone_curve_2.
   */ 
  template<class XCurveInputIterator>
  void insert_curves (XCurveInputIterator begin,
                      XCurveInputIterator end)
  {
    // Copy the x-montone curves.
    std::list<X_monotone_curve_2>      x_curves;

    std::copy (begin, end, std::back_inserter(x_curves));

    // Add the existing curves in the arrangement.
    Edge_iterator                      eit;

    for (eit = m_arr->edges_begin(); eit != m_arr->edges_end(); ++eit) 
      x_curves.push_back (eit->curve());

    // Add the existing isolated vertices in the arrangement.
    std::list<Point_2>                 iso_points;
    Vertex_iterator                    vit;

    for (vit = m_arr->vertices_begin(); vit != m_arr->vertices_end(); ++vit)
    {
      if (vit->is_isolated())
        iso_points.push_back (vit->point());
    }

    // Perform the sweep.
    m_arr->clear();
    m_sweep_line.sweep (x_curves.begin(),
			                  x_curves.end(),
			                  iso_points.begin(),
			                  iso_points.end());

    return;
  }
              
protected:

  Arrangement_2       *m_arr;
  Traits_2            *m_traits;
  Visitor              m_visitor;
  Basic_sweep_line_2   m_sweep_line;
};

CGAL_END_NAMESPACE

#endif
