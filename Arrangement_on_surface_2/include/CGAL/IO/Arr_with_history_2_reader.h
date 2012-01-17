// Copyright (c) 2005,2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
// Author(s)     : Ron Wein           <wein@post.tau.ac.il>
#ifndef CGAL_ARR_WITH_HISTORY_2_READER_H
#define CGAL_ARR_WITH_HISTORY_2_READER_H

/*! \file
 * The header file for the Arr_with_history_2_reader<Arrangement> class.
 */

#include <CGAL/IO/Arrangement_2_reader.h>
#include <CGAL/Arrangement_2/Arr_with_history_accessor.h>

namespace CGAL {

/*! \class
 * An auxiliary class for reading an arrangement with history from an
 * input stream.
 */
template <class ArrWithHistory_>
class Arr_with_history_2_reader : private Arrangement_2_reader<ArrWithHistory_>
{
public:

  typedef ArrWithHistory_                                 Arr_with_history_2;
  typedef Arr_with_history_2_reader<Arr_with_history_2>   Self;

protected:
 
  typedef Arrangement_2_reader<Arr_with_history_2>        Base;
  typedef typename Arr_with_history_2::Size               Size;
  typedef typename Arr_with_history_2::Curve_handle       Curve_handle;
  typedef typename Arr_with_history_2::Halfedge_handle    Halfedge_handle;

  typedef Arr_with_history_accessor<Arr_with_history_2>   Arr_with_hist_access;
  typedef typename Arr_with_history_2::Curve_2            Curve_2;

protected:

  // Data members:
  Curve_2               m_in_curve;
  Arr_with_hist_access  m_arr_with_hist_access;

private:

  // Copy constructor and assignment operator - not supported.
  Arr_with_history_2_reader (const Self& );
  Self& operator= (const Self& );

public:

  /*! Constructor. */
  Arr_with_history_2_reader (Arr_with_history_2& arr) :
    Base (arr),
    m_arr_with_hist_access (arr)
  {}

  /*! Read the arrangement. */
  template <class Formatter>
  void operator()(Formatter& formatter)
  {
    // Read the arrangement (without history).
    Base::operator() (formatter);
    
    // Read the inducing curves.
    formatter.read_curves_begin();

    const Size  number_of_curves = formatter.read_size("number_of_curves");
    Size        k;

    for (k = 0; k < number_of_curves; k++)
      _read_curve (formatter);

    formatter.read_curves_end();
    return;
  }

protected:

  /*! Read a curve with its induced edges. */
  template <class Formatter>
  void _read_curve (Formatter& formatter)
  {
    formatter.read_curve_begin();

    // Read the curve.
    formatter.read_curve (m_in_curve);

    // Insert the curve to the list of inducing curves of the arrangement.
    Curve_handle     new_cv = m_arr_with_hist_access.new_curve (m_in_curve);

    // Read the induced edges.
    formatter.read_induced_edges_begin();

    const Size       number_of_edges = formatter.read_size("induced_edges");
    std::size_t      curr_idx;
    Halfedge_handle  curr_he;
    Size             k;

    for (k = 0; k < number_of_edges; k++)
    {
      curr_idx = formatter.read_halfedge_index();
      curr_he = Halfedge_handle (this->m_halfedges[curr_idx]);

      // Connect the curve and the edge it induces.
      m_arr_with_hist_access.connect_curve_edge (new_cv, curr_he);
    }
    formatter.read_induced_edges_end();

    formatter.read_curve_end();
    return;
  }
   
};

} //namespace CGAL

#endif
