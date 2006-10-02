// Copyright (c) 2005  Tel-Aviv University (Israel).
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
// Author(s)     : Efi Fogel          <efif@post.tau.ac.il>

#ifndef CGAL_ARR_POSTSCRIPT_FILE_STREAM_H
#define CGAL_ARR_POSTSCRIPT_FILE_STREAM_H

/*! \file
 * Postscript output stream for the Arrangement_2<Traits,Dcel> class.
 */

#include <CGAL/IO/Postscript_file_stream.h>
#include <CGAL/Arrangement_2.h>

CGAL_BEGIN_NAMESPACE

/*! Global exporter of an Arrangement_2 object to a Postscript stream */
class Arr_postscript_file_stream : public Postscript_file_stream {
public:
  /*! Constructor */
   Arr_postscript_file_stream(double w,double h,
                              leda_string name = "CGAL_arr.ps") :
     Postscript_file_stream(w, h, name),
     m_point_color(RED), m_curve_color(BLUE)
  {}

  /*! Constructor */
  Arr_postscript_file_stream(leda_string name = "CGAL_arr.ps") :
    Postscript_file_stream(name),
    m_point_color(RED), m_curve_color(BLUE)
  {}

  /*! Set the color of the points */
  void set_point_color(Color color) { m_point_color = color; }

  /*! Obtain the color of the points */
  Color point_color() const { return m_point_color; }
  
  /*! Set the color of the curves */
  void set_curve_color(Color color) { m_curve_color = color; }

  /*! Obtain the color of the curves */
  Color curve_color() const { return m_curve_color; }
  
private:
  /*! The color of the points */
  Color m_point_color;

  /*! The color of the curves */
  Color m_curve_color;
};

/*! Export an Arrangement_2 object to a Postscript stream
 * \param ps_stream the Postscript stream
 * \param arr the arrangement
 * \return the Postscript stream
 */
template<typename Traits, typename Dcel> Arr_postscript_file_stream &
operator<<(Arr_postscript_file_stream & ps_stream,
           const Arrangement_2<Traits, Dcel> & arr)
{
  // Draw the curves:
  ps_stream << ps_stream.curve_color();
  typename Arrangement_2<Traits, Dcel>::Edge_const_iterator ei;
  for (ei = arr.edges_begin(); ei != arr.edges_end(); ++ei) {
    ps_stream << ei->curve();
  }

  // Draw the points:
  ps_stream << ps_stream.point_color();
  typename Arrangement_2<Traits, Dcel>::Vertex_const_iterator vi;
  for (vi = arr.vertices_begin(); vi != arr.vertices_end(); ++vi) {
    ps_stream << vi->point();
  }
  return ps_stream;
}

CGAL_END_NAMESPACE

#endif
