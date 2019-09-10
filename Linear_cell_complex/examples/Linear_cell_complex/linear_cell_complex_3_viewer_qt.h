// Copyright (c) 2011 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0+
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_LCC_3_VIEWER_QT_H
#define CGAL_LCC_3_VIEWER_QT_H

#include "basic_viewer.h"
#include <CGAL/Linear_cell_complex.h>
#include <CGAL/Cartesian_converter.h>
#include <CGAL/Random.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3  Local_point;
typedef Local_kernel::Vector_3 Local_vector;

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctor
{
  template<typename LCC>
  static CGAL::Color run(const LCC& alcc,
                         typename LCC::Dart_const_handle dh)
  {
    if (dh==alcc.null_handle) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    // Here dh is the smaller dart of its face
    CGAL::Random random(alcc.darts().index(dh));
    CGAL::Color res;
    do
    {
      res=CGAL::Color(random.get_int(0,256),
                      random.get_int(0,256),
                      random.get_int(0,256));
    }
    while(res.red()==255 && res.green()==255 && res.blue()==255);
    return res;
  }
};

template<class LCC, int dim=LCC::ambient_dimension>
struct Geom_utils;

template<class LCC>
struct Geom_utils<LCC,3>
{
  Local_point get_point(const LCC& lcc,
                        typename LCC::Vertex_attribute_const_handle vh)
  { return converter(lcc.point_of_vertex_attribute(vh)); }

  Local_point get_point(const LCC& lcc, typename LCC::Dart_const_handle dh)
  { return converter(lcc.point(dh)); }

  Local_vector get_vertex_normal(const LCC& lcc,
                                 typename LCC::Dart_const_handle dh)
  {
    Local_vector n = converter(CGAL::compute_normal_of_cell_0<LCC>(lcc,dh));
    n = n/(CGAL::sqrt(n*n));
    return n;
  }
protected:
  CGAL::Cartesian_converter<typename LCC::Traits, Local_kernel> converter;
};

template<class LCC>
struct Geom_utils<LCC,2>
{
  Local_point get_point(const LCC& lcc,
                        typename LCC::Vertex_attribute_const_handle vh)
  {
    Local_point p(converter(lcc.point_of_vertex_attribute(vh).x()),0,
                  converter(lcc.point_of_vertex_attribute(vh).y()));
    return p;
  }

  Local_point get_point(const LCC& lcc, typename LCC::Dart_const_handle dh)
  { return get_point(lcc, lcc.vertex_attribute(dh)); }

  Local_vector get_vertex_normal(const LCC&, typename LCC::Dart_const_handle)
  {
    Local_vector n(0,-1,0);
    return n;
  }
protected:
  CGAL::Cartesian_converter<typename LCC::Traits, Local_kernel> converter;
};

// Viewer class for LCC 
template<class LCC, class ColorFunctor>
class SimpleLCCViewerQt : public Basic_viewer
{
  typedef Basic_viewer Base;
  typedef typename LCC::Dart_const_handle Dart_const_handle;
  
public:
  /// Construct the viewer.
  /// @param alcc the lcc to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big LCC where this time could be long)
  SimpleLCCViewerQt(const LCC& alcc, const char* title="", bool anofaces=false) :
    Base(title),
    lcc(alcc),
    m_nofaces(anofaces)
  {
    compute_elements();
  }

protected:
  
  void compute_face(Dart_const_handle dh)
  {
    // We fill only closed faces.
    Dart_const_handle cur=dh;
    Dart_const_handle min=dh;
    do
    {
      if (!lcc.is_next_exist(cur)) return; // open face=>not filled
      if (cur<min) min=cur;
      cur=lcc.next(cur);
    }
    while(cur!=dh);
    
    CGAL::Color c=ColorFunctor::run(lcc, dh);

    if (c.red()<60 || c.green()<60 || c.blue()<60)
      mono_face_begin();
    else
      colored_face_begin(c);

    cur=dh;
    do
    {
      add_point_in_face(geomutils.get_point(lcc, cur),
                        geomutils.get_vertex_normal(lcc, cur));
      cur=lcc.next(cur);
    }
    while(cur!=dh);

    face_end();
  }

  void compute_edge(Dart_const_handle dh)
  {
    Local_point p1 = geomutils.get_point(lcc, dh);
    Dart_const_handle d2 = lcc.other_extremity(dh);
    if ( d2!=NULL )
    {
      Local_point p2 = geomutils.get_point(lcc, d2);
      add_mono_segment(p1, p2);
    }
  }

  void compute_vertex(Dart_const_handle dh, bool& empty)
  {
    Local_point p = geomutils.get_point(lcc, dh);
    add_mono_point(p);
  }

  void compute_elements()
  {
    clear();
    
    unsigned int markfaces    = lcc.get_new_mark();
    unsigned int markedges    = lcc.get_new_mark();
    unsigned int markvertices = lcc.get_new_mark();

    bool empty = true;
    
    for (typename LCC::Dart_range::const_iterator it=lcc.darts().begin(),
         itend=lcc.darts().end(); it!=itend; ++it )
    {
      if ( !m_nofaces && !lcc.is_marked(it, markfaces) )
      {
        compute_face(it);
        CGAL::mark_cell<LCC, 2>(lcc, it, markfaces);
      }

      if ( !lcc.is_marked(it, markedges) )
      {
        compute_edge(it);
        CGAL::mark_cell<LCC, 1>(lcc, it, markedges);
      }

      if ( !lcc.is_marked(it, markvertices) )
      {
        compute_vertex(it, empty);
        CGAL::mark_cell<LCC, 0>(lcc, it, markvertices);
      }
    }

    lcc.free_mark(markfaces);
    lcc.free_mark(markedges);
    lcc.free_mark(markvertices);
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    Base::keyPressEvent(e);
  }

protected:
  const LCC& lcc;
  bool m_nofaces;
  Geom_utils<LCC> geomutils;
};

  
template<class LCC, class ColorFunctor=DefaultColorFunctor>
void display_lcc(const LCC& alcc,
                 const char* title="",
                 bool nofill=false)
{
  int argc=1;

  const char* argv[2]={"lccviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimpleLCCViewerQt<LCC, ColorFunctor> mainwindow(alcc, title, nofill);
  mainwindow.show();

  app.exec();
}

#endif // CGAL_LCC_3_VIEWER_QT_H
