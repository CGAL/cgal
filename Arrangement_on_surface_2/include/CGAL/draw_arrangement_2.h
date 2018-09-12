// Copyright (c) 2006,2007,2009,2010,2011 Tel-Aviv University (Israel).
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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_ARRANGEMENT_2_H
#define CGAL_DRAW_ARRANGEMENT_2_H

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Random.h>

namespace CGAL
{
  
// Viewer class for Arr 
template <class Arr>
class SimpleArrangementViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename Arr::Dart_const_handle Dart_const_handle;
  typedef typename Arr::Traits Kernel;
  typedef typename Kernel::Point Point;
  typedef typename Kernel::Vector Vector;
  
public:
  /// Construct the viewer.
  /// @param alcc the lcc to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleLCCViewerQt(QWidget* parent,
                    const Arr& a_arr,
                    const char* title="Basic Arrangement Viewer",
                    bool anofaces=false) :
    // First draw: vertices; edges, faces; multi-color; inverse normal
    Base(parent, title, true, true, true, true, true), 
    arr(a_arr),
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
    
    CGAL::Color c=m_fcolor.run(lcc, dh);
    face_begin(c);

    cur=dh;
    do
    {
      add_point_in_face(lcc.point(cur),
                        Geom_utils<LCC>::get_vertex_normal(lcc, cur));
      cur=lcc.next(cur);
    }
    while(cur!=dh);

    face_end();
  }

  void compute_edge(Dart_const_handle dh)
  {
    Point p1 = lcc.point(dh);
    Dart_const_handle d2 = lcc.other_extremity(dh);
    if (d2!=NULL)
    { add_segment(p1, lcc.point(d2)); }
  }

  void compute_vertex(Dart_const_handle dh)
  { add_point(lcc.point(dh)); }

  void compute_elements()
  {
    clear();

    typename LCC::size_type markfaces    = lcc.get_new_mark();
    typename LCC::size_type markedges    = lcc.get_new_mark();
    typename LCC::size_type markvertices = lcc.get_new_mark();

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
        compute_vertex(it);
        CGAL::mark_cell<LCC, 0>(lcc, it, markvertices);
      }
    }

    lcc.free_mark(markfaces);
    lcc.free_mark(markedges);
    lcc.free_mark(markvertices);
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    // Test key pressed:
    //    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    //    if ((e->key()==Qt::Key_PageUp) && (modifiers==Qt::NoButton)) { ... }
    
    // Call: * compute_elements() if the model changed, followed by
    //       * redraw() if some viewing parameters changed that implies some
    //                  modifications of the buffers
    //                  (eg. type of normal, color/mono)
    //       * update() just to update the drawing

    // Call the base method to process others/classicals key
    Base::keyPressEvent(e);
  }

protected:
  const Arr& arr;
  bool m_nofaces;
};
  
template<class GeomTraits_, class TopTraits_>
void draw(const Arrangement_on_surface_2<GeomTraits_, TopTraits_>& a_arr,
          const char* title,
          bool nofill)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"arr_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimpleArrangementViewerQt<Arrangement_on_surface_2<GeomTraits_, TopTraits_> >
      mainwindow(app.activeWindow(), a_arr, title);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_DRAW_ARRANGEMENT_2_H

#endif // CGAL_DRAW_LCC_H
