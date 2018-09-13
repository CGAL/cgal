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

//#include <CGAL/Arrangement_on_surface_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Random.h>

namespace CGAL
{
  
// Viewer class for Arr 
template <class Arr>
class SimpleArrangementViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename Arr::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Arr::Face_const_handle Face_const_handle;
  typedef typename Arr::Geometry_traits_2 Kernel;
  typedef typename Kernel::Point_2 Point;
  typedef typename Kernel::Vector_2 Vector;
  
public:
  /// Construct the viewer.
  /// @param a_a the arrangement to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleArrangementViewerQt(QWidget* parent,
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
  void print_ccb (typename Arr::Ccb_halfedge_const_circulator circ)
  {
    typename Arr::Ccb_halfedge_const_circulator curr = circ;
    std::cout << "(" << curr->source()->point() << ")";
    do
    {
      Halfedge_const_handle he = curr;
      /*    std::cout << " [" << he->curve() << "] "
            << "(" << he->target()->point() << ")";*/
    }
    while (++curr != circ);
  }
  void compute_face(Face_const_handle fh)
  {
    if (fh->is_unbounded())
    { return; }

    // // Print the isolated vertices.
    // typename Arr_2::Isolated_vertex_const_iterator iv;
    // for (iv = fh->isolated_vertices_begin();
    //      iv != fh->isolated_vertices_end(); ++iv)
    // {
    //   iv->point();
    // }
    
    CGAL::Random random((unsigned long)(&*fh));
    CGAL::Color c=get_random_color(random);
    
    face_begin(c);

    print_ccb (fh->outer_ccb());
    typename Arr::Hole_const_iterator hi;
    for (hi=fh->holes_begin(); hi!=fh->holes_end(); ++hi)
    { print_ccb (*hi); }

    /*    cur=dh;
    do
    {
      add_point_in_face(lcc.point(cur));
      cur=lcc.next(cur);
    }
    while(cur!=dh);*/

    face_end();
  }

  /*  void compute_edge(Dart_const_handle dh)
  {
    add_segment(p1, p2);
    } */

  void compute_elements()
  {
    clear();
    
    // Draw the arrangement vertices.
    typename Arr::Vertex_const_iterator vit;    
    for (vit=arr.vertices_begin(); vit!=arr.vertices_end(); ++vit)
    {
      add_point(vit->point());
    }

    // Draw the arrangement edges.
    typename Arr::Edge_const_iterator eit;
    for (eit=arr.edges_begin(); eit!=arr.edges_end(); ++eit)
    { std::cout << "[" << eit->curve() << "]" << std::endl; }

    // Draw the arrangement faces.
    typename Arr::Face_const_iterator fit;
    for (fit=arr.faces_begin(); fit!=arr.faces_end(); ++fit)
    { compute_face(fit); }
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
void draw(const Arrangement_2<GeomTraits_, TopTraits_>& a_arr,
          const char* title="Basic Arrangement Viewer",
          bool nofill=false)
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
    SimpleArrangementViewerQt<Arrangement_2<GeomTraits_, TopTraits_> >
      mainwindow(app.activeWindow(), a_arr, title);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_DRAW_ARRANGEMENT_2_H

#endif // CGAL_DRAW_LCC_H
