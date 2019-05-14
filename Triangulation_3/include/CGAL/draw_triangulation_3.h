// Copyright (c) 2018  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_T3_H
#define CGAL_DRAW_T3_H

#include <CGAL/license/Triangulation_3.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Random.h>

namespace CGAL
{
  
// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorT3
{
  template<typename T3>
  static CGAL::Color run(const T3&,
                         const typename T3::Finite_facets_iterator* fh)
  {
    if (fh==NULL) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    CGAL::Random random((unsigned int)((std::size_t)(&*((*fh)->first))+
                                       (std::size_t)((*fh)->second)));
    return get_random_color(random);
  }
};

// Viewer class for T3 
template<class T3, class ColorFunctor>
class SimpleTriangulation3ViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt                 Base;
  typedef typename T3::Vertex_handle Vertex_const_handle;
  typedef typename T3::Finite_edges_iterator  Edge_const_handle;
  typedef typename T3::Finite_facets_iterator Facet_const_handle;
  typedef typename T3::Cell_handle         Cell_handle;
  typedef typename T3::Point               Point;
 
public:
  /// Construct the viewer.
  /// @param at3 the t3 to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleTriangulation3ViewerQt(QWidget* parent,
                               const T3& at3,
                               const char* title="Basic T3 Viewer",
                               bool anofaces=false,
                               const ColorFunctor& fcolor=ColorFunctor()) :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false), 
    t3(at3),
    m_nofaces(anofaces),
    m_fcolor(fcolor)
  {
    compute_elements();
  }

protected:
  void compute_face(Facet_const_handle fh)
  {
    CGAL::Color c=m_fcolor.run(t3, &fh);
    face_begin(c);

    add_point_in_face(fh->first->vertex((fh->second+1)%4)->point());
    add_point_in_face(fh->first->vertex((fh->second+2)%4)->point());
    add_point_in_face(fh->first->vertex((fh->second+3)%4)->point());
    
    face_end();
  }

  void compute_edge(Edge_const_handle eh)
  {
    add_segment(eh->first->vertex(eh->second)->point(),
                eh->first->vertex(eh->third)->point());
  }

  void compute_vertex(Vertex_const_handle vh)
  { add_point(vh->point()); }

  void compute_elements()
  {
    clear();

    if (!m_nofaces)
    {
      for (typename T3::Finite_facets_iterator it=t3.finite_facets_begin();
           it!=t3.finite_facets_end(); ++it)
      { compute_face(it); } 
    }
    
    for (typename T3::Finite_edges_iterator it=t3.finite_edges_begin();
         it!=t3.finite_edges_end(); ++it)
    { compute_edge(it); } 

    for (typename T3::Finite_vertices_iterator it=t3.finite_vertices_begin();
         it!=t3.finite_vertices_end(); ++it)
    { compute_vertex(it); } 
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
  const T3& t3;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
};
  
template<class T3, class ColorFunctor>
void draw(const T3& at3,
          const char* title,
          bool nofill,
          const ColorFunctor& fcolor)
{

#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=false;
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"t3_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimpleTriangulation3ViewerQt<T3, ColorFunctor> mainwindow(app.activeWindow(),
                                                              at3,
                                                              title,
                                                              nofill,
                                                              fcolor);
    mainwindow.show();
    app.exec();
  }
}

template<class T3>
void draw(const T3& at3, const char* title, bool nofill)
{
  DefaultColorFunctorT3 c;
  draw(at3, title, nofill, c);
}

template<class T3>
void draw(const T3& at3, const char* title)
{ draw(at3, title, false); }

template<class T3>
void draw(const T3& at3)
{ draw(at3, "Basic T3 Viewer"); }

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_T3_H
