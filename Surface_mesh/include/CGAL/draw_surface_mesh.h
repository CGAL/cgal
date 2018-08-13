// Copyright (c) 2018 GeometryFactory (France)
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

#ifndef CGAL_DRAW_SURFACE_MESH_H
#define CGAL_DRAW_SURFACE_MESH_H

#ifdef DOXYGEN_RUNNING

/*!
\ingroup PkgDrawSurfaceMesh

Open a new window and draw `asm`, an instance of the `CGAL::Surface_mesh` class. The function is blocking, that is the program continues as soon as the user closes the window. This function requires CGAL_Qt5, and is only available if the flag CGAL_USE_BASIC_VIEWER is defined at compile time.
\tparam SM an instance of the `CGAL::Surface_mesh` class.
\param asm the surface mesh to draw.

*/
template<class SM>
void draw(const SM& asm);

#else // DOXYGEN_RUNNING
  
#include <CGAL/license/Surface_mesh.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Random.h>

namespace CGAL
{

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorSM
{
  template<typename SM>
  static CGAL::Color run(const SM&,
                         typename SM::Face_index fh)
  {
    if (fh==boost::graph_traits<SM>::null_face()) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    CGAL::Random random((unsigned int)fh);
    return get_random_color(random);
  }
};

template<class SM, class ColorFunctor>
class SimpleSurfaceMeshViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename SM::Point Point;
  typedef typename CGAL::Kernel_traits<Point>::Kernel Kernel;
  typedef typename SM::Vertex_index vertex_descriptor;
  typedef typename SM::Face_index face_descriptor;
  typedef typename SM::Edge_index edge_descriptor;
  typedef typename SM::Halfedge_index halfedge_descriptor;
  
public:
  /// Construct the viewer.
  /// @param amesh the surface mesh to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleSurfaceMeshViewerQt(QWidget* parent,
                            const SM& amesh,
                            const char* title="Basic Surface_mesh Viewer",
                            bool anofaces=false,
                            const ColorFunctor& fcolor=ColorFunctor()) :
    // First draw: no vertex; edges, faces; mono-color; inverse normal
    Base(parent, title, false, true, true, true, false),
    sm(amesh),
    m_nofaces(anofaces),
    m_fcolor(fcolor)
  {
    compute_elements();
  }

protected:
  void compute_face(face_descriptor fh)
  {
    CGAL::Color c=m_fcolor.run(sm, fh);
    face_begin(c);
    halfedge_descriptor hd=sm.halfedge(fh);
    do
    {
      add_point_in_face(sm.point(sm.source(hd)), get_vertex_normal(hd));
      hd=sm.next(hd);
    }
    while(hd!=sm.halfedge(fh));
    face_end();
  }

  void compute_edge(edge_descriptor e)
  {
    add_segment(sm.point(sm.source(sm.halfedge(e))),
                sm.point(sm.target(sm.halfedge(e))));
  } 

  void compute_vertex(vertex_descriptor vh)
  { add_point(sm.point(vh)); }

  void compute_elements()
  {
    clear();

    if (!m_nofaces)
    {
      for (typename SM::Face_range::iterator f=sm.faces().begin();
           f!=sm.faces().end(); ++f)
      {
        if (*f!=boost::graph_traits<SM>::null_face())
        { compute_face(*f); }
      }
    }
    
    for (typename SM::Edge_range::iterator e=sm.edges().begin();
         e!=sm.edges().end(); ++e)
    { compute_edge(*e); }

    for (typename SM::Vertex_range::iterator v=sm.vertices().begin();
         v!=sm.vertices().end(); ++v)
    { compute_vertex(*v); }
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
  typename Kernel::Vector_3 get_face_normal(halfedge_descriptor he)
  {
    typename Kernel::Vector_3 normal=CGAL::NULL_VECTOR;
    halfedge_descriptor end=he;
    unsigned int nb=0;
    do
    {
      internal::newell_single_step_3(sm.point(sm.source(he)),
                                     sm.point(sm.target(he)), normal);
      ++nb;
      he=sm.next(he);
    }
    while (he!=end);
    assert(nb>0);
    return (typename Kernel::Construct_scaled_vector_3()(normal, 1.0/nb));
  }
  
  typename Kernel::Vector_3 get_vertex_normal(halfedge_descriptor he)
  {
    typename Kernel::Vector_3 normal=CGAL::NULL_VECTOR;
    halfedge_descriptor end=he;
    do
    {
      if (!sm.is_border(he))
      {
        typename Kernel::Vector_3 n=get_face_normal(he);
        normal=typename Kernel::Construct_sum_of_vectors_3()(normal, n);
      }
      he=sm.next(sm.opposite(he));
    }
    while (he!=end);
    
    if (!typename Kernel::Equal_3()(normal, CGAL::NULL_VECTOR))
    { normal=(typename Kernel::Construct_scaled_vector_3()
              (normal, 1.0/CGAL::sqrt(normal.squared_length()))); }
    
    return normal;
  }

protected:
  const SM& sm;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
};

template<class SM, class ColorFunctor>
void draw(const SM& amesh,
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
    const char* argv[2]={"surface_mesh_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimpleSurfaceMeshViewerQt<SM, ColorFunctor> mainwindow(app.activeWindow(),
                                                           amesh,
                                                           title,
                                                           nofill,
                                                           fcolor);
    mainwindow.show();
    app.exec();
  }
}

template<class SM>
void draw(const SM& amesh, const char* title, bool nofill)
{
  DefaultColorFunctorSM c;
  draw(amesh, title, nofill, c);
}

template<class SM>
void draw(const SM& amesh, const char* title)
{ draw(amesh, title, false); }

template<class SM>
void draw(const SM& amesh)
{ draw(amesh, "Basic Surface_mesh Viewer"); }

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DOXYGEN_RUNNING

#endif // CGAL_DRAW_SURFACE_MESH_H
