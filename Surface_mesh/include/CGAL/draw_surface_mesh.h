// Copyright (c) 2018 GeometryFactory (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
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

#include <CGAL/Surface_mesh.h>
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
  Local_vector get_face_normal(halfedge_descriptor he)
  {
    Local_vector normal=CGAL::NULL_VECTOR;
    halfedge_descriptor end=he;
    unsigned int nb=0;
    do
    {
      internal::newell_single_step_3
        (this->get_local_point(sm.point(sm.source(he))),
         this->get_local_point(sm.point(sm.target(he))), normal);
      ++nb;
      he=sm.next(he);
    }
    while (he!=end);
    assert(nb>0);
    return (typename Local_kernel::Construct_scaled_vector_3()(normal, 1.0/nb));
  }
  
  Local_vector get_vertex_normal(halfedge_descriptor he)
  {
    Local_vector normal=CGAL::NULL_VECTOR;
    halfedge_descriptor end=he;
    do
    {
      if (!sm.is_border(he))
      {
        Local_vector n=get_face_normal(he);
        normal=typename Local_kernel::Construct_sum_of_vectors_3()(normal, n);
      }
      he=sm.next(sm.opposite(he));
    }
    while (he!=end);
    
    if (!typename Local_kernel::Equal_3()(normal, CGAL::NULL_VECTOR))
    { normal=(typename Local_kernel::Construct_scaled_vector_3()
              (normal, 1.0/CGAL::sqrt(normal.squared_length()))); }
    
    return normal;
  }

protected:
  const SM& sm;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
};

// Specialization of draw function.
template<class K>
void draw(const Surface_mesh<K>& amesh,
          const char* title="Surface_mesh Basic Viewer",
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
    const char* argv[2]={"surface_mesh_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    DefaultColorFunctorSM fcolor;
    SimpleSurfaceMeshViewerQt<Surface_mesh<K>, DefaultColorFunctorSM>
      mainwindow(app.activeWindow(), amesh, title, nofill, fcolor);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DOXYGEN_RUNNING

#endif // CGAL_DRAW_SURFACE_MESH_H
