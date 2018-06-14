// Copyright (c) 2018  ETH Zurich (Switzerland).
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

#ifndef CGAL_DRAW_POLYHEDRON_H
#define CGAL_DRAW_POLYHEDRON_H

#include <CGAL/license/Polyhedron.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Random.h>

namespace CGAL
{
  
// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorPolyhedron
{
  template<typename Polyhedron>
  static CGAL::Color run(const Polyhedron&,
                         typename Polyhedron::Facet_const_handle fh)
  {
    if (fh==boost::graph_traits<Polyhedron>::null_face()) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    CGAL::Random random((unsigned int)(std::size_t)(&(*fh)));
    return get_random_color(random);
  }
};

// Viewer class for Polyhedron
template<class Polyhedron, class ColorFunctor>
class SimplePolyhedronViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename Polyhedron::Traits Kernel;
  typedef typename Polyhedron::Halfedge_const_handle Halfedge_const_handle;
  typedef typename Polyhedron::Vertex_const_handle Vertex_const_handle;
  typedef typename Polyhedron::Facet_const_handle Facet_const_handle;
  
public:
  /// Construct the viewer.
  /// @param apoly the polyhedron to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimplePolyhedronViewerQt(QWidget* parent,
                           const Polyhedron& apoly,
                           const char* title="Basic Polyhedron Viewer",
                           bool anofaces=false,
                           const ColorFunctor& fcolor=ColorFunctor()) :
    // First draw: no vertex; edges, faces; mono-color; inverse normal
    Base(parent, title, false, true, true, true, false),
    poly(apoly),
    m_nofaces(anofaces),
    m_fcolor(fcolor)
  {
    compute_elements();
  }

protected:
  void compute_face(Facet_const_handle fh)
  {
    CGAL::Color c=m_fcolor.run(poly, fh);
    face_begin(c);
    Halfedge_const_handle he=fh->facet_begin();
    do
    {
      add_point_in_face(he->vertex()->point(),
                        get_vertex_normal(he));
      he=he->next();
    }
    while (he!=fh->facet_begin());
    face_end();
  }

  void compute_edge(Halfedge_const_handle he)
  {
    add_segment(he->vertex()->point(),
                he->opposite()->vertex()->point());
    // We can use add_segment(p1, p2, c) with c a CGAL::Color to add a colored segment
  } 

  void compute_vertex(Vertex_const_handle vh)
  {
    add_point(vh->point());
    // We can use add_point(p, c) with c a CGAL::Color to add a colored point
  }

  void compute_elements()
  {
    clear();

    if (!m_nofaces)
    {
      for(typename Polyhedron::Facet_const_iterator f=poly.facets_begin();
          f!=poly.facets_end(); f++)
      {
        if (f!=boost::graph_traits<Polyhedron>::null_face())
        { compute_face(f); }
      }
    }

    for ( typename Polyhedron::Halfedge_const_iterator e=poly.halfedges_begin();
          e!=poly.halfedges_end(); ++e)
    {
      if (e<e->opposite())
      { compute_edge(e); }
    }

    for ( typename Polyhedron::Vertex_const_iterator v=poly.vertices_begin();
          v!=poly.vertices_end(); ++v)
    { compute_vertex(v); }
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
  typename Kernel::Vector_3 get_face_normal(Halfedge_const_handle he)
  {
    typename Kernel::Vector_3 normal=CGAL::NULL_VECTOR;
    Halfedge_const_handle end=he;
    unsigned int nb=0;
    do
    {
      internal::newell_single_step_3(he->vertex()->point(),
                                     he->next()->vertex()->point(),
                                     normal);
      ++nb;
      he=he->next();
    }
    while (he!=end);
    assert(nb>0);
    return (typename Kernel::Construct_scaled_vector_3()(normal, 1.0/nb));
  }
  
  typename Kernel::Vector_3 get_vertex_normal(Halfedge_const_handle he)
  {
    typename Kernel::Vector_3 normal=CGAL::NULL_VECTOR;
    Halfedge_const_handle end=he;
    do
    {
      if (!he->is_border())
      {
        typename Kernel::Vector_3 n=get_face_normal(he);
        normal=typename Kernel::Construct_sum_of_vectors_3()(normal, n);
      }
      he=he->next()->opposite();
    }
    while (he!=end);
    
    if (!typename Kernel::Equal_3()(normal, CGAL::NULL_VECTOR))
    { normal=(typename Kernel::Construct_scaled_vector_3()
              (normal, 1.0/CGAL::sqrt(normal.squared_length()))); }
    
    return normal;
  }

protected:
  const Polyhedron& poly;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
};
  
template<class Polyhedron, class ColorFunctor>
void draw(const Polyhedron& apoly,
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
    const char* argv[2]={"polyhedron_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimplePolyhedronViewerQt<Polyhedron, ColorFunctor>
      mainwindow(app.activeWindow(), apoly, title, nofill, fcolor);
    mainwindow.show();
    app.exec();
  }
}

template<class Polyhedron>
void draw(const Polyhedron& apoly, const char* title, bool nofill)
{
  DefaultColorFunctorPolyhedron c;
  draw(apoly, title, nofill, c);
}

template<class Polyhedron>
void draw(const Polyhedron& apoly, const char* title)
{ draw(apoly, title, false); }

template<class Polyhedron>
void draw(const Polyhedron& apoly)
{ draw(apoly, "Basic Polyhedron Viewer"); }

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_POLYHEDRON_H
