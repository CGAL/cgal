// Copyright (c) 2017 CNRS and LIRIS' Establishments (France).
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
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_SURFACE_MESH_VIEWER_QT_H
#define CGAL_SURFACE_MESH_VIEWER_QT_H

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Random.h>

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctor
{
  template<typename SM>
  static CGAL::Color run(const SM& amesh,
                         typename SM::Face_index fh)
  {
    if (fh==boost::graph_traits<SM>::null_face()) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    CGAL::Random random(fh);
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
  SimpleSurfaceMeshViewerQt(const SM& amesh, const char* title="", bool anofaces=false) :
    Base(title),
    sm(amesh),
    m_nofaces(anofaces)
  {
    compute_elements();
  }

protected:
  void compute_face(face_descriptor fh)
  {
    CGAL::Color c=ColorFunctor::run(sm, fh);
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
      for (typename SM::Face_range::iterator f=sm.faces().begin(); f!=sm.faces().end(); ++f)
      {
        if (*f!=boost::graph_traits<SM>::null_face())
        { compute_face(*f); }
      }
    }
    
    for (typename SM::Edge_range::iterator e=sm.edges().begin(); e!=sm.edges().end(); ++e)
    { compute_edge(*e); }

    for (typename SM::Vertex_range::iterator v=sm.vertices().begin(); v!=sm.vertices().end(); ++v)
    { compute_vertex(*v); }
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();
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
    { normal=(typename Kernel::Construct_scaled_vector_3()(normal,
                                                           1.0/CGAL::sqrt(normal.squared_length()))); }
    
    return normal;
  }

protected:
  const SM& sm;
  bool m_nofaces;
};

  
template<class SM, class ColorFunctor>
void display(const SM& amesh,
             const char* title="",
             bool nofill=false)
{
  int argc=1;

  const char* argv[2]={"Surface mesh viewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimpleSurfaceMeshViewerQt<SM, ColorFunctor> mainwindow(amesh, title, nofill);
  mainwindow.show();

  app.exec();
}

template<class SM>
void display(const SM& amesh,
             const char* title="",
             bool nofill=false)
{ return display<SM, DefaultColorFunctor>(amesh, title, nofill); }

#endif // CGAL_SURFACE_MESH_VIEWER_QT_H
