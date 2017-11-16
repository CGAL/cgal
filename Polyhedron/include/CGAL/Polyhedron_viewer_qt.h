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

#ifndef CGAL_POLYHEDRON_VIEWER_QT_H
#define CGAL_POLYHEDRON_VIEWER_QT_H

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Random.h>

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctor
{
  template<typename Polyhedron>
  static CGAL::Color run(const Polyhedron& alcc,
                         typename Polyhedron::Facet_const_handle fh)
  {
    if (fh==boost::graph_traits<Polyhedron>::null_face()) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    CGAL::Random random((unsigned long int)&(*fh));
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

// Viewer class for Polyhedron
template<class Polyhedron, class ColorFunctor>
class SimplePolyhedronViewerQt : public Basic_viewer
{
  typedef Basic_viewer Base;
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
  SimplePolyhedronViewerQt(const Polyhedron& apoly, const char* title="", bool anofaces=false) :
    Base(title),
    poly(apoly),
    m_nofaces(anofaces)
  {
    compute_elements();
  }

protected:
  void compute_face(Facet_const_handle fh)
  {
    CGAL::Color c=ColorFunctor::run(poly, fh);
    colored_face_begin(c);
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
    add_mono_segment(he->vertex()->point(),
                     he->opposite()->vertex()->point());
  } 

  void compute_vertex(Vertex_const_handle vh)
  { add_mono_point(vh->point()); }

  void compute_elements()
  {
    clear();

    if (!m_nofaces)
    {
      for(typename Polyhedron::Facet_const_iterator f=poly.facets_begin(); f!=poly.facets_end(); f++)
      {
        if (f!=boost::graph_traits<Polyhedron>::null_face())
        { compute_face(f); }
      }
    }

    for ( typename Polyhedron::Halfedge_const_iterator e=poly.halfedges_begin(); e!=poly.halfedges_end(); ++e)
    {
      if (e<e->opposite())
      { compute_edge(e); }
    }

    for ( typename Polyhedron::Vertex_const_iterator v=poly.vertices_begin(); v!=poly.vertices_end(); ++v)
    { compute_vertex(v); }
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();
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
      internal::newell_single_step_3(he->vertex()->point(), he->next()->vertex()->point(), normal);
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
    { normal=(typename Kernel::Construct_scaled_vector_3()(normal,
                                                           1.0/CGAL::sqrt(normal.squared_length()))); }
    
    return normal;
  }

protected:
  const Polyhedron& poly;
  bool m_nofaces;
};

  
template<class Polyhedron, class ColorFunctor>
void display(const Polyhedron& apoly,
             const char* title="",
             bool nofill=false)
{
  int argc=1;

  const char* argv[2]={"polyhedron viewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimplePolyhedronViewerQt<Polyhedron, ColorFunctor> mainwindow(apoly, title, nofill);
  mainwindow.show();

  app.exec();
}

template<class Polyhedron>
void display(const Polyhedron& apoly,
             const char* title="",
             bool nofill=false)
{ return display<Polyhedron, DefaultColorFunctor>(apoly, title, nofill); }

#endif // CGAL_POLYHEDRON_VIEWER_QT_H
