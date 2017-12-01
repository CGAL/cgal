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

#ifndef CGAL_T3_VIEWER_QT_H
#define CGAL_T3_VIEWER_QT_H

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Random.h>

namespace CGAL
{
  
// Default color functor; user can change it to have its own face color
struct DefaultColorFunctor
{
  template<typename T3>
  static CGAL::Color run(const T3& at3,
                         const typename T3::Finite_facets_iterator* fh)
  {
    if (fh==NULL) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    // Here dh is the smaller dart of its face
    CGAL::Random random((unsigned long int)(&*((*fh)->first))+(*fh)->second);
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
  SimpleTriangulation3ViewerQt(const T3& at3, const char* title="",
                               bool anofaces=false,
                               const ColorFunctor& fcolor=ColorFunctor()) :
    Base(title),
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

    for (typename T3::Finite_facets_iterator it=t3.finite_facets_begin();
         it!=t3.finite_facets_end(); ++it)
    { compute_face(it); } 

    for (typename T3::Finite_edges_iterator it=t3.finite_edges_begin();
         it!=t3.finite_edges_end(); ++it)
    { compute_edge(it); } 

    for (typename T3::Finite_vertices_iterator it=t3.finite_vertices_begin();
         it!=t3.finite_vertices_end(); ++it)
    { compute_vertex(it); } 
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    // const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    Base::keyPressEvent(e);
  }

protected:
  const T3& t3;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
};
  
template<class T3, class ColorFunctor>
void display(const T3& at3,
             const char* title="T3 Viewer",
             bool nofill=false,
             const ColorFunctor& fcolor=ColorFunctor())
{
  int argc=1;

  const char* argv[2]={"t3_viewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimpleTriangulation3ViewerQt<T3, ColorFunctor> mainwindow(at3, title,
                                                            nofill, fcolor);
  mainwindow.show();

  app.exec();
}

template<class T3>
void display(const T3& at3,
             const char* title="t3_viewer",
             bool nofill=false)
{ return display<T3, DefaultColorFunctor>(at3, title, nofill); }

} // End namespace CGAL

#endif // CGAL_T3_VIEWER_QT_H
