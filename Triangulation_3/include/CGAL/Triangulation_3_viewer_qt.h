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

#include "CGAL/Qt/Basic_viewer_qt.h"
#include <CGAL/Random.h>

CGAL::Color get_random_color(CGAL::Random& random)
{
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

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctor
{
  template<typename T3>
  static CGAL::Color run(const T3& at3,
                         typename LCC::Dart_const_handle dh)
  {
    if (dh==alcc.null_handle) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    // Here dh is the smaller dart of its face
    CGAL::Random random(alcc.darts().index(dh));
    return get_random_color(random);
  }
};

// Viewer class for T3 
template<class T3, class ColorFunctor>
class SimpleTriangulation3ViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename LCC::Traits Kernel;
  typedef typename T3::Cell_handle    Cell_handle;
  typedef typename T3::Point          Point;
  
public:
  /// Construct the viewer.
  /// @param at3 the t3 to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleLCCViewerQt(const T3& at3, const char* title="", bool anofaces=false) :
    Base(title),
    t3(at3),
    m_nofaces(anofaces)
  {
    compute_elements();
  }

protected:
  void compute_face(Facet_const_handle fh)
  {
    /*    // We fill only closed faces.
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
      face_begin(); // TODO REMOVE LATER
    else
    {
      // c=CGAL::Color(100,255,100);
      face_begin(c);
    }

    cur=dh;
    do
    {
      add_point_in_face(lcc.point(cur),
                        Geom_utils<LCC>::get_vertex_normal(lcc, cur));
      cur=lcc.next(cur);
    }
    while(cur!=dh);
    */
    face_end();
  }

  void compute_edge(Edge_const_handle eh)
  {
    /*    CGAL::Random random((unsigned long int)&*dh);
    CGAL::Color c=get_random_color(random); // TODO REMOVE LATER

    Point p1 = lcc.point(dh);
    Dart_const_handle d2 = lcc.other_extremity(dh);
    if (d2!=NULL)
    {
      Point p2 = lcc.point(d2);
      if (c.red()<60 || c.green()<60 || c.blue()<60)
        add_segment(p1, p2);
      else
      {
        add_segment(p1, p2, c); // TODO REMOVE LATER
      }
      }*/
  }

  void compute_vertex(Vertex_const_handle vh)
  {
    /*    CGAL::Random random((unsigned long int)&*dh);  // TODO REMOVE LATER
    CGAL::Color c=get_random_color(random);

    Point p = lcc.point(dh);
    if (c.red()<60 || c.green()<60 || c.blue()<60)
      add_point(p);
    else
    {
      add_point(p, c); // TODO REMOVE LATER
      }*/
  }

  void compute_elements()
  {
    clear();

    for (typename T3::Finite_facets_iterator it=t3.finite_facets_begin(); it!=t3.finite_facets_end(); ++it)
    { compute_face(it); } 

    for (typename T3::Finite_edges_iterator it=t3.finite_edges_begin(); it!=t3.finite_edges_end(); ++it)
    { compute_edge(it); } 

    for (typename T3::Finite_vertices_iterator it=t3.finitevertices_begin(); it!=t3.finite_vertices_end(); ++it)
    { compute_vertex(it); } 
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    Base::keyPressEvent(e);
  }

protected:
  const T3& t3;
  bool m_nofaces;
};

  
template<class T3, class ColorFunctor>
void display(const T3& at3,
             const char* title="",
             bool nofill=false)
{
  int argc=1;

  const char* argv[2]={"T3 viewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimpleTriangulation3ViewerQt<T3, ColorFunctor> mainwindow(at3, title, nofill);
  mainwindow.show();

  app.exec();
}

template<class T3>
void display(const T3& at3,
             const char* title="",
             bool nofill=false)
{ return display<T3, DefaultColorFunctor>(at3, title, nofill); }

#endif // CGAL_T3_VIEWER_QT_H
