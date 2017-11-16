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

#ifndef CGAL_LCC_VIEWER_QT_H
#define CGAL_LCC_VIEWER_QT_H

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
  template<typename LCC>
  static CGAL::Color run(const LCC& alcc,
                         typename LCC::Dart_const_handle dh)
  {
    if (dh==alcc.null_handle) // use to get the mono color
      return CGAL::Color(100, 125, 200); // R G B between 0-255

    // Here dh is the smaller dart of its face
    CGAL::Random random(alcc.darts().index(dh));
    return get_random_color(random);
  }
};

template<class LCC, int dim=LCC::ambient_dimension>
struct Geom_utils;

template<class LCC>
struct Geom_utils<LCC, 3>
{
  static typename LCC::Vector get_vertex_normal(const LCC& lcc,
                                                typename LCC::Dart_const_handle dh)
  {
    typename LCC::Vector n = CGAL::compute_normal_of_cell_0<LCC>(lcc,dh);
    n = n/(CGAL::sqrt(n*n));
    return n;
  }
};

template<class LCC>
struct Geom_utils<LCC, 2>
{
  static typename LCC::Vector get_vertex_normal(const LCC&,
                                                typename LCC::Dart_const_handle)
  {
    typename LCC::Vector res=CGAL::NULL_VECTOR;
    return res;
  }
};

// Viewer class for LCC 
template<class LCC, class ColorFunctor>
class SimpleLCCViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename LCC::Dart_const_handle Dart_const_handle;
  typedef typename LCC::Traits Kernel;
  typedef typename Kernel::Point Point;
  typedef typename Kernel::Vector Vector;
  
public:
  /// Construct the viewer.
  /// @param alcc the lcc to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimpleLCCViewerQt(const LCC& alcc, const char* title="", bool anofaces=false) :
    Base(title),
    lcc(alcc),
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

    face_end();
  }

  void compute_edge(Dart_const_handle dh)
  {
    CGAL::Random random((unsigned long int)&*dh);
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
    }
  }

  void compute_vertex(Dart_const_handle dh)
  {
    CGAL::Random random((unsigned long int)&*dh);  // TODO REMOVE LATER
    CGAL::Color c=get_random_color(random);

    Point p = lcc.point(dh);
    if (c.red()<60 || c.green()<60 || c.blue()<60)
      add_point(p);
    else
    {
      add_point(p, c); // TODO REMOVE LATER
    }
  }

  void compute_elements()
  {
    clear();

    unsigned int markfaces    = lcc.get_new_mark();
    unsigned int markedges    = lcc.get_new_mark();
    unsigned int markvertices = lcc.get_new_mark();

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
    const Qt::KeyboardModifiers modifiers = e->modifiers();
    Base::keyPressEvent(e);
  }

protected:
  const LCC& lcc;
  bool m_nofaces;
};

  
template<class LCC, class ColorFunctor>
void display(const LCC& alcc,
             const char* title="",
             bool nofill=false)
{
  int argc=1;

  const char* argv[2]={"lccviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimpleLCCViewerQt<LCC, ColorFunctor> mainwindow(alcc, title, nofill);
  mainwindow.show();

  app.exec();
}

template<class LCC>
void display(const LCC& alcc,
             const char* title="",
             bool nofill=false)
{ return display<LCC, DefaultColorFunctor>(alcc, title, nofill); }

#endif // CGAL_LCC_VIEWER_QT_H
