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

#ifndef CGAL_DRAW_LCC_H
#define CGAL_DRAW_LCC_H

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Random.h>

namespace CGAL
{
  
// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorLCC
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
  SimpleLCCViewerQt(const LCC& alcc, const char* title="", bool anofaces=false,
                    const ColorFunctor& fcolor=ColorFunctor()) :
    Base(title),
    lcc(alcc),
    m_nofaces(anofaces),
    m_fcolor(fcolor)
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
    
    CGAL::Color c=m_fcolor.run(lcc, dh);
    face_begin(c);

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
    Point p1 = lcc.point(dh);
    Dart_const_handle d2 = lcc.other_extremity(dh);
    if (d2!=NULL)
    { add_segment(p1, lcc.point(d2)); }
  }

  void compute_vertex(Dart_const_handle dh)
  { add_point(lcc.point(dh)); }

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
  const LCC& lcc;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
};
  
template<class LCC, class ColorFunctor>
void draw(const LCC& alcc,
          const char* title="LCC Viewer",
          bool nofill=false,
          const ColorFunctor& fcolor=ColorFunctor())
{
  int argc=1;

  const char* argv[2]={"lccviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  SimpleLCCViewerQt<LCC, ColorFunctor> mainwindow(alcc, title, nofill, fcolor);

#if !defined(CGAL_TEST_SUITE)
  mainwindow.show();
  app.exec();
#endif
}

template<class LCC>
void draw(const LCC& alcc,
          const char* title="LCC Viewer",
          bool nofill=false)
{ return draw<LCC, DefaultColorFunctorLCC>(alcc, title, nofill); }

} // End namespace CGAL

#else // CGAL_USE_BASIC_VIEWER

namespace CGAL 
{
  
template<class LCC, class ColorFunctor>
void draw(const LCC&,
          const char* ="LCC Viewer",
          bool=false,
          const ColorFunctor& =ColorFunctor())
{
  std::cerr<<"Impossible to draw a Linear_cell_complex because CGAL_USE_BASIC_VIEWER is not defined."
           <<std::endl;
}

template<class LCC>
void draw(const LCC&,
          const char* ="LCC Viewer",
          bool=false)
{
  std::cerr<<"Impossible to draw a Linear_cell_complex because CGAL_USE_BASIC_VIEWER is not defined."
           <<std::endl;
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_LCC_H
