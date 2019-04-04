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

#ifndef CGAL_DRAW_LCC_WITH_PATHS_H
#define CGAL_DRAW_LCC_WITH_PATHS_H

#include <CGAL/draw_linear_cell_complex.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Random.h>
#include <CGAL/Path_on_surface.h>

namespace CGAL
{
  
// Viewer class for LCC with paths.
template<class Mesh, class DrawingFunctorLCC>
class LCC_with_path_viewer : public Basic_viewer_qt
{
  typedef Basic_viewer_qt Base;
  typedef typename Get_map<Mesh, Mesh>::type LCC;
  typedef typename LCC::Dart_const_handle Dart_const_handle;
  typedef typename CGAL::Kernel_traits<typename LCC::Point_3>::Kernel Kernel;
  typedef typename Kernel::Point_3 Point;
  typedef typename Kernel::Vector_3 Vector;
  
public:
  /// Construct the viewer.
  /// @param alcc the lcc to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  LCC_with_path_viewer(QWidget* parent,
                       const Mesh& amesh,
                       const std::vector<Path_on_surface<LCC> >* paths=NULL,
                       std::size_t amark=LCC::INVALID_MARK,
                       const char* title="", bool anofaces=false,
                       const DrawingFunctorLCC& drawing_functor=DrawingFunctorLCC()) :
    Base(parent, title, true, true, true, false, true),
    lcc(amesh),
    m_nofaces(anofaces),
    m_drawing_functor(drawing_functor),
    m_paths(paths),
    m_current_path(m_paths->size()),
    m_current_dart(0),
    m_draw_marked_darts(true),
    m_amark(amark==-1?LCC::INVALID_MARK:amark)
  {
    m_current_dart=lcc.number_of_darts(); compute_elements();
  }

protected:
  
  void compute_elements()
  {
    clear();

    unsigned int markfaces    = lcc.get_new_mark();
    unsigned int markedges    = lcc.get_new_mark();
    unsigned int markvertices = lcc.get_new_mark();

    if (m_current_dart!=lcc.number_of_darts())
    { // We want to draw only one dart
      Dart_const_handle selected_dart=lcc.darts().iterator_to(lcc.darts()[m_current_dart]);
      if (lcc.is_dart_used(selected_dart))
      {
        compute_edge(selected_dart, CGAL::Color(255,0,0));
        lcc.template mark_cell<1>(selected_dart, markedges);
        compute_vertex(selected_dart);

        if ( !m_nofaces )
        { compute_face(selected_dart); }
      }

      for (typename LCC::Dart_range::const_iterator it=lcc.darts().begin(),
           itend=lcc.darts().end(); it!=itend; ++it )
      {
        if ( !lcc.is_marked(it, markedges) )
        {
          compute_edge(it);
          lcc.template mark_cell<1>(it, markedges);
        }
      }
    }
    else
    {
      if (m_current_path==m_paths->size())
      {
        for (unsigned int i=0; i<m_paths->size(); ++i)
        { compute_path(i, markedges); }
      }
      else if (m_current_path!=m_paths->size()+1)
      { compute_path(m_current_path, markedges); }

      for (typename LCC::Dart_range::const_iterator it=lcc.darts().begin(),
           itend=lcc.darts().end(); it!=itend; ++it )
      {
        if ( !m_nofaces && !lcc.is_marked(it, markfaces) )
        {
          compute_face(it);
          lcc.template mark_cell<2>(it, markfaces);
        }

        if ( !lcc.is_marked(it, markedges) )
        {
          compute_edge(it);
          lcc.template mark_cell<1>(it, markedges);
        }

        /*if ( !lcc.is_marked(it, markvertices) )
        {
          compute_vertex(it);
          lcc.template mark_cell<0>(it, markvertices);
        }*/
      }
    }

    lcc.free_mark(markfaces);
    lcc.free_mark(markedges);
    lcc.free_mark(markvertices);
  }

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
    
    // CGAL::Color c=m_fcolor.run(*lcc, dh);
    face_begin(); //c);

    cur=dh;
    do
    {
      add_point_in_face(lcc.point(cur),
                        LCC_geom_utils<LCC, Local_kernel>::
                        get_vertex_normal(lcc, cur));
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
    {
      if (m_draw_marked_darts && m_amark!=LCC::INVALID_MARK &&
          (lcc.is_marked(dh, m_amark) || lcc.is_marked(lcc.beta(dh, 2), m_amark)))
      { add_segment(p1, lcc.point(d2), CGAL::Color(0, 0, 255)); }
      else
      { add_segment(p1, lcc.point(d2)); }
    }
  }

  void compute_edge(Dart_const_handle dh, const CGAL::Color& color)
  {
    Point p1 = lcc.point(dh);
    Dart_const_handle d2 = lcc.other_extremity(dh);
    if (d2!=NULL)
    { add_segment(p1, lcc.point(d2), color); }
  }

  void compute_vertex(Dart_const_handle dh)
  { add_point(lcc.point(dh)); }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    const ::Qt::KeyboardModifiers modifiers = e->modifiers();

    if ((e->key()==::Qt::Key_D) && (modifiers==::Qt::NoButton))
    {
      m_current_dart=(m_current_dart+1)%(lcc.number_of_darts()+1);
      if (m_current_dart==lcc.number_of_darts())
      {
        displayMessage(QString("Draw all darts."));
      }
      else
      {
        displayMessage(QString("Draw dart=%1.").arg((m_current_dart)));
      }
      compute_elements();
      redraw();
    }
    else if ((e->key()==::Qt::Key_M) && (modifiers==::Qt::NoButton))
    {
      m_draw_marked_darts=!m_draw_marked_darts;

      if (m_draw_marked_darts)
      { displayMessage(QString("Draw marked darts in blue.")); }
      else
      {
        displayMessage(QString("Do not draw marked darts in different color."));
      }
      compute_elements();
      redraw();
    }
    else if ((e->key()==::Qt::Key_N) && (modifiers==::Qt::NoButton))
    {
      m_current_path=(m_current_path+1)%(m_paths->size()+2);
      if (m_current_path==m_paths->size())
      { displayMessage(QString("Draw all paths.")); }
      else if (m_current_path==m_paths->size()+1)
      { displayMessage(QString("Do not draw paths.")); }
      else
      { displayMessage(QString("Draw path=%1, nb_darts=%2.").
                       arg(m_current_path).
                       arg((*m_paths)[m_current_path].length())); }
      compute_elements();
      redraw();
    }
    else if ((e->key()==::Qt::Key_P) && (modifiers==::Qt::NoButton))
    {
      m_current_dart=(m_current_dart==0?lcc.number_of_darts():m_current_dart-1);
      if (m_current_dart==lcc.number_of_darts())
      {
        displayMessage(QString("Draw all darts."));
      }
      else
      {
        displayMessage(QString("Draw dart=%1.").arg((m_current_dart)));
      }
      compute_elements();
      redraw();
    }
    else
    { Base::keyPressEvent(e); }
  }

  void compute_path(unsigned int i, unsigned int amark)
  {
    if ((*m_paths)[i].is_empty())
    { return; }

    CGAL::Random random(i);
    CGAL::Color color=get_random_color(random);
    
    add_point(lcc.point((*m_paths)[i].get_ith_dart(0)), color);
    for (unsigned int j=0; j<(*m_paths)[i].length(); ++j)
    {
      if ( !lcc.is_marked( (*m_paths)[i].get_ith_dart(j), amark) )
      {
        compute_edge((*m_paths)[i].get_ith_dart(j), color);
        lcc.template mark_cell<1>((*m_paths)[i].get_ith_dart(j), amark);
      }
    }
  }
  
protected:
  const typename Get_map<LCC, LCC>::storage_type lcc;
  bool m_nofaces;
  const DrawingFunctorLCC& m_drawing_functor;
  const std::vector<Path_on_surface<LCC> >* m_paths;
  unsigned int m_current_path;
  unsigned int m_current_dart;
  bool m_draw_marked_darts;
  std::size_t m_amark; // If !=INVALID_MARK, show darts marked with this mark
};
  
template<class LCC, class DrawingFunctor>
void draw(const LCC& alcc,
          const std::vector<Path_on_surface<LCC> >& paths,
          const char* title="LCC Viewer",
          std::size_t amark=-1,
          bool nofill=false,
          const DrawingFunctor& drawing_functor=DrawingFunctor())
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=false;
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"lccviewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    LCC_with_path_viewer<LCC, DrawingFunctor> mainwindow(app.activeWindow(),
                                                         alcc, &paths, amark,
                                                         title, nofill,
                                                         drawing_functor);
    mainwindow.show();
    app.exec();
  }
}

template<class LCC>
void draw(const LCC& alcc,
          const std::vector<Path_on_surface<LCC> >& paths,
          const char* title="LCC Viewer",
          std::size_t amark=-1,
          bool nofill=false)
{
  DefaultDrawingFunctorLCC f;
  return draw<LCC, DefaultDrawingFunctorLCC>(alcc, paths, title,
                                             amark, nofill, f);
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_LCC_WITH_PATHS_H
