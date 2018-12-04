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

#ifndef CGAL_LCC_WITH_PATHS_H
#define CGAL_LCC_WITH_PATHS_H

#include <CGAL/Basic_viewer_qt.h>
#include <CGAL/Random.h>

namespace CGAL
{
  
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
class LCC_with_path_viewer : public Basic_viewer_qt
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
  LCC_with_path_viewer(const LCC& alcc,
                       const std::vector<const Path_on_surface<LCC>*>& paths,
                       std::size_t amark=LCC::INVALID_MARK,
                       const char* title="", bool anofaces=false,
                       const ColorFunctor& fcolor=ColorFunctor()) :
    Base(title),
    lcc(alcc),
    m_nofaces(anofaces),
    m_fcolor(fcolor),
    m_paths(paths),
    m_current_path(m_paths.size()),
    m_current_dart(alcc.number_of_darts()),
    m_draw_marked_darts(true),
    m_amark(amark)
  { compute_elements(); }

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
        CGAL::mark_cell<LCC, 1>(lcc, selected_dart, markedges);
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
          CGAL::mark_cell<LCC, 1>(lcc, it, markedges);
        }
      }
    }
    else
    {
      if (m_current_path==m_paths.size())
      {
        for (unsigned int i=0; i<m_paths.size(); ++i)
        { compute_path(i, markedges); }
      }
      else if (m_current_path!=m_paths.size()+1)
      { compute_path(m_current_path, markedges); }

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

        /*if ( !lcc.is_marked(it, markvertices) )
        {
          compute_vertex(it);
          CGAL::mark_cell<LCC, 0>(lcc, it, markvertices);
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
    
    // CGAL::Color c=m_fcolor.run(lcc, dh);
    face_begin(); //c);

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
      initialize_buffers();
      compile_shaders();
      updateGL();
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
      initialize_buffers();
      compile_shaders();
      updateGL();
    }
    else if ((e->key()==::Qt::Key_N) && (modifiers==::Qt::NoButton))
    {
      m_current_path=(m_current_path+1)%(m_paths.size()+2);
      if (m_current_path==m_paths.size())
      { displayMessage(QString("Draw all paths.")); }
      else if (m_current_path==m_paths.size()+1)
      { displayMessage(QString("Do not draw paths.")); }
      else
      { displayMessage(QString("Draw path=%1.").arg((m_current_path))); }
      compute_elements();
      initialize_buffers();
      compile_shaders();
      updateGL();
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
      initialize_buffers();
      compile_shaders();
      updateGL();
    }
    else
    { Base::keyPressEvent(e); }
  }

  void compute_path(unsigned int i, unsigned int amark)
  {
    if (m_paths[i]->is_empty())
    { return; }

    CGAL::Random random(i);
    CGAL::Color color=get_random_color(random);
    
    add_point(lcc.point(m_paths[i]->get_ith_dart(0)), color);
    for (unsigned int j=0; j<m_paths[i]->length(); ++j)
    {
      if ( !lcc.is_marked( m_paths[i]->get_ith_dart(j), amark) )
      {
        compute_edge(m_paths[i]->get_ith_dart(j), color);
        CGAL::mark_cell<LCC, 1>(lcc, m_paths[i]->get_ith_dart(j), amark);
      }
    }
  }
  
protected:
  const LCC& lcc;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
  const std::vector<const Path_on_surface<LCC>*>& m_paths;
  unsigned int m_current_path;
  unsigned int m_current_dart;
  bool m_draw_marked_darts;
  std::size_t m_amark; // If !=INVALID_MARK, show darts marked with this mark
};
  
template<class LCC, class ColorFunctor>
void display(const LCC& alcc,
             std::vector<const Path_on_surface<LCC>*> paths,
             const char* title="LCC Viewer",
             std::size_t amark=LCC::INVALID_MARK,
             bool nofill=false,
             const ColorFunctor& fcolor=ColorFunctor())
{
  int argc=1;

  const char* argv[2]={"lccviewer","\0"};
  QApplication app(argc,const_cast<char**>(argv));

  LCC_with_path_viewer<LCC, ColorFunctor> mainwindow(alcc, paths, amark,
                                                     title, nofill, fcolor);
  mainwindow.show();

  app.exec();
}

template<class LCC>
void display(const LCC& alcc,
             std::vector<const Path_on_surface<LCC>*> paths,
             const char* title="LCC Viewer",
             std::size_t amark=LCC::INVALID_MARK,
             bool nofill=false)
{ return display<LCC, DefaultColorFunctor>(alcc, paths, title, amark, nofill); }

} // End namespace CGAL

#endif // CGAL_LCC_WITH_PATHS_H
