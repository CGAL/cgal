// Copyright (c) 2019 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//
#ifndef CGAL_DRAW_FACE_GRAPH_WITH_PATHS_H
#define CGAL_DRAW_FACE_GRAPH_WITH_PATHS_H

#include <CGAL/license/Surface_mesh_topology.h>

#include <iostream>
#include <initializer_list>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/Path_on_surface.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Random.h>

namespace CGAL {

// Specialisation for face graph; otherwise use the LCC_geom_utils of LCC.
template<class Mesh, class Local_kernel>
struct LCC_geom_utils<CGAL::Face_graph_wrapper<Mesh>, Local_kernel, 3>
{
  static typename Get_traits<Mesh>::Vector
  get_face_normal(const CGAL::Face_graph_wrapper<Mesh>& mesh,
                  typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_handle dh)
  {
    typename Get_traits<Mesh>::Vector normal(CGAL::NULL_VECTOR);
    const typename Get_traits<Mesh>::Point*
        curr=&Get_traits<Mesh>::get_point(mesh.get_fg(), dh);
    typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_handle adart=dh;
    unsigned int nb=0;

    do
    {
      const typename Get_traits<Mesh>::Point*
          next=&Get_traits<Mesh>::get_point(mesh.get_fg(),
                                            mesh.other_extremity(adart));
      internal::newell_single_step_3_for_lcc(*curr, *next, normal);
      ++nb;
      curr=next;
      adart=mesh.next(adart);
    }
    while(adart!=dh);

    assert(nb>0);
    return typename Get_traits<Mesh>::Kernel::Construct_scaled_vector_3()
      (normal, 1.0/nb);
  }
  static typename Local_kernel::Vector_3
  get_vertex_normal(const CGAL::Face_graph_wrapper<Mesh>& mesh,
                    typename CGAL::Face_graph_wrapper<Mesh>::Dart_const_handle dh)
  {
    typename Get_traits<Mesh>::Vector normal(CGAL::NULL_VECTOR);
    unsigned int nb = 0;

    for ( typename CGAL::Face_graph_wrapper<Mesh>::template Dart_of_cell_range<0>::
          const_iterator it=mesh.template darts_of_cell<0>(dh).begin(),
          itend=mesh.template darts_of_cell<0>(dh).end(); it!=itend; ++it )
    {
      normal=typename Get_traits<Mesh>::Kernel::Construct_sum_of_vectors_3()
          (normal, get_face_normal(mesh, it));
      ++nb;
    }

    if ( nb<2 ) return internal::Geom_utils
                  <typename Get_traits<Mesh>::Kernel, Local_kernel>::
                  get_local_vector(normal);

    return internal::Geom_utils
      <typename Get_traits<Mesh>::Kernel, Local_kernel>::
      get_local_vector(typename Get_traits<Mesh>::Kernel::
                       Construct_scaled_vector_3()(normal, 1.0/nb));
  }
};

// Viewer class for Face_graph with paths.
template<class Mesh, class DrawingFunctorLCC>
class Face_graph_with_path_viewer : public Basic_viewer_qt
{
  typedef Basic_viewer_qt                         Base;
  typedef typename Get_map<Mesh, Mesh>::type      LCC;
  typedef typename LCC::Dart_const_handle         Dart_const_handle;
  typedef typename LCC::size_type                 size_type;
  typedef typename CGAL::Get_traits<Mesh>::Kernel Kernel;
  typedef typename CGAL::Get_traits<Mesh>::Point  Point;
  typedef typename CGAL::Get_traits<Mesh>::Vector Vector;

public:
  /// Construct the viewer.
  /// @param alcc the lcc to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed;
  ///     this can be usefull for very big object where this time could be long)
  Face_graph_with_path_viewer(QWidget* parent,
                              const Mesh& amesh,
                              const std::vector
                              <Surface_mesh_topology::Path_on_surface<Mesh> >
                              *paths=nullptr,
                              size_type amark=LCC::INVALID_MARK,
                              const char* title="", bool anofaces=false,
                              const DrawingFunctorLCC&
                              drawing_functor=DrawingFunctorLCC()) :
    Base(parent, title, true, true, true, false, true),
    mesh(amesh),
    lcc(amesh),
    m_oriented_mark(lcc.get_new_mark()),
    m_nofaces(anofaces),
    m_drawing_functor(drawing_functor),
    m_paths(paths),
    m_current_path(m_paths->size()),
    m_current_dart(lcc.darts().end()),
    m_draw_marked_darts(true),
    m_amark(amark==(std::numeric_limits<size_type>::max)()?
            LCC::INVALID_MARK:amark)
  {
    lcc.orient(m_oriented_mark);
    compute_elements();
  }

  ~Face_graph_with_path_viewer()
  {
    lcc.free_mark(m_oriented_mark);
  }

protected:

  const Point& get_point(Dart_const_handle dh) const
  { return CGAL::Get_traits<Mesh>::get_point(mesh, dh); }

  void compute_elements()
  {
    clear();

    typename LCC::size_type markfaces    = lcc.get_new_mark();
    typename LCC::size_type markedges    = lcc.get_new_mark();
    typename LCC::size_type markvertices = lcc.get_new_mark();

    if (m_current_dart!=lcc.darts().end())
    { // We want to draw only one dart
      Dart_const_handle selected_dart=m_current_dart; //lcc.dart_handle(m_current_dart);
      compute_edge(selected_dart, CGAL::Color(255,0,0));
      lcc.template mark_cell<1>(selected_dart, markedges);
      compute_vertex(selected_dart);

      if ( !m_nofaces )
      { compute_face(selected_dart); }

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
        for (std::size_t i=0; i<m_paths->size(); ++i)
        { compute_path(i, markedges); }
      }
      else if (m_current_path!=m_paths->size()+1)
      { compute_path(m_current_path, markedges); }

      for (typename LCC::Dart_range::const_iterator it=lcc.darts().begin(),
           itend=lcc.darts().end(); it!=itend; ++it )
      {
        if (!m_nofaces && !lcc.is_marked(it, markfaces) &&
            !lcc.is_perforated(it) && lcc.is_marked(it, m_oriented_mark))
        {
          compute_face(it);
          lcc.template mark_cell<2>(it, markfaces);
        }

        if ( !lcc.is_marked(it, markedges) )
        {
          compute_edge(it);
          lcc.template mark_cell<1>(it, markedges);
        }

        if ( !lcc.is_marked(it, markvertices) )
        {
          compute_vertex(it);
          lcc.template mark_cell<0>(it, markvertices);
        }
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
      add_point_in_face(get_point(cur),
                        LCC_geom_utils<LCC, Local_kernel>::
                        get_vertex_normal(lcc, cur));
      cur=lcc.next(cur);
    }
    while(cur!=dh);

    face_end();
  }

  void compute_edge(Dart_const_handle dh)
  {
    Point p1 = get_point(dh);
    Dart_const_handle d2 = lcc.other_extremity(dh);
    if (d2!=LCC::null_handle)
    {
      if (m_draw_marked_darts && m_amark!=LCC::INVALID_MARK &&
          (lcc.is_marked(dh, m_amark) || lcc.is_marked(lcc.opposite2(dh), m_amark)))
      { add_segment(p1, get_point(d2), CGAL::Color(0, 0, 255)); }
      else
      { add_segment(p1, get_point(d2)); }
    }
  }

  void compute_edge(Dart_const_handle dh, const CGAL::Color& color)
  {
    Point p1 = get_point(dh);
    Dart_const_handle d2 = lcc.other_extremity(dh);
    if (d2!=LCC::null_handle)
    { add_segment(p1, get_point(d2), color); }
  }

  void compute_vertex(Dart_const_handle dh)
  { add_point(get_point(dh)); }

  virtual void init()
  {
    Base::init();
    setKeyDescription(::Qt::Key_D, "Increase current dart drawing");
    setKeyDescription(::Qt::Key_D+::Qt::ControlModifier, "Decrease current dart drawing");
    setKeyDescription(::Qt::Key_D+::Qt::ShiftModifier, "Draw all darts");

    setKeyDescription(::Qt::Key_X, "Toggles marked darts display");

    setKeyDescription(::Qt::Key_P, "Increase current path drawing");
    setKeyDescription(::Qt::Key_P+::Qt::ControlModifier, "Decrease current path drawing");
    setKeyDescription(::Qt::Key_P+::Qt::ShiftModifier, "Draw all paths");
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    const ::Qt::KeyboardModifiers modifiers = e->modifiers();

    if ((e->key()==::Qt::Key_D) && (modifiers==::Qt::NoButton))
    {
      if (m_current_dart==lcc.darts().end())
      { m_current_dart=lcc.darts().begin(); }
      else
      { ++m_current_dart; }
      if (m_current_dart==lcc.darts().end())
      { displayMessage(QString("Draw all darts.")); }
      else
      { displayMessage(QString("Draw dart=%1.").arg(lcc.darts().index(m_current_dart))); }
      compute_elements();
      redraw();
    }
    else if ((e->key()==::Qt::Key_D) && (modifiers==::Qt::ControlModifier))
    {
      if (m_current_dart==lcc.darts().begin())
      { m_current_dart=lcc.darts().end(); }
      else
      { --m_current_dart; }
      if (m_current_dart==lcc.darts().end())
      { displayMessage(QString("Draw all darts.")); }
      else
      { displayMessage(QString("Draw dart=%1.").arg(lcc.darts().index(m_current_dart))); }
      compute_elements();
      redraw();
    }
    else if ((e->key()==::Qt::Key_D) && (modifiers==::Qt::ShiftModifier))
    {
      m_current_dart=lcc.darts().end();
      displayMessage(QString("Draw all darts."));
      compute_elements();
      redraw();
    }
    else if ((e->key()==::Qt::Key_P) && (modifiers==::Qt::NoButton))
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
    else if ((e->key()==::Qt::Key_P) && (modifiers==::Qt::ControlModifier))
    {
      m_current_path=(m_current_path==0?m_paths->size()+1:m_current_path-1);
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
    else if ((e->key()==::Qt::Key_P) && (modifiers==::Qt::ShiftModifier))
    {
      m_current_path=m_paths->size();
      displayMessage(QString("Draw all paths."));
      compute_elements();
      redraw();
    }
    else if ((e->key()==::Qt::Key_X) && (modifiers==::Qt::NoButton))
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
    else
    { Base::keyPressEvent(e); }
  }

  void compute_path(std::size_t i, typename LCC::size_type amark)
  {
    if ((*m_paths)[i].is_empty())
    { return; }

    CGAL::Random random(static_cast<unsigned int>(i));
    CGAL::Color color=get_random_color(random);

    add_point(get_point((*m_paths)[i].get_ith_dart(0)), color);
    for (std::size_t j=0; j<(*m_paths)[i].length(); ++j)
    {
      if ( !lcc.is_marked( (*m_paths)[i].get_ith_dart(j), amark) )
      {
        compute_edge((*m_paths)[i].get_ith_dart(j), color);
        lcc.template mark_cell<1>((*m_paths)[i].get_ith_dart(j), amark);
      }
    }
  }

protected:
  const Mesh& mesh;
  const typename Get_map<Mesh, Mesh>::storage_type lcc;
  typename LCC::size_type m_oriented_mark;
  bool m_nofaces;
  const DrawingFunctorLCC& m_drawing_functor;
  const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >* m_paths;
  std::size_t m_current_path;
  typename LCC::Dart_range::const_iterator m_current_dart;
  bool m_draw_marked_darts;
  typename LCC::size_type m_amark; // If !=INVALID_MARK, show darts marked with this mark
};

template<class Mesh, class DrawingFunctor>
void draw(const Mesh& alcc,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& paths,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          bool nofill=false,
          const DrawingFunctor& drawing_functor=DrawingFunctor())
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"lccviewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    Face_graph_with_path_viewer<Mesh, DrawingFunctor> mainwindow(app.activeWindow(),
                                                                 alcc, &paths, amark,
                                                                 title, nofill,
                                                                 drawing_functor);
    mainwindow.show();
    app.exec();
  }
}

template<class Mesh>
void draw(const Mesh& alcc,
          const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& paths,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          bool nofill=false)
{
  DefaultDrawingFunctorLCC f;
  draw<Mesh, DefaultDrawingFunctorLCC>(alcc, paths, title, amark, nofill, f);
}

template<class Mesh>
void draw(const Mesh& alcc,
          std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>> l,
          const char* title="Mesh Viewer With Path",
          typename Get_map<Mesh, Mesh>::type::size_type amark=
          (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
          bool nofill=false)
{
  std::vector<Surface_mesh_topology::Path_on_surface<Mesh> > paths=l;
  draw(alcc, paths, title, amark, nofill);
}

} // End namespace CGAL

#else  // CGAL_USE_BASIC_VIEWER

namespace CGAL
{

  template<class Mesh>
  void draw(const Mesh&,
            const std::vector<Surface_mesh_topology::Path_on_surface<Mesh> >& ,
            const char* ="",
            typename Get_map<Mesh, Mesh>::type::size_type=
            (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
            bool=false)
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }

  template<class Mesh>
  void draw(const Mesh&,
            std::initializer_list<Surface_mesh_topology::Path_on_surface<Mesh>>,
            const char* ="",
            typename Get_map<Mesh, Mesh>::type::size_type=
            (std::numeric_limits<typename Get_map<Mesh, Mesh>::type::size_type>::max)(),
            bool=false)
  {
    std::cerr<<"Impossible to draw, CGAL_USE_BASIC_VIEWER is not defined."<<std::endl;
  }

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_FACE_GRAPH_WITH_PATHS_H
