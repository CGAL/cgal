// Copyright (c) 2018 CNRS and LIRIS' Establishments (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_LCC_H
#define CGAL_DRAW_LCC_H

#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Linear_cell_complex_operations.h>
#include <CGAL/Random.h>

namespace CGAL
{

// Default color functor; user can change it to have its own face color
struct DefaultDrawingFunctorLCC
{
  /// @return true iff the volume containing dh is drawn.
  template<typename LCC>
  bool draw_volume(const LCC&,
                   typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the face containing dh is drawn.
  template<typename LCC>
  bool draw_face(const LCC&,
                 typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the edge containing dh is drawn.
  template<typename LCC>
  bool draw_edge(const LCC&,
                 typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the vertex containing dh is drawn.
  template<typename LCC>
  bool draw_vertex(const LCC&,
                   typename LCC::Dart_const_handle) const
  { return true; }

  /// @return true iff the volume containing dh is drawn in wireframe.
  template<typename LCC>
  bool volume_wireframe(const LCC&,
                        typename LCC::Dart_const_handle) const
  { return false; }
  /// @return true iff the face containing dh is drawn in wireframe.
  template<typename LCC>
  bool face_wireframe(const LCC&,
                        typename LCC::Dart_const_handle) const
  { return false; }

  /// @return true iff the volume containing dh is colored.
  template<typename LCC>
  bool colored_volume(const LCC&,
                      typename LCC::Dart_const_handle) const
  { return true; }
  /// @return true iff the face containing dh is colored.
  ///  if we have also colored_volume(alcc, dh), the volume color is
  ///  ignored and only the face color is considered.
  template<typename LCC>
  bool colored_face(const LCC&,
                    typename LCC::Dart_const_handle) const
  { return false; }
  /// @return true iff the edge containing dh is colored.
  template<typename LCC>
  bool colored_edge(const LCC&,
                    typename LCC::Dart_const_handle) const
  { return false; }
  /// @return true iff the vertex containing dh is colored.
  template<typename LCC>
  bool colored_vertex(const LCC&,
                      typename LCC::Dart_const_handle) const
  { return false; }

  /// @return the color of the volume containing dh
  ///  used only if colored_volume(alcc, dh) and !colored_face(alcc, dh)
  template<typename LCC>
  CGAL::Color volume_color(const LCC& alcc,
                           typename LCC::Dart_const_handle dh) const
  {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the face containing dh
  ///  used only if colored_face(alcc, dh)
  template<typename LCC>
  CGAL::Color face_color(const LCC& alcc,
                         typename LCC::Dart_const_handle dh) const
  {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the edge containing dh
  ///  used only if colored_edge(alcc, dh)
  template<typename LCC>
  CGAL::Color edge_color(const LCC& alcc,
                         typename LCC::Dart_const_handle dh) const
  {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
  /// @return the color of the vertex containing dh
  ///  used only if colored_vertex(alcc, dh)
  template<typename LCC>
  CGAL::Color vertex_color(const LCC& alcc,
                           typename LCC::Dart_const_handle dh) const
  {
    CGAL::Random random((unsigned int)(alcc.darts().index(dh)));
    return get_random_color(random);
  }
};

template<class LCC, class Local_kernel, int dim=LCC::ambient_dimension>
struct LCC_geom_utils;

template<class LCC, class Local_kernel>
struct LCC_geom_utils<LCC, Local_kernel, 3>
{
  static typename Local_kernel::Vector_3
  get_vertex_normal(const LCC& lcc, typename LCC::Dart_const_handle dh)
  {
    typename Local_kernel::Vector_3 n = internal::Geom_utils
      <typename LCC::Traits, Local_kernel>::
      get_local_vector(CGAL::compute_normal_of_cell_0<LCC>(lcc,dh));
    n = n/(CGAL::sqrt(n*n));
    return n;
  }
};

template<class LCC, class Local_kernel>
struct LCC_geom_utils<LCC, Local_kernel, 2>
{
  static typename Local_kernel::Vector_3
  get_vertex_normal(const LCC&, typename LCC::Dart_const_handle)
  {
    typename Local_kernel::Vector_3 n=CGAL::NULL_VECTOR;
    return n;
  }
};

// Viewer class for LCC
template<class LCC, class DrawingFunctorLCC>
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
  SimpleLCCViewerQt(QWidget* parent,
                    const LCC* alcc=nullptr,
                    const char* title="Basic LCC Viewer",
                    bool anofaces=false,
                    const DrawingFunctorLCC& drawing_functor=DrawingFunctorLCC()) :
    // First draw: vertices; edges, faces; multi-color; inverse normal
    Base(parent, title, true, true, true, false, false),
    lcc(alcc),
    m_oriented_mark(lcc->get_new_mark()),
    m_nofaces(anofaces),
    m_random_face_color(false),
    m_drawing_functor(drawing_functor)
  {
    lcc->orient(m_oriented_mark);
    compute_elements();
  }

  ~SimpleLCCViewerQt()
  { lcc->free_mark(m_oriented_mark); }

protected:
  void set_lcc(const LCC* alcc, bool doredraw=true)
  {
    if (lcc!=nullptr)
    { lcc->free_mark(m_oriented_mark); }

    lcc=alcc;
    m_oriented_mark=lcc->get_new_mark();
    lcc->orient(m_oriented_mark);

    compute_elements();
    if (doredraw) { redraw(); }
  }

  void compute_face(Dart_const_handle dh, Dart_const_handle voldh)
  {
    if (m_nofaces || !m_drawing_functor.draw_face(*lcc, dh)) return;

    // We fill only closed faces.
    Dart_const_handle cur=dh;
    Dart_const_handle min=dh;
    do
    {
      if (!lcc->is_next_exist(cur)) return; // open face=>not filled
      if (cur<min) min=cur;
      cur=lcc->next(cur);
    }
    while(cur!=dh);

    if (m_random_face_color)
    {
      CGAL::Random random((unsigned int)(lcc->darts().index(dh)));
      CGAL::Color c=get_random_color(random);
      face_begin(c);
    }
    else if (m_drawing_functor.colored_face(*lcc, dh))
    {
      CGAL::Color c=m_drawing_functor.face_color(*lcc, dh);
      face_begin(c);
    }
    else if (m_drawing_functor.colored_volume(*lcc, voldh))
    {
      CGAL::Color c=m_drawing_functor.volume_color(*lcc, voldh);
      face_begin(c);
    }
    else
    { face_begin(); }

    cur=dh;
    do
    {
      add_point_in_face(lcc->point(cur), LCC_geom_utils<LCC, Local_kernel>::
                        get_vertex_normal(*lcc, cur));
      cur=lcc->next(cur);
    }
    while(cur!=dh);

    face_end();
  }

  void compute_edge(Dart_const_handle dh)
  {
    if (!m_drawing_functor.draw_edge(*lcc, dh)) return;

    Point p1 = lcc->point(dh);
    Dart_const_handle d2 = lcc->other_extremity(dh);
    if (d2!=nullptr)
    {
      if (m_drawing_functor.colored_edge(*lcc, dh))
      { add_segment(p1, lcc->point(d2), m_drawing_functor.edge_color(*lcc, dh)); }
      else
      { add_segment(p1, lcc->point(d2)); }
    }
  }

  void compute_vertex(Dart_const_handle dh)
  {
    if (!m_drawing_functor.draw_vertex(*lcc, dh)) return;

    if (m_drawing_functor.colored_vertex(*lcc, dh))
    { add_point(lcc->point(dh), m_drawing_functor.vertex_color(*lcc, dh)); }
    else
    { add_point(lcc->point(dh)); }
  }

  void compute_elements()
  {
    clear();
    if (lcc==nullptr) return;

    typename LCC::size_type markvolumes  = lcc->get_new_mark();
    typename LCC::size_type markfaces    = lcc->get_new_mark();
    typename LCC::size_type markedges    = lcc->get_new_mark();
    typename LCC::size_type markvertices = lcc->get_new_mark();

    for (typename LCC::Dart_range::const_iterator it=lcc->darts().begin(),
         itend=lcc->darts().end(); it!=itend; ++it )
    {
      if (!lcc->is_marked(it, markvolumes) &&
          m_drawing_functor.draw_volume(*lcc, it))
      {
        for (typename LCC::template Dart_of_cell_basic_range<3>::
             const_iterator itv=lcc->template darts_of_cell_basic<3>(it, markvolumes).begin(),
             itvend=lcc->template darts_of_cell_basic<3>(it, markvolumes).end();
             itv!=itvend; ++itv)
        {
          lcc->mark(itv, markvolumes); // To be sure that all darts of the basic iterator will be marked
          if (!lcc->is_marked(itv, markfaces) &&
              lcc->is_marked(itv, m_oriented_mark) &&
              m_drawing_functor.draw_face(*lcc, itv))
          {
            if (!m_drawing_functor.volume_wireframe(*lcc, itv) &&
                !m_drawing_functor.face_wireframe(*lcc, itv))
            { compute_face(itv, it); }
            for (typename LCC::template Dart_of_cell_basic_range<2>::
                 const_iterator itf=lcc->template darts_of_cell_basic<2>(itv, markfaces).begin(),
                 itfend=lcc->template darts_of_cell_basic<2>(itv, markfaces).end();
                 itf!=itfend; ++itf)
            {
              if (!m_drawing_functor.volume_wireframe(*lcc, itv) &&
                !m_drawing_functor.face_wireframe(*lcc, itv))
              { lcc->mark(itf, markfaces); } // To be sure that all darts of the basic iterator will be marked
              if ( !lcc->is_marked(itf, markedges)  &&
                   m_drawing_functor.draw_edge(*lcc, itf))
              {
                compute_edge(itf);
                for (typename LCC::template Dart_of_cell_basic_range<1>::
                     const_iterator ite=lcc->template darts_of_cell_basic<1>(itf, markedges).begin(),
                     iteend=lcc->template darts_of_cell_basic<1>(itf, markedges).end();
                     ite!=iteend; ++ite)
                {
                  lcc->mark(ite, markedges); // To be sure that all darts of the basic iterator will be marked
                  if ( !lcc->is_marked(ite, markvertices)  &&
                       m_drawing_functor.draw_vertex(*lcc, ite))
                  {
                    compute_vertex(ite);
                    CGAL::mark_cell<LCC, 0>(*lcc, ite, markvertices);
                  }
                }
              }
            }
          }
        }
      }
    }

    for (typename LCC::Dart_range::const_iterator it=lcc->darts().begin(),
         itend=lcc->darts().end(); it!=itend; ++it )
    {
      lcc->unmark(it, markvertices);
      lcc->unmark(it, markedges);
      lcc->unmark(it, markfaces);
      lcc->unmark(it, markvolumes);

    }

    lcc->free_mark(markvolumes);
    lcc->free_mark(markfaces);
    lcc->free_mark(markedges);
    lcc->free_mark(markvertices);
  }

  virtual void init()
  {
    Base::init();
    setKeyDescription(::Qt::Key_R, "Toggles random face colors");
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    if ((e->key()==::Qt::Key_R) && (modifiers==::Qt::NoButton))
    {
      m_random_face_color=!m_random_face_color;
      displayMessage(QString("Random face color=%1.").arg(m_random_face_color?"true":"false"));
      compute_elements();
      redraw();
    }
    else
    { Base::keyPressEvent(e); } // Call the base method to process others/classicals key

    // Call: * compute_elements() if the model changed, followed by
    //       * redraw() if some viewing parameters changed that implies some
    //                  modifications of the buffers
    //                  (eg. type of normal, color/mono)
    //       * update() just to update the drawing
  }

protected:
  const LCC* lcc;
  typename LCC::size_type m_oriented_mark;
  bool m_nofaces;
  bool m_random_face_color;
  const DrawingFunctorLCC& m_drawing_functor;
};

// Specialization of draw function.
#define CGAL_LCC_TYPE CGAL::Linear_cell_complex_base                    \
  <d_, ambient_dim, Traits_, Items_, Alloc_, Map, Refs, Storage_>

template < unsigned int d_, unsigned int ambient_dim,
           class Traits_,
           class Items_,
           class Alloc_,
           template<unsigned int,class,class,class,class>
           class Map,
           class Refs,
           class Storage_,
           class DrawingFunctorLCC=DefaultDrawingFunctorLCC>
void draw(const CGAL_LCC_TYPE& alcc,
          const char* title="LCC for CMap Basic Viewer",
          bool nofill=false,
          const DrawingFunctorLCC& drawing_functor=DrawingFunctorLCC())
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
    SimpleLCCViewerQt<CGAL_LCC_TYPE, DrawingFunctorLCC>
      mainwindow(app.activeWindow(), &alcc, title, nofill, drawing_functor);
    mainwindow.show();
    app.exec();
  }
}

// Todo a function taking a const DrawingFunctorLCC& drawing_functor as parameter
#undef CGAL_LCC_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_LCC_H
