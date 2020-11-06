// Copyright(c) 2018  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_SS2_H
#define CGAL_DRAW_SS2_H

#include <CGAL/license/Straight_skeleton_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Straight_skeleton_2.h>
#include <CGAL/Random.h>

#include <sstream>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorSS2
{
  template<typename SS2>
  static CGAL::Color run(const SS2&,
                         const typename SS2::Finite_faces_iterator fh)
  {
    CGAL::Random random((unsigned int)(std::size_t)(&*fh));
    return get_random_color(random);
  }
};

// Viewer class for SS2
template<class SS2, class ColorFunctor>
class SimpleStraightSkeleton2ViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt                    Base;
  typedef typename SS2::Vertex_const_handle   Vertex_const_handle;
  typedef typename SS2::Halfedge_const_handle Halfedge_const_handle;

  //  typedef typename SS2::Point                 Point;

public:
  /// Construct the viewer.
  /// @param ass2 the ss2 to view
  /// @param title the title of the window
  SimpleStraightSkeleton2ViewerQt(QWidget* parent, const SS2& ass2,
                               const char* title="Basic SS2 Viewer",
                               const ColorFunctor& fcolor=ColorFunctor()) :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    ss2(ass2),
    m_fcolor(fcolor)
  {
    compute_elements();
  }

protected:
  /*
  void compute_face(Facet_const_handle fh)
  {
    CGAL::Color c=m_fcolor.run(ss2, fh);
    face_begin(c);

    add_point_in_face(fh->vertex(0)->point());
    add_point_in_face(fh->vertex(1)->point());
    add_point_in_face(fh->vertex(2)->point());

    face_end();
  }
  */
  void compute_edge(Halfedge_const_handle eh)
  {
    if(eh->is_bisector())
      add_segment(eh->opposite()->vertex()->point(), eh->vertex()->point(), CGAL::red());
    else
      add_segment(eh->opposite()->vertex()->point(), eh->vertex()->point(), CGAL::black());
  }
  void print_halfedge_labels(Halfedge_const_handle h)
  {
    std::stringstream label;
    label << "H" << h->id() << " (V" << h->vertex()->id() << ") ";
    label << "H" << h->opposite()->id() << " (V" << h->opposite()->vertex()->id() << ") ";
    add_text(CGAL::midpoint(h->opposite()->vertex()->point(), h->vertex()->point()), label.str());
  }

  void compute_vertex(Vertex_const_handle vh)
  {
    if(vh->is_split())
      add_point(vh->point(), CGAL::Color(10,10,180)); // blue, but not flashy
    else if(vh->has_infinite_time())
      add_point(vh->point(), CGAL::orange());
    else
      add_point(vh->point(), CGAL::Color(10,180,10)); // green, but not flashy
  }
  void print_vertex_label(Vertex_const_handle vh)
  {
    std::stringstream label;
    label << "V" << vh->id() << std::ends;
    add_text(vh->point(), label.str());
  }

  void compute_elements()
  {
    clear();

    for (typename SS2::Halfedge_const_iterator it=ss2.halfedges_begin(); it!=ss2.halfedges_end(); ++it)
    {
      if(it->id() < it->opposite()->id())
      {
        compute_edge(it);
        print_halfedge_labels(it);
      }
    }
    for (typename SS2::Vertex_const_iterator it=ss2.vertices_begin(); it!=ss2.vertices_end(); ++it)
    {
      compute_vertex(it);
      print_vertex_label(it);
    }
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
  const SS2& ss2;
  const ColorFunctor& m_fcolor;
};

// Specialization of draw function.
#define CGAL_SS_TYPE CGAL::Straight_skeleton_2<K>

template<class K>
void draw(const CGAL_SS_TYPE& ass2,
          const char* title="Straight Skeleton Basic Viewer")
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"ss2_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    DefaultColorFunctorSS2 fcolor;
    SimpleStraightSkeleton2ViewerQt<CGAL_SS_TYPE, DefaultColorFunctorSS2>
      mainwindow(app.activeWindow(), ass2, title, fcolor);
    mainwindow.show();
    app.exec();
  }
}

#undef CGAL_SS_TYPE

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_SS2_H
