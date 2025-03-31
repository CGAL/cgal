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

#include <CGAL/point_generators_2.h>
#include <CGAL/Random.h>

#include <cmath>
#include <sstream>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorSS2
{
  template<typename SS2>
  static CGAL::IO::Color run(const SS2&,
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
  typedef Basic_viewer_qt                     Base;

  typedef typename SS2::Vertex_const_handle   Vertex_const_handle;
  typedef typename SS2::Halfedge_const_handle Halfedge_const_handle;
  typedef typename SS2::Face_const_handle     Face_const_handle;

  typedef typename SS2::Traits::FT            FT;
  typedef typename SS2::Traits::Point_2       Point;

public:
  /// Construct the viewer.
  /// @param ass2 the ss2 to view
  /// @param title the title of the window
  SimpleStraightSkeleton2ViewerQt(QWidget* parent, const SS2& ass2,
                                  const char* title="Basic SS2 Viewer",
                                  const ColorFunctor& fcolor = ColorFunctor()) :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    ss2(ass2),
    m_fcolor(fcolor)
  {
    compute_elements();
  }

protected:
  void compute_face(Halfedge_const_handle eh)
  {
    // Skip faces with infinite halfedges
    Halfedge_const_handle h = eh;
    do {
      if(h->vertex()->has_infinite_time()) {
        return;
      }
      h = h->next();
    } while(h != eh);

    CGAL::IO::Color c(150,170,170);

    face_begin(c);
    h = eh;
    do {
      add_point_in_face(h->vertex()->point());
      h = h->next();
    } while(h != eh);
    face_end();
  }

  void compute_edge(Halfedge_const_handle eh)
  {
    const bool src_inf = eh->opposite()->vertex()->has_infinite_time();
    const bool dst_inf = eh->vertex()->has_infinite_time();

    // if there is an inf vertex, put it at the destination
    const Point& src_p = src_inf ? eh->vertex()->point() : eh->opposite()->vertex()->point();
    Point dst_p = src_inf ? eh->opposite()->vertex()->point() : eh->vertex()->point();

    if(src_inf || dst_inf)
      dst_p = src_p + 0.1 * (dst_p - src_p);

    if(eh->is_bisector())
    {
      if (src_inf || dst_inf)
        add_segment(src_p, dst_p, CGAL::IO::orange());
      else
        add_segment(src_p, dst_p, CGAL::IO::red());
    }
    else
    {
      CGAL_assertion(!src_inf && !dst_inf);
      if (eh->weight() == 0)
        add_segment(src_p, dst_p, CGAL::IO::violet());
      else
        add_segment(src_p, dst_p, CGAL::IO::black());
    }
  }
  void print_halfedge_labels(Halfedge_const_handle h)
  {
    // Calculate length of the halfedge
    double edge_length = std::sqrt(CGAL::squared_distance(h->opposite()->vertex()->point(),
                                                          h->vertex()->point()));

    double radius = 0.01 * edge_length;

    // Generate random point on circle centered at the midpoint
    typedef CGAL::Random_points_on_circle_2<Point, CGAL::Creator_uniform_2<FT, Point> > Random_points;
    Point midpoint = CGAL::barycenter(h->opposite()->vertex()->point(), 0.75,
                                      h->vertex()->point(), 0.25);

    CGAL::Random rnd(h->id());
    Random_points random_point(radius, rnd);
    Point label_pos = midpoint + (*random_point - CGAL::ORIGIN);

    std::stringstream label;
    label << "H" << h->id() << " (V" << h->vertex()->id() << ") ";
    label << "H" << h->opposite()->id() << " (V" << h->opposite()->vertex()->id() << ") ";
    add_text(label_pos, label.str());
  }

  void compute_vertex(Vertex_const_handle vh)
  {
    if(vh->is_split())
      add_point(vh->point(), CGAL::IO::Color(10,10,180)); // blue, but not flashy
    else if(vh->has_infinite_time())
      add_point(vh->point(), CGAL::IO::orange());
    else if(vh->is_contour())
      add_point(vh->point(), CGAL::IO::black());
    else
      add_point(vh->point(), CGAL::IO::Color(10,180,10)); // green, but not flashy
  }
  void print_vertex_label(Vertex_const_handle vh)
  {
    // Calculate minimum incident edge length
    double min_length = std::numeric_limits<double>::max();
    auto h = vh->halfedge();
    do {
      if(!h->vertex()->has_infinite_time() && !h->opposite()->vertex()->has_infinite_time()) {
        double len = CGAL::squared_distance(h->vertex()->point(),
                                          h->opposite()->vertex()->point());
        min_length = std::min(min_length, std::sqrt(len));
      }
      h = h->next();
    } while(h != vh->halfedge());

    double radius = 0.001 * min_length;

    // Generate random point on circle
    typedef CGAL::Random_points_on_circle_2<Point, CGAL::Creator_uniform_2<FT, Point> > Random_points;

    CGAL::Random rnd(vh->id());
    Random_points random_point(radius, rnd);
    Point label_pos = vh->point() + (*random_point - CGAL::ORIGIN);

    std::stringstream label;
    label << "V" << vh->id() << std::ends;

    add_text(label_pos, label.str());
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

      if(!it->is_bisector() && !it->is_border())
        compute_face(it);
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
    const char* argv[2]={"ss2_viewer", nullptr};
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
