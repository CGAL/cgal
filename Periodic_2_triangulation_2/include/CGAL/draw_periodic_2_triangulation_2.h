// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Jasmeet Singh <jasmeet.singh.mec11@iitbhu.ac.in>

#ifndef DRAW_PERIODIC_2_TRIANGULATION_2_H
#define DRAW_PERIODIC_2_TRIANGULATION_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Random.h>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorP2T2 {
  template <typename P2T2>
  static CGAL::Color run(const P2T2 &,
                         const typename P2T2::Periodic_triangle_iterator ti) {
    CGAL::Random random((unsigned int)(std::size_t)(&*ti));
    return get_random_color(random);
  }
};

// Viewer class for P2T2
template <class P2T2, class ColorFunctor>
class SimplePeriodic2Triangulation2ViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt                         Base;

  typedef typename P2T2::Iterator_type            Iterator_type;

  // Vertex iterators
  typedef typename P2T2::Periodic_point_iterator      Periodic_point_iterator;

  // Edge iterators
  typedef typename P2T2::Periodic_segment_iterator    Periodic_segment_iterator;
  typedef typename P2T2::Segment                      Segment;

  // Face iterators
  typedef typename P2T2::Periodic_triangle_iterator   Periodic_triangle_iterator;
  typedef typename P2T2::Triangle                     Triangle;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

public:
  /// Construct the viewer.
  /// @param ap2t2 the p2t2 to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this can be
  ///        usefull for very big object where this time could be long)
  SimplePeriodic2Triangulation2ViewerQt(QWidget* parent, const P2T2& ap2t2,
                               const char* title="Basic P2T2 Viewer",
                               bool anofaces=false,
                               const ColorFunctor& fcolor=ColorFunctor()) :
    // First draw: vertices; edges, faces; multi-color; no inverse normal
    Base(parent, title, true, true, true, false, false),
    p2t2(ap2t2),
    m_nofaces(anofaces),
    m_fcolor(fcolor),
    m_display_type(STORED_COVER_DOMAIN),
    m_domain(true)
  {
    // Add custom key description (see keyPressEvent).
    setKeyDescription(::Qt::Key_1, "STORED: Display all geometric primitives as they are stored in "
                                   "Triangulation_data_structure_2");
    setKeyDescription(::Qt::Key_2, "UNIQUE: Display only one representative of each geometric primitive even "
                                   "if the triangulation is computed in multiply sheeted covering space.");
    setKeyDescription(::Qt::Key_3, "STORED_COVER_DOMAIN: Same as STORED but also display "
                                   "all primitives whose intersection with the original "
                                   "domain of the current covering space is non-empty");
    setKeyDescription(::Qt::Key_4, "UNIQUE_COVER_DOMAIN: Same as UNIQUE but also display "
                                   "all primitives whose intersection with the original "
                                   "domain of the current covering space is non-empty");
    setKeyDescription(::Qt::Key_D, "Toggle 9-sheeted domain display");

    compute_elements();
  }
protected:
  void compute_vertex(Periodic_point_iterator pi)
  {
    // Construct the point in 9-sheeted covering space and add to viewer
    add_point(p2t2.point(*pi));
  }

  void compute_edge(Periodic_segment_iterator si)
  {
    // Construct the segment in 9-sheeted covering space and add to viewer
    Segment s(p2t2.segment(*si));
    add_segment(s[0], s[1]);
  }

  void compute_face(Periodic_triangle_iterator ti)
  {
    // Construct the triangle in 9-sheeted covering space and add to viewer
    Triangle t(p2t2.triangle(*ti));

    //CGAL::Color c=m_fcolor.run(p2t2, ti);
    face_begin(CGAL::Color(73, 250, 117));
    add_point_in_face(t[0]);
    add_point_in_face(t[1]);
    add_point_in_face(t[2]);
    face_end();

    // Display the edges of the faces as segments with a
    // light gray color for better visualization
    add_segment(t[0], t[1], CGAL::Color(207, 213, 211));
    add_segment(t[1], t[2], CGAL::Color(207, 213, 211));
    add_segment(t[2], t[0], CGAL::Color(207, 213, 211));
  }

  void compute_domain()
  {
    Kernel::Iso_rectangle_2 orig_domain =  p2t2.domain();
    std::array<int, 2> covering_sheets = p2t2.number_of_sheets();

    for(int i = 0; i < covering_sheets[0]; i++){
      for(int j = 0; j < covering_sheets[1]; j++){
        Kernel::Vector_2 shift(i * (orig_domain.xmax() - orig_domain.xmin()),
                               j * orig_domain.ymax() - orig_domain.ymin());
        Kernel::Point_2 p1(orig_domain.min());
        Kernel::Point_2 p2(orig_domain.xmin(), orig_domain.ymax());
        Kernel::Point_2 p3(orig_domain.xmax(), orig_domain.ymin());
        Kernel::Point_2 p4(orig_domain.max());

        add_segment(p1 + shift, p2 + shift, CGAL::Color(96, 104, 252));
        add_segment(p1 + shift, p3 + shift, CGAL::Color(96, 104, 252));
        add_segment(p2 + shift, p4 + shift, CGAL::Color(96, 104, 252));
        add_segment(p3 + shift, p4 + shift, CGAL::Color(96, 104, 252));
      }
    }
  }

  void compute_elements() {
    // Clear the buffers
    clear();

    // Get the display type, iterate through periodic elements according
    // to the display type
    Iterator_type it_type = (Iterator_type)m_display_type;

    // Iterate through vertices, edges and faces, add elements to buffer
    for (Periodic_point_iterator it =
             p2t2.periodic_points_begin(it_type);
         it != p2t2.periodic_points_end(it_type); it++)
    {
      compute_vertex(it);
    }

    for (Periodic_segment_iterator it =
             p2t2.periodic_segments_begin(it_type);
         it != p2t2.periodic_segments_end(it_type); it++)
    {
      compute_edge(it);
    }

    for (Periodic_triangle_iterator it =
             p2t2.periodic_triangles_begin(it_type);
         it != p2t2.periodic_triangles_end(it_type); it++)
    {
      compute_face(it);
    }

    if(m_domain){
      // Compute the (9-sheet covering space) domain of the periodic triangulation
      compute_domain();
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
    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    if (e->text()=="1")
    {
      m_display_type = Display_type::UNIQUE;
      displayMessage(QString("Display type= UNIQUE"));
      compute_elements();
      redraw();
    } else if (e->text()=="2")
    {
      m_display_type = Display_type::UNIQUE_COVER_DOMAIN;
      displayMessage(QString("Display type= UNIQUE_COVER_DOMAIN"));
      compute_elements();
      redraw();
    } else if (e->text()=="3")
    {
      m_display_type = Display_type::STORED;
      displayMessage(QString("Display type= STORED"));
      compute_elements();
      redraw();
    } else if (e->text()=="4")
    {
      m_display_type = Display_type::STORED_COVER_DOMAIN;
      displayMessage(QString("Display type= STORED_COVER_DOMAIN"));
      compute_elements();
      redraw();
    } else if ((e->key() == ::Qt::Key_D) && (modifiers == ::Qt::NoButton))
    {
      m_domain=!m_domain;
      displayMessage(QString("Draw domain=%1.").arg(m_domain?"true":"false"));
      compute_elements();
      redraw();
    } else {
      Base::keyPressEvent(e);
    }
  }

protected:
  const P2T2& p2t2;
  bool m_nofaces;
  const ColorFunctor& m_fcolor;
  enum Display_type
  {
    STORED = 0,
    UNIQUE, // 1
    STORED_COVER_DOMAIN, // 2
    UNIQUE_COVER_DOMAIN // 3
  };
  Display_type m_display_type;
  bool m_domain;
};

template<class P2T2, class ColorFunctor>
void draw(const P2T2& ap2t2,
          const char* title,
          bool nofill,
          const ColorFunctor& fcolor)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite=true;
#else
  bool cgal_test_suite=false;
#endif

  if (!cgal_test_suite)
  {
    int argc=1;
    const char* argv[2]={"p2t2_viewer","\0"};
    QApplication app(argc,const_cast<char**>(argv));
    SimplePeriodic2Triangulation2ViewerQt<P2T2, ColorFunctor> mainwindow(app.activeWindow(),
                                                              ap2t2,
                                                              title,
                                                              nofill,
                                                              fcolor);
    mainwindow.show();
    app.exec();
  }
}

template<class P2T2>
void draw(const P2T2& ap2t2, const char* title, bool nofill)
{
  DefaultColorFunctorP2T2 c;
  draw(ap2t2, title, nofill, c);
}

template<class P2T2>
void draw(const P2T2& ap2t2, const char* title)
{ draw(ap2t2, title, false); }

template<class P2T2>
void draw(const P2T2& ap2t2)
{ draw(ap2t2, "Basic 2D Periodic Triangulation Viewer"); }

} // namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DRAW_PERIODIC_2_TRIANGULATION_2_H
