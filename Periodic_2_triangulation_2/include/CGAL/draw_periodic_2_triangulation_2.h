// Copyright (c) 2019 INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Jasmeet Singh <jasmeet.singh.mec11@iitbhu.ac.in>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef DRAW_PERIODIC_2_TRIANGULATION_2_H
#define DRAW_PERIODIC_2_TRIANGULATION_2_H

#include <CGAL/license/Periodic_2_triangulation_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/Graphic_storage.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_2.h>

namespace CGAL {

// We need a specific functor for periodic_2_triangulation_2 for the parameter
// of the domain.
template <typename DS,
          typename vertex_handle,
          typename edge_handle,
          typename face_handle>
struct Drawing_functor_periodic_2_triangulation_2 :
    public CGAL::Drawing_functor<DS, vertex_handle, edge_handle, face_handle>
{
  Drawing_functor_periodic_2_triangulation_2():
    domain_color(CGAL::IO::Color(96, 104, 252)),
    m_draw_domain(true),
    m_current_display_type(DS::STORED)
  {}

  bool get_draw_domain() const
  { return m_draw_domain; }
  void set_draw_domain(bool b)
  { m_draw_domain=b; }
  void negate_draw_domain()
  { m_draw_domain=!m_draw_domain; }

  typename DS::Iterator_type current_display_type() const
  { return m_current_display_type; }

  void increase_current_display_type()
  {
    if(m_current_display_type==DS::UNIQUE_COVER_DOMAIN)
    { m_current_display_type=DS::STORED; }
    else
    { m_current_display_type=typename DS::Iterator_type(static_cast<int>(m_current_display_type)+1); }
  }

public:
  CGAL::IO::Color domain_color;

protected:
  bool m_draw_domain;
  typename DS::Iterator_type m_current_display_type;
};

namespace draw_function_for_P2T2
{

template <typename BufferType=float, class P2T2, class DrawingFunctor>
void compute_vertex(const P2T2 &p2t2,
                    typename P2T2::Periodic_point_iterator pi,
                    CGAL::Graphic_storage<BufferType>& graphic_storage,
                    const DrawingFunctor& drawing_functor)
{
  // Construct the point in 9-sheeted covering space and add to viewer
  if(!drawing_functor.draw_vertex(p2t2, pi))
  { return; }

  if(drawing_functor.colored_vertex(p2t2, pi))
  { graphic_storage.add_point(p2t2.point(*pi),
                             drawing_functor.vertex_color(p2t2, pi)); }
  else
  { graphic_storage.add_point(p2t2.point(*pi)); }
}

template <typename BufferType=float, class P2T2, class DrawingFunctor>
void compute_edge(const P2T2 &p2t2,
                  typename P2T2::Periodic_segment_iterator si,
                  CGAL::Graphic_storage<BufferType>& graphic_storage,
                  const DrawingFunctor& drawing_functor)
{
  if(!drawing_functor.draw_edge(p2t2, si))
  { return; }

  // Construct the segment in 9-sheeted covering space and add to viewer
  typename P2T2::Segment s(p2t2.segment(*si));
  if(drawing_functor.colored_edge(p2t2, si))
  { graphic_storage.add_segment(s[0], s[1],
                               drawing_functor.edge_color(p2t2, si)); }
  else
  { graphic_storage.add_segment(s[0], s[1]); }
}

template <typename BufferType=float, class P2T2, class DrawingFunctor>
void compute_face(const P2T2 &p2t2,
                  typename P2T2::Periodic_triangle_iterator ti,
                  CGAL::Graphic_storage<BufferType>& graphic_storage,
                  const DrawingFunctor& drawing_functor)
{
  if(!drawing_functor.draw_face(p2t2, ti))
  { return; }

  // Construct the triangle in 9-sheeted covering space and add to viewer
  typename P2T2::Triangle t(p2t2.triangle(*ti));

  if(drawing_functor.colored_face(p2t2, ti))
  { graphic_storage.face_begin(drawing_functor.face_color(p2t2, ti)); }
  else
  { graphic_storage.face_begin(); }

  graphic_storage.add_point_in_face(t[0]);
  graphic_storage.add_point_in_face(t[1]);
  graphic_storage.add_point_in_face(t[2]);
  graphic_storage.face_end();
}

template <typename BufferType=float, class P2T2, class DrawingFunctor>
void compute_domain(const P2T2& p2t2,
                    CGAL::Graphic_storage<BufferType>& graphic_storage,
                    const DrawingFunctor& drawing_functor)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

  Kernel::Iso_rectangle_2 orig_domain =  p2t2.domain();
  std::array<int, 2> covering_sheets = p2t2.number_of_sheets();

  for(int i = 0; i < covering_sheets[0]; i++)
  {
    for(int j = 0; j < covering_sheets[1]; j++)
    {
      Kernel::Vector_2 shift(i * (orig_domain.xmax() - orig_domain.xmin()),
                             j * orig_domain.ymax() - orig_domain.ymin());
      Kernel::Point_2 p1((orig_domain.min)());
      Kernel::Point_2 p2(orig_domain.xmin(), orig_domain.ymax());
      Kernel::Point_2 p3(orig_domain.xmax(), orig_domain.ymin());
      Kernel::Point_2 p4((orig_domain.max)());

      graphic_storage.add_segment(p1 + shift, p2 + shift, drawing_functor.domain_color);
      graphic_storage.add_segment(p1 + shift, p3 + shift, drawing_functor.domain_color);
      graphic_storage.add_segment(p2 + shift, p4 + shift, drawing_functor.domain_color);
      graphic_storage.add_segment(p3 + shift, p4 + shift, drawing_functor.domain_color);
    }
  }
}

template <typename BufferType=float, class P2T2, class DrawingFunctor>
void compute_elements(const P2T2& p2t2,
                      CGAL::Graphic_storage<BufferType>& graphic_storage,
                      const DrawingFunctor& drawing_functor)
{
  // Get the display type, iterate through periodic elements according
  // to the display type
  typedef typename P2T2::Iterator_type Iterator_type;
  Iterator_type it_type = (Iterator_type)drawing_functor.current_display_type();

  // Iterate through vertices, edges and faces, add elements to buffer
  if(drawing_functor.are_vertices_enabled())
  {
    for (typename P2T2::Periodic_point_iterator it=p2t2.periodic_points_begin(it_type);
         it!=p2t2.periodic_points_end(it_type); ++it)
    { compute_vertex(p2t2, it, graphic_storage, drawing_functor); }
  }

  if(drawing_functor.are_edges_enabled())
  {
    for (typename P2T2::Periodic_segment_iterator it=p2t2.periodic_segments_begin(it_type);
         it!=p2t2.periodic_segments_end(it_type); ++it)
    { compute_edge(p2t2, it, graphic_storage, drawing_functor); }
  }

  if (drawing_functor.are_faces_enabled())
  {
    for (typename P2T2::Periodic_triangle_iterator it=p2t2.periodic_triangles_begin(it_type);
         it!=p2t2.periodic_triangles_end(it_type); ++it)
    { compute_face(p2t2, it, graphic_storage, drawing_functor); }
  }

  if(drawing_functor.get_draw_domain())
  {
    // Compute the (9-sheet covering space) domain of the periodic triangulation
    compute_domain(p2t2, graphic_storage, drawing_functor);
  }
}

} // namespace draw_function_for_P2T2

#define CGAL_P2T2_TYPE CGAL::Periodic_2_triangulation_2<Gt, Tds >

template <typename BufferType=float, class Gt, class Tds, class DrawingFunctor>
void add_in_graphic_storage(const CGAL_P2T2_TYPE& p2t2,
                           CGAL::Graphic_storage<BufferType>& graphic_storage,
                           const DrawingFunctor& drawing_functor)
{
  draw_function_for_P2T2::compute_elements(p2t2, graphic_storage, drawing_functor);
}

template <typename BufferType=float, class Gt, class Tds>
void add_in_graphic_storage(const CGAL_P2T2_TYPE& p2t2,
                           CGAL::Graphic_storage<BufferType>& graphic_storage)
{
  CGAL::Drawing_functor_periodic_2_triangulation_2
    <CGAL_P2T2_TYPE,
     typename CGAL_P2T2_TYPE::Periodic_point_iterator,
     typename CGAL_P2T2_TYPE::Periodic_segment_iterator,
     typename CGAL_P2T2_TYPE::Periodic_triangle_iterator> drawing_functor;

  add_in_graphic_storage(p2t2, graphic_storage, drawing_functor);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function
template<class Gt, class Tds, class BufferType=float, class DrawingFunctor>
void draw(const CGAL_P2T2_TYPE& ap2t2,
          const DrawingFunctor& drawing_functor,
          const char* title="2D Periodic Triangulation Viewer")
{
  CGAL::Graphic_storage<BufferType> buffer;
  add_in_graphic_storage(ap2t2, buffer, drawing_functor);
  draw_graphic_storage(buffer);
}

template<class Gt, class Tds, class BufferType=float>
void draw(const CGAL_P2T2_TYPE& ap2t2,
          const char* title="2D Periodic Triangulation Viewer")
{
  CGAL::Graphic_storage<BufferType> buffer;
  CGAL::Drawing_functor_periodic_2_triangulation_2
    <CGAL_P2T2_TYPE,
     typename CGAL_P2T2_TYPE::Periodic_point_iterator,
     typename CGAL_P2T2_TYPE::Periodic_segment_iterator,
     typename CGAL_P2T2_TYPE::Periodic_triangle_iterator> drawing_functor;

  add_in_graphic_storage(ap2t2, buffer, drawing_functor);
  QApplication_and_basic_viewer app(buffer, title);
  if(app)
  {
    // Here we define the std::function to capture key pressed.
    app.basic_viewer().on_key_pressed=
      [&ap2t2, &drawing_functor] (QKeyEvent* e, CGAL::Basic_viewer_qt<float>* basic_viewer) -> bool
      {
        const ::Qt::KeyboardModifiers modifiers = e->modifiers();
        if ((e->key() == ::Qt::Key_D) && (modifiers == ::Qt::NoButton))
        {
          drawing_functor.increase_current_display_type();
          basic_viewer->displayMessage
            (QString("Display type=%1.").arg(drawing_functor.current_display_type()==0?"Stored":
                                             (drawing_functor.current_display_type()==1?"Unique":
                                              (drawing_functor.current_display_type()==2?"Stored cover":
                                               "Unique cover"))));
          basic_viewer->clear();
          draw_function_for_P2T2::compute_elements(ap2t2,
                                                   basic_viewer->get_graphic_storage(),
                                                   drawing_functor);
          basic_viewer->redraw();
        }
        else
        {
          // Return false will call the base method to process others/classicals key
          return false;
        }
        return true; // the key was captured
      };

    // Here we add shortcut descriptions
    app.basic_viewer().setKeyDescription(::Qt::Key_D, "Next display type");

    // Then we run the app
    app.run();
  }
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_P2T2_TYPE

} // namespace CGAL

#endif // DRAW_PERIODIC_2_TRIANGULATION_2_H
