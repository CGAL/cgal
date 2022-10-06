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

#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/license/Periodic_2_triangulation_2.h>
#include <CGAL/Qt/Basic_viewer_qt.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Random.h>

namespace CGAL {

namespace draw_function_for_P2T2
{

enum Display_type
{
  STORED = 0,
  UNIQUE, // 1
  STORED_COVER_DOMAIN, // 2
  UNIQUE_COVER_DOMAIN // 3
};

template <typename BufferType = float, class P2T2>
void compute_vertex(typename P2T2::Periodic_point_iterator pi, const P2T2 &p2t2, 
                    CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  // Construct the point in 9-sheeted covering space and add to viewer
  graphic_buffer.add_point(p2t2.point(*pi));
}

template <typename BufferType = float, class P2T2>
void compute_edge(typename P2T2::Periodic_segment_iterator si, const P2T2 &p2t2, 
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  typedef typename P2T2::Segment  Segment;

  // Construct the segment in 9-sheeted covering space and add to viewer
  Segment s(p2t2.segment(*si));
  graphic_buffer.add_segment(s[0], s[1]);
}

template <typename BufferType = float, class P2T2, class DrawingFunctor>
void compute_face(typename P2T2::Periodic_triangle_iterator ti, const P2T2 &p2t2, 
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer, 
                  const DrawingFunctor &m_drawing_functor)
{
  typedef typename P2T2::Triangle    Triangle;

  // Construct the triangle in 9-sheeted covering space and add to viewer
  Triangle t(p2t2.triangle(*ti));

  if(m_drawing_functor.colored_face(p2t2, ti)) {

    // Need CGAL::IO::Color(73, 250, 117);
    CGAL::IO::Color c= m_drawing_functor.face_color(p2t2, ti);
    graphic_buffer.face_begin(c);
    graphic_buffer.add_point_in_face(t[0]);
    graphic_buffer.add_point_in_face(t[1]);
    graphic_buffer.add_point_in_face(t[2]);
    graphic_buffer.face_end();
  }


  // Display the edges of the faces as segments with a
  // light gray color for better visualization
  if(m_drawing_functor.colored_face(p2t2, ti)) {

    // Need CGAL::IO::Color(207, 213, 211);
    CGAL::IO::Color c = m_drawing_functor.face_color(p2t2, ti);

    graphic_buffer.add_segment(t[0], t[1], c);
    graphic_buffer.add_segment(t[1], t[2], c);
    graphic_buffer.add_segment(t[2], t[0], c);
  }
}

template <typename BufferType = float, class P2T2>
void compute_domain(const P2T2 &p2t2, CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

  Kernel::Iso_rectangle_2 orig_domain =  p2t2.domain();
  std::array<int, 2> covering_sheets = p2t2.number_of_sheets();

  for(int i = 0; i < covering_sheets[0]; i++){
    for(int j = 0; j < covering_sheets[1]; j++){
      Kernel::Vector_2 shift(i * (orig_domain.xmax() - orig_domain.xmin()),
                              j * orig_domain.ymax() - orig_domain.ymin());
      Kernel::Point_2 p1((orig_domain.min)());
      Kernel::Point_2 p2(orig_domain.xmin(), orig_domain.ymax());
      Kernel::Point_2 p3(orig_domain.xmax(), orig_domain.ymin());
      Kernel::Point_2 p4((orig_domain.max)());

      graphic_buffer.add_segment(p1 + shift, p2 + shift, CGAL::IO::Color(96, 104, 252));
      graphic_buffer.add_segment(p1 + shift, p3 + shift, CGAL::IO::Color(96, 104, 252));
      graphic_buffer.add_segment(p2 + shift, p4 + shift, CGAL::IO::Color(96, 104, 252));
      graphic_buffer.add_segment(p3 + shift, p4 + shift, CGAL::IO::Color(96, 104, 252));
    }
  }
}

template <typename BufferType = float, class P2T2, class DrawingFunctor>
void compute_elements(const P2T2 &p2t2, CGAL::Graphic_buffer<BufferType> &graphic_buffer, const DrawingFunctor &m_drawing_functor, bool m_domain = true) {

  typedef typename P2T2::Iterator_type   Iterator_type;

  Display_type m_display_type(Display_type::STORED_COVER_DOMAIN);

  // Get the display type, iterate through periodic elements according
  // to the display type
  Iterator_type it_type = (Iterator_type)m_display_type;

  // Iterate through vertices, edges and faces, add elements to buffer
  for (typename P2T2::Periodic_point_iterator it =
            p2t2.periodic_points_begin(it_type);
        it != p2t2.periodic_points_end(it_type); it++)
  {
    compute_vertex(it, p2t2, graphic_buffer);
  }

  for (typename P2T2::Periodic_segment_iterator it =
            p2t2.periodic_segments_begin(it_type);
        it != p2t2.periodic_segments_end(it_type); it++)
  {
    compute_edge(it, p2t2, graphic_buffer);
  }

  for (typename P2T2::Periodic_triangle_iterator it =
            p2t2.periodic_triangles_begin(it_type);
        it != p2t2.periodic_triangles_end(it_type); it++)
  {
    compute_face(it, p2t2, graphic_buffer, m_drawing_functor);
  }

  if(m_domain){
    // Compute the (9-sheet covering space) domain of the periodic triangulation
    compute_domain(p2t2, graphic_buffer);
  }
}

} // namespace draw_function_for_P2T2

template <typename BufferType = float, class P2T2, class DrawingFunctor>
void add_in_graphic_buffer(const P2T2 &p2t2, CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                               const DrawingFunctor &m_drawing_functor) {
  draw_function_for_P2T2::compute_elements(p2t2, graphic_buffer, m_drawing_functor);
}

template <typename BufferType = float, class P2T2>
void add_in_graphic_buffer(const P2T2 &p2t2, CGAL::Graphic_buffer<BufferType> &graphic_buffer) {

  Drawing_functor<P2T2,typename P2T2::Periodic_triangle_iterator,
                  typename P2T2::Periodic_triangle_iterator,
                  typename P2T2::Periodic_triangle_iterator> drawing_functor;

  drawing_functor.colored_face = [](const P2T2&,
                      typename P2T2::Periodic_triangle_iterator) -> bool
  { return true; };

  // TODO: I think we need to add std::function like this:
  // drawing_functor.face_color =  [] (const P2T2&,
  //                          typename P2T2::Periodic_triangle_iterator,
  // int R, int G, int B)  -> CGAL::IO::Color
  // {
  //     return CGAL::IO::Color(R, G, B);
  // };  // What do you think?

  drawing_functor.face_color =  [] (const P2T2& alcc,
                           typename P2T2::Periodic_triangle_iterator dh) -> CGAL::IO::Color
  {
      return CGAL::IO::Color(207, 213, 211);
  };

  add_in_graphic_buffer(p2t2, graphic_buffer, drawing_functor);
}

// Specialization of draw function
#define CGAL_P2T2_TYPE CGAL::Periodic_2_triangulation_2<Gt, Tds >

template < class Gt,
           class Tds,
           class BufferType = float,
           class DrawingFunctor >
void draw(const CGAL_P2T2_TYPE& ap2t2,
          const DrawingFunctor &drawing_functor)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(ap2t2, buffer, drawing_functor);
  draw_buffer(buffer);
}

template < class Gt,
           class Tds,
           class BufferType = float >
void draw(const CGAL_P2T2_TYPE& ap2t2,
          const char* title = "2D Periodic Triangulation Viewer",
          bool nofill = false)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(ap2t2, buffer);
  draw_buffer(buffer);
}

} // namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // DRAW_PERIODIC_2_TRIANGULATION_2_H
