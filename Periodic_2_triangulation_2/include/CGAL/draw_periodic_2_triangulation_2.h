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
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_2_triangulation_2.h>

namespace CGAL {

// We need a specific graphics view options for periodic_2_triangulation_2 for the parameters
// of the domain.
template <typename DS,
          typename vertex_handle,
          typename edge_handle,
          typename face_handle>
struct Graphics_scene_options_periodic_2_triangulation_2 :
    public CGAL::Graphics_scene_options<DS, vertex_handle, edge_handle, face_handle>
{
  Graphics_scene_options_periodic_2_triangulation_2():
    m_domain_color(CGAL::IO::Color(96, 104, 252)),
    m_draw_domain(true),
    m_display_type(DS::STORED)
  {}

  bool draw_domain() const
  { return m_draw_domain; }
  void draw_domain(bool b)
  { m_draw_domain=b; }
  void toggle_draw_domain()
  { m_draw_domain=!m_draw_domain; }

  typename DS::Iterator_type display_type() const
  { return m_display_type; }

  void increase_display_type()
  {
    if(m_display_type==DS::UNIQUE_COVER_DOMAIN)
    { m_display_type=DS::STORED; }
    else
    { m_display_type=typename DS::Iterator_type(static_cast<int>(m_display_type)+1); }
  }

  const CGAL::IO::Color& domain_color() const
  { return m_domain_color; }
  void domain_color(const CGAL::IO::Color& c)
  { m_domain_color=c; }

protected:
  CGAL::IO::Color m_domain_color;
  bool m_draw_domain;
  typename DS::Iterator_type m_display_type;
};

namespace draw_function_for_P2T2
{

template <class P2T2, class GSOptions>
void compute_vertex(const P2T2 &p2t2,
                    typename P2T2::Periodic_point_iterator pi,
                    CGAL::Graphics_scene& graphics_scene,
                    const GSOptions& gs_options)
{
  // Construct the point in 9-sheeted covering space and add to viewer
  if(!gs_options.draw_vertex(p2t2, pi))
  { return; }

  if(gs_options.colored_vertex(p2t2, pi))
  { graphics_scene.add_point(p2t2.point(*pi),
                             gs_options.vertex_color(p2t2, pi)); }
  else
  { graphics_scene.add_point(p2t2.point(*pi)); }
}

template <class P2T2, class GSOptions>
void compute_edge(const P2T2 &p2t2,
                  typename P2T2::Periodic_segment_iterator si,
                  CGAL::Graphics_scene& graphics_scene,
                  const GSOptions& gs_options)
{
  if(!gs_options.draw_edge(p2t2, si))
  { return; }

  // Construct the segment in 9-sheeted covering space and add to viewer
  typename P2T2::Segment s(p2t2.segment(*si));
  if(gs_options.colored_edge(p2t2, si))
  { graphics_scene.add_segment(s[0], s[1],
                               gs_options.edge_color(p2t2, si)); }
  else
  { graphics_scene.add_segment(s[0], s[1]); }
}

template <class P2T2, class GSOptions>
void compute_face(const P2T2 &p2t2,
                  typename P2T2::Periodic_triangle_iterator ti,
                  CGAL::Graphics_scene& graphics_scene,
                  const GSOptions& gs_options)
{
  if(!gs_options.draw_face(p2t2, ti))
  { return; }

  // Construct the triangle in 9-sheeted covering space and add to viewer
  typename P2T2::Triangle t(p2t2.triangle(*ti));

  if(gs_options.colored_face(p2t2, ti))
  { graphics_scene.face_begin(gs_options.face_color(p2t2, ti)); }
  else
  { graphics_scene.face_begin(); }

  graphics_scene.add_point_in_face(t[0]);
  graphics_scene.add_point_in_face(t[1]);
  graphics_scene.add_point_in_face(t[2]);
  graphics_scene.face_end();
}

template <class P2T2, class GSOptions>
void compute_domain(const P2T2& p2t2,
                    CGAL::Graphics_scene& graphics_scene,
                    const GSOptions& gs_options)
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

      graphics_scene.add_segment(p1 + shift, p2 + shift, gs_options.domain_color());
      graphics_scene.add_segment(p1 + shift, p3 + shift, gs_options.domain_color());
      graphics_scene.add_segment(p2 + shift, p4 + shift, gs_options.domain_color());
      graphics_scene.add_segment(p3 + shift, p4 + shift, gs_options.domain_color());
    }
  }
}

template <class P2T2, class GSOptions>
void compute_elements(const P2T2& p2t2,
                      CGAL::Graphics_scene& graphics_scene,
                      const GSOptions& gs_options)
{
  // Get the display type, iterate through periodic elements according
  // to the display type
  typedef typename P2T2::Iterator_type Iterator_type;
  Iterator_type it_type = (Iterator_type)gs_options.display_type();

  // Iterate through vertices, edges and faces, add elements to buffer
  if(gs_options.are_vertices_enabled())
  {
    for (typename P2T2::Periodic_point_iterator it=p2t2.periodic_points_begin(it_type);
         it!=p2t2.periodic_points_end(it_type); ++it)
    { compute_vertex(p2t2, it, graphics_scene, gs_options); }
  }

  if(gs_options.are_edges_enabled())
  {
    for (typename P2T2::Periodic_segment_iterator it=p2t2.periodic_segments_begin(it_type);
         it!=p2t2.periodic_segments_end(it_type); ++it)
    { compute_edge(p2t2, it, graphics_scene, gs_options); }
  }

  if (gs_options.are_faces_enabled())
  {
    for (typename P2T2::Periodic_triangle_iterator it=p2t2.periodic_triangles_begin(it_type);
         it!=p2t2.periodic_triangles_end(it_type); ++it)
    { compute_face(p2t2, it, graphics_scene, gs_options); }
  }

  if(gs_options.draw_domain())
  {
    // Compute the (9-sheet covering space) domain of the periodic triangulation
    compute_domain(p2t2, graphics_scene, gs_options);
  }
}

} // namespace draw_function_for_P2T2

#define CGAL_P2T2_TYPE CGAL::Periodic_2_triangulation_2<Gt, Tds >

template <class Gt, class Tds, class GSOptions>
void add_to_graphics_scene(const CGAL_P2T2_TYPE& p2t2,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& gs_options)
{
  draw_function_for_P2T2::compute_elements(p2t2, graphics_scene, gs_options);
}

template <class Gt, class Tds>
void add_to_graphics_scene(const CGAL_P2T2_TYPE& p2t2,
                           CGAL::Graphics_scene& graphics_scene)
{
  CGAL::Graphics_scene_options_periodic_2_triangulation_2
    <CGAL_P2T2_TYPE,
     typename CGAL_P2T2_TYPE::Periodic_point_iterator,
     typename CGAL_P2T2_TYPE::Periodic_segment_iterator,
     typename CGAL_P2T2_TYPE::Periodic_triangle_iterator> gs_options;

  add_to_graphics_scene(p2t2, graphics_scene, gs_options);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function
template<class Gt, class Tds, class GSOptions>
void draw(const CGAL_P2T2_TYPE& ap2t2,
          GSOptions& gs_options,
          const char* title="2D Periodic Triangulation Viewer")
{
  CGAL::Graphics_scene gs;
  add_to_graphics_scene(ap2t2, gs, gs_options);
  CGAL::Qt::QApplication_and_basic_viewer app(gs, title);
  if(app)
  {
    // Here we define the std::function to capture key pressed.
    app.basic_viewer().on_key_pressed=
      [&ap2t2, &gs, &gs_options] (QKeyEvent* e, CGAL::Qt::Basic_viewer* basic_viewer) -> bool
      {
        const ::Qt::KeyboardModifiers modifiers = e->modifiers();
        if ((e->key() == ::Qt::Key_D) && (modifiers == ::Qt::NoButton))
        {
          gs_options.increase_display_type();
          basic_viewer->displayMessage
            (QString("Display type=%1.").arg(gs_options.display_type()==0?"Stored":
                                             (gs_options.display_type()==1?"Unique":
                                              (gs_options.display_type()==2?"Stored cover":
                                               "Unique cover"))));
          gs.clear();
          add_to_graphics_scene(ap2t2, gs, gs_options);
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

template<class Gt, class Tds>
void draw(const CGAL_P2T2_TYPE& ap2t2,
          const char* title="2D Periodic Triangulation Viewer")
{
  CGAL::Graphics_scene_options_periodic_2_triangulation_2
    <CGAL_P2T2_TYPE,
     typename CGAL_P2T2_TYPE::Periodic_point_iterator,
     typename CGAL_P2T2_TYPE::Periodic_segment_iterator,
     typename CGAL_P2T2_TYPE::Periodic_triangle_iterator> gs_options;
  draw(ap2t2, gs_options, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_P2T2_TYPE

} // namespace CGAL

#endif // DRAW_PERIODIC_2_TRIANGULATION_2_H
