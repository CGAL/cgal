// Copyright (c) 1997
// Utrecht University (The Netherlands),
// ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France),
// Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).  All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>

#ifndef CGAL_DRAW_POLYGON_SET_2_H
#define CGAL_DRAW_POLYGON_SET_2_H

#include <CGAL/license/Boolean_set_operations_2.h>

#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Polygon_set_2.h>

namespace CGAL {

template <typename DS,
          typename VertexDescriptor,
          typename EdgeDescriptor,
          typename FaceDescriptor>
struct Graphics_scene_options_polygon_set_2 :
    public CGAL::Graphics_scene_options<DS, VertexDescriptor, EdgeDescriptor, FaceDescriptor>
{
  void unbounded_face_color(const CGAL::IO::Color& c)
  { m_unbounded_face_color=c; }
  const CGAL::IO::Color& unbounded_face_color() const
  { return m_unbounded_face_color; }

  bool draw_unbounded() const
  { return m_draw_unbounded; }
  void draw_unbounded(bool b)
  { m_draw_unbounded=b; }
  void toggle_draw_unbounded()
  { m_draw_unbounded=!m_draw_unbounded; }

  int height() const
  { return m_height; }
  int width() const
  { return m_width; }

  void height(int i)
  { m_height=i; }
  void width(int i)
  { m_width=i; }

protected:
  bool m_draw_unbounded=true;
  CGAL::IO::Color m_unbounded_face_color=CGAL::IO::Color(75,160,255);
  int m_width=0, m_height=0;
};

namespace draw_function_for_boolean_set_2 {

template <typename PS2, class GSOptions>
void compute_loop(const typename PS2::Polygon_2& p, bool hole,
                  CGAL::Graphics_scene& gs,
                  const GSOptions& gso)
{
  if (gso.are_faces_enabled() && hole)
  { gs.add_point_in_face(p.vertex(p.size()-1)); }

  auto prev = p.vertices_begin();
  auto it = prev;

  // First vertex
  if (gso.are_vertices_enabled() && gso.draw_vertex(p, it))
  {
    if(gso.colored_vertex(p, it))
    { gs.add_point(*it, gso.vertex_color(p, it)); }
    else
    { gs.add_point(*it); }
  }

  if (gso.are_faces_enabled())
  { gs.add_point_in_face(*it); }

  // Other vertices
  for (++it; it!=p.vertices_end(); ++it)
  {
    if (gso.are_vertices_enabled() && gso.draw_vertex(p, it))
    { // Add point
      if(gso.colored_vertex(p, it))
      { gs.add_point(*it, gso.vertex_color(p, it)); }
      else
      { gs.add_point(*it); }
    }

    if (gso.are_edges_enabled() && gso.draw_edge(p, prev))
    { // Add segment with previous point
      if(gso.colored_edge(p, prev))
      { gs.add_segment(*prev, *it, gso.edge_color(p, prev)); }
      else
      { gs.add_segment(*prev, *it); }
    }

    if (gso.are_faces_enabled())
    { gs.add_point_in_face(*it); } // add point in face
    prev = it;
  }

  // Last segment (between last and first point)
  if (gso.are_edges_enabled() && gso.draw_edge(p, prev))
  {
    if(gso.colored_edge(p, prev))
    { gs.add_segment(*prev, *(p.vertices_begin()), gso.edge_color(p, prev)); }
    else
    { gs.add_segment(*prev, *(p.vertices_begin())); }
  }
}

/// Compute the elements of a polygon with holes.
template <typename PWH, class GSOptions>
void compute_elements(const PWH& pwh,
                      CGAL::Graphics_scene& gs,
                      const GSOptions& gso)
{
  if (!gso.draw_unbounded() && pwh.outer_boundary().is_empty()) return;

  if (gso.are_faces_enabled())
  { gs.face_begin(gso.unbounded_face_color()); }

  using Pnt=typename PWH::Polygon_2::Point_2;
  const Pnt* point_in_face=nullptr;
  if (pwh.outer_boundary().is_empty())
  {
    if(gso.width()!=0 && gso.height()!=0)
    {
      typename PWH::Polygon_2 pgn;
      pgn.push_back(Pnt(-gso.width(), -gso.height()));
      pgn.push_back(Pnt(gso.width(), -gso.height()));
      pgn.push_back(Pnt(gso.width(), gso.height()));
      pgn.push_back(Pnt(-gso.width(), gso.height()));
      draw_function_for_boolean_set_2::compute_loop<PWH>(pgn, false, gs, gso);
      point_in_face = &(pgn.vertex(pgn.size()-1));
    }
  }
  else
  {
    const auto& outer_boundary = pwh.outer_boundary();
    draw_function_for_boolean_set_2::compute_loop<PWH>(outer_boundary, false, gs, gso);
    point_in_face = &(outer_boundary.vertex(outer_boundary.size()-1));
  }

  for (auto it = pwh.holes_begin(); it != pwh.holes_end(); ++it)
  {
    draw_function_for_boolean_set_2::compute_loop<PWH>(*it, true, gs, gso);
    if (gso.are_faces_enabled() && point_in_face && point_in_face!=nullptr)
    { gs.add_point_in_face(*point_in_face); }
  }

  if (gso.are_faces_enabled())
  { gs.face_end(); }
}

} // End namespace draw_function_for_boolean_set_2

#ifdef CGAL_USE_BASIC_VIEWER

template <typename PolygonSet_2, typename GSOptions>
class Polygon_set_2_basic_viewer_qt : public Basic_viewer
{
  using Base = Basic_viewer;
  using Ps = PolygonSet_2;
  using Pwh = typename Ps::Polygon_with_holes_2;
  using Pgn = typename Ps::Polygon_2;
  using Pnt = typename Pgn::Point_2;

public:
  Polygon_set_2_basic_viewer_qt(QWidget* parent, const Ps& ps,
                                GSOptions& gs_options,
                                const char* title = "Basic Polygon_set_2 Viewer") :
    Base(parent, graphics_scene, title),
    m_ps(ps),
    gso(gs_options)
  {
    gso.width(CGAL_BASIC_VIEWER_INIT_SIZE_X);
    gso.height(CGAL_BASIC_VIEWER_INIT_SIZE_Y);

    if (ps.is_empty()) return;

    // mimic the computation of Camera::pixelGLRatio()
    auto bbox = bounding_box();
    CGAL::qglviewer::Vec minv(bbox.xmin(), bbox.ymin(), 0);
    CGAL::qglviewer::Vec maxv(bbox.xmax(), bbox.ymax(), 0);
    auto diameter = (maxv - minv).norm();
    m_pixel_ratio = diameter / gso.height();
  }

  /*! Intercept the resizing of the window.
   */
  virtual void resizeGL(int width, int height) {
    CGAL::QGLViewer::resizeGL(width, height);
    gso.width(width);
    gso.height(height);
    CGAL::qglviewer::Vec p;
    auto ratio = this->camera()->pixelGLRatio(p);
    if (ratio != m_pixel_ratio)
    {
      m_pixel_ratio = ratio;
      add_elements();
    }
  }

  /*! Obtain the pixel ratio.
   */
  double pixel_ratio() const { return m_pixel_ratio; }

  /*! Compute the elements of a polygon set.
   */
  virtual void add_elements()
  {
    graphics_scene.clear();
    std::vector<Pwh> pwhs;
    m_ps.polygons_with_holes(std::back_inserter(pwhs));
    for (const auto& pwh : pwhs)
    { draw_function_for_boolean_set_2::compute_elements(pwh, graphics_scene, gso); }
  }

  /*! Compute the bounding box.
   */
  CGAL::Bbox_2 bounding_box() {
    Bbox_2 bbox;
    std::vector<Pwh> pwhs;
    m_ps.polygons_with_holes(std::back_inserter(pwhs));
    for (const auto& pwh : pwhs) {
      if (! pwh.is_unbounded()) {
        bbox += pwh.outer_boundary().bbox();
        continue;
      }
      for (auto it = pwh.holes_begin(); it != pwh.holes_end(); ++it)
        bbox += it->bbox();
    }
    return bbox;
  }

private:
  //! The ratio between pixel and opengl units (in world coordinate system).
  double m_pixel_ratio = 1;

  //! The polygon set to draw.
  const Ps& m_ps;

  //! Indicates whether to draw unbounded polygons with holes.
  bool m_draw_unbounded = false;

  Graphics_scene graphics_scene;
  GSOptions& gso;
};

#endif // CGAL_USE_BASIC_VIEWER

#define CGAL_PS2_TYPE CGAL::Polygon_set_2<T, C, D>

// Specializations of add_to_graphics_scene function
template<class T, class C, class D, class GSOptions>
void add_to_graphics_scene(const CGAL_PS2_TYPE& ap2,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& gso)
{ draw_function_for_boolean_set_2::compute_elements(ap2, graphics_scene, gso); }

template<class T, class C, class D>
void add_to_graphics_scene(const CGAL_PS2_TYPE& ap2,
                           CGAL::Graphics_scene &graphics_scene)
{
  CGAL::Graphics_scene_options_polygon_set_2<typename CGAL_PS2_TYPE::Polygon_2,
                        typename CGAL_PS2_TYPE::Vertex_const_iterator,
                        typename CGAL_PS2_TYPE::Vertex_const_iterator,
                        void*> gso;
  draw_function_for_boolean_set_2::compute_elements(ap2, graphics_scene, gso);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function.
template<class T, class C, class D, class GSOptions>
void draw(const CGAL_PS2_TYPE& ps, GSOptions& gso,
          const char* title = "Polygon_set_2 Basic Viewer")
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (! cgal_test_suite)
  {
    using Ps = CGAL::Polygon_set_2<T, C, D>;
    using Viewer = Polygon_set_2_basic_viewer_qt<Ps, GSOptions>;
    CGAL::Qt::init_ogl_context(4,3);
    int argc = 1;
    const char* argv[2] = {"t2_viewer", nullptr};
    QApplication app(argc, const_cast<char**>(argv));
    Viewer basic_viewer(app.activeWindow(), ps, gso, title);
    basic_viewer.show();
    app.exec();
  }
}

template<class T, class C, class D>
void draw(const CGAL_PS2_TYPE& ps,
          const char* title = "Polygon_set_2 Basic Viewer")
{
  CGAL::Graphics_scene_options_polygon_set_2<typename CGAL_PS2_TYPE::Polygon_2,
                                             typename CGAL_PS2_TYPE::Polygon_2::Vertex_const_iterator,
                                             typename CGAL_PS2_TYPE::Polygon_2::Vertex_const_iterator,
                               void*> gso;
  draw(ps, gso, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_PS2_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_POLYGON_SET_2_H
