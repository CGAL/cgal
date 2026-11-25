// Copyright (c) 2018  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
// Author(s)     : Guillaume Damiand <guillaume.damiand@liris.cnrs.fr>
//                 Mostafa Ashraf <mostaphaashraf1996@gmail.com>

#ifndef CGAL_DRAW_SS2_H
#define CGAL_DRAW_SS2_H

#include <CGAL/license/Straight_skeleton_2.h>

#include <CGAL/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Straight_skeleton_2.h>

#include <CGAL/Bbox_2.h>
#include <CGAL/number_utils.h>
#include <CGAL/Random.h>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <sstream>
#include <string>

namespace CGAL {

namespace draw_function_for_ss2 {

struct Drawing_params
{
  double incomplete_bisector_length = 1.;
  double move_radius = 1e-3;
};

template <class SS2>
Drawing_params compute_drawing_parameters(const SS2& ss2)
{
  Bbox_2 bbox;
  for(typename SS2::Vertex_const_iterator it=ss2.vertices_begin(); it!=ss2.vertices_end(); ++it)
  {
    if(it->has_infinite_time())
      continue;

    bbox += it->point().bbox();
  }

  Drawing_params ctx;
  const double dx = bbox.xmax() - bbox.xmin();
  const double dy = bbox.ymax() - bbox.ymin();
  const double diag = std::sqrt(dx*dx + dy*dy);
  if(diag > 0)
  {
    ctx.incomplete_bisector_length = 0.1 * diag;
    ctx.move_radius = 1e-3 * diag;
  }

  return ctx;
}

// Direction (not normalized) followed in time by the wavefront vertex tracing the bisector `h`.
//
// The two faces incident to a bisector are each defined by a contour edge, whose supporting
// line moves inward (i.e., to its left) at a speed equal to its weight. The direction `d`
// of the bisector thus satisfies `n_i . d = w_i` for both defining contour edges,
// with `n_i` the inward unit normal and `w_i` the weight of the edge.
//
// Returns `false` if the direction cannot be determined (e.g., parallel defining edges).
template <class SS2>
bool bisector_direction(typename SS2::Halfedge_const_handle h,
                        double& dx, double& dy)
{
  CGAL_precondition(h->is_bisector());

  typename SS2::Halfedge_const_handle c1 = h->defining_contour_edge();
  typename SS2::Halfedge_const_handle c2 = h->opposite()->defining_contour_edge();

  auto inward_normal = [](typename SS2::Halfedge_const_handle c, double& nx, double& ny) -> bool
  {
    const auto& s = c->opposite()->vertex()->point();
    const auto& t = c->vertex()->point();
    const double ex = CGAL::to_double(t.x()) - CGAL::to_double(s.x());
    const double ey = CGAL::to_double(t.y()) - CGAL::to_double(s.y());
    const double l = std::sqrt(ex*ex + ey*ey);
    if(!(l > 0))
      return false;

    // the face is on the left of its defining contour halfedge
    nx = -ey / l;
    ny =  ex / l;
    return true;
  };

  double n1x, n1y, n2x, n2y;
  if(!inward_normal(c1, n1x, n1y) || !inward_normal(c2, n2x, n2y))
    return false;

  const double w1 = CGAL::to_double(c1->weight());
  const double w2 = CGAL::to_double(c2->weight());

  const double det = n1x * n2y - n1y * n2x;
  if(std::abs(det) < 1e-12)
    return false;

  dx = ( w1 * n2y - w2 * n1y) / det;
  dy = (-w1 * n2x + w2 * n1x) / det;

  return (dx != 0. || dy != 0.);
}

template <class SS2>
bool vertex_position(typename SS2::Vertex_const_handle vh,
                     const Drawing_params& ctx,
                     double& x, double& y)
{
  if(!vh->has_infinite_time())
  {
    x = CGAL::to_double(vh->point().x());
    y = CGAL::to_double(vh->point().y());
    return true;
  }

  typename SS2::Halfedge_const_handle h = vh->halfedge(); // incoming halfedge
  if(h == typename SS2::Halfedge_const_handle() || !h->is_bisector())
    return false;

  typename SS2::Vertex_const_handle src = h->opposite()->vertex();
  if(src->has_infinite_time())
    return false;

  double dx, dy;
  if(!bisector_direction<SS2>(h, dx, dy))
    return false;

  const double l = std::sqrt(dx*dx + dy*dy);
  x = CGAL::to_double(src->point().x()) + ctx.incomplete_bisector_length * dx / l;
  y = CGAL::to_double(src->point().y()) + ctx.incomplete_bisector_length * dy / l;

  return true;
}

template <class SS2>
typename SS2::Traits::Point_2 make_point(double x, double y)
{
  using FT = typename SS2::Traits::FT;
  return typename SS2::Traits::Point_2{FT(x), FT(y)};
}

// clunky try to avoid labels on top of one another
inline void perturb_position(int seed, double radius, double& x, double& y)
{
  CGAL::Random rnd(static_cast<unsigned int>(seed));
  const double angle = rnd.get_double(0., 2. * CGAL_PI);
  x += radius * std::cos(angle);
  y += radius * std::sin(angle);
}

template <class SS2, class GSOptions>
void compute_face(const SS2& ss2, typename SS2::Face_const_handle fh,
                  CGAL::Graphics_scene& graphics_scene,
                  const GSOptions& gs_options)
{
  if(!gs_options.draw_face(ss2, fh))
    return;

  typename SS2::Halfedge_const_handle start = fh->halfedge();
  if(start == typename SS2::Halfedge_const_handle())
    return;

  // only finite faces are drawn
  typename SS2::Halfedge_const_handle h = start;
  do
  {
    if(h->vertex()->has_infinite_time())
      return;
    h = h->next();
  }
  while(h != start);

  if(gs_options.colored_face(ss2, fh))
    graphics_scene.face_begin(gs_options.face_color(ss2, fh));
  else
    graphics_scene.face_begin();

  h = start;
  do
  {
    graphics_scene.add_point_in_face(h->vertex()->point());
    h = h->next();
  }
  while(h != start);

  graphics_scene.face_end();
}

template <class SS2, class GSOptions>
void compute_edge(const SS2& ss2, typename SS2::Halfedge_const_handle eh,
                  CGAL::Graphics_scene& graphics_scene,
                  const GSOptions& gs_options,
                  const Drawing_params& ctx)
{
  if(!gs_options.draw_edge(ss2, eh))
    return;

  CGAL_precondition_code(const bool src_inf = eh->opposite()->vertex()->has_infinite_time();)
  CGAL_precondition_code(const bool tgt_inf = eh->vertex()->has_infinite_time();)
  CGAL_precondition(eh->is_bisector() || (!src_inf && !tgt_inf));

  double sx, sy, tx, ty;
  if(!vertex_position<SS2>(eh->opposite()->vertex(), ctx, sx, sy) ||
     !vertex_position<SS2>(eh->vertex(), ctx, tx, ty))
    return;

  const auto s = make_point<SS2>(sx, sy);
  const auto t = make_point<SS2>(tx, ty);

  if(gs_options.colored_edge(ss2, eh))
    graphics_scene.add_segment(s, t, gs_options.edge_color(ss2, eh));
  else
    graphics_scene.add_segment(s, t);
}

template<class SS2, class GSOptions>
void print_halfedge_labels(const SS2& ss2,
                           typename SS2::Halfedge_const_handle h,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& gs_options,
                           const Drawing_params& ctx)
{
  if(!gs_options.draw_edge(ss2, h))
    return;

  double sx, sy, tx, ty;
  if(!vertex_position<SS2>(h->opposite()->vertex(), ctx, sx, sy) ||
     !vertex_position<SS2>(h->vertex(), ctx, tx, ty))
    return;

  // Place the label closer to the target than to the source, so that the labels
  // of edges sharing a vertex are less likely to overlap, and slightly
  // perturb it so that labels of overlapping edges can be distinguished.
  double x = 0.25 * sx + 0.75 * tx;
  double y = 0.25 * sy + 0.75 * ty;
  const double length = std::sqrt((tx-sx)*(tx-sx) + (ty-sy)*(ty-sy));
  perturb_position(h->id(), (std::max)(0.01 * length, ctx.move_radius), x, y);

  std::stringstream label;
  label << "H" << h->id() << " (V" << h->vertex()->id() << ") ";
  label << "H" << h->opposite()->id() << " (V" << h->opposite()->vertex()->id() << ") ";

  graphics_scene.add_text(make_point<SS2>(x, y), label.str());
}

template <class SS2, class GSOptions>
void compute_vertex(const SS2& ss2, typename SS2::Vertex_const_handle vh,
                    CGAL::Graphics_scene& graphics_scene,
                    const GSOptions& gs_options,
                    const Drawing_params& ctx)
{
  if(!gs_options.draw_vertex(ss2, vh))
    return;

  double x, y;
  if(!vertex_position<SS2>(vh, ctx, x, y))
    return;

  const auto p = make_point<SS2>(x, y);
  if(gs_options.colored_vertex(ss2, vh))
    graphics_scene.add_point(p, gs_options.vertex_color(ss2, vh));
  else
    graphics_scene.add_point(p);
}

template <class SS2, class GSOptions>
void print_vertex_label(const SS2& ss2,
                        typename SS2::Vertex_const_handle vh,
                        CGAL::Graphics_scene& graphics_scene,
                        const GSOptions& gs_options,
                        const Drawing_params& ctx)
{
  // TODO? an option different from draw_vertex to allow to show only some labels ??
  if(!gs_options.draw_vertex(ss2, vh))
  { return; }

  double x, y;
  if(!vertex_position<SS2>(vh, ctx, x, y))
  { return; }

  // Minimum (non-zero) length of the finite edges incident to the vertex
  double min_length = (std::numeric_limits<double>::max)();
  typename SS2::Halfedge_const_handle start = vh->halfedge(), h = start;
  if(!vh->has_infinite_time() && start != typename SS2::Halfedge_const_handle())
  {
    std::size_t guard = 0;
    do
    {
      if(++guard > ss2.size_of_halfedges()) // safety: should never happen in a valid HDS
        break;

      if(!h->opposite()->vertex()->has_infinite_time())
      {
        const auto& p = h->opposite()->vertex()->point();
        const double dx = CGAL::to_double(p.x()) - x;
        const double dy = CGAL::to_double(p.y()) - y;
        const double l = std::sqrt(dx*dx + dy*dy);
        if(l > 0)
          min_length = (std::min)(min_length, l);
      }
      h = h->next()->opposite(); // next incoming halfedge around the vertex
    }
    while(h != start);
  }

  // Slightly perturb the label so that coincident vertices can be distinguished
  double radius = ctx.move_radius;
  if(min_length != (std::numeric_limits<double>::max)())
    radius = (std::min)(radius, 0.01 * min_length);
  perturb_position(vh->id(), radius, x, y);

  std::stringstream label;
  label << "V" << vh->id();
  graphics_scene.add_text(make_point<SS2>(x, y), label.str());
}

template <class SS2, class GSOptions>
void compute_elements(const SS2& ss2,
                      CGAL::Graphics_scene& graphics_scene,
                      const GSOptions& gs_options)
{
  const Drawing_params ctx = compute_drawing_parameters(ss2);

  if(gs_options.are_faces_enabled())
  {
    for(typename SS2::Face_const_iterator it=ss2.faces_begin(); it!=ss2.faces_end(); ++it)
      compute_face(ss2, it, graphics_scene, gs_options);
  }

  if(gs_options.are_edges_enabled())
  {
    for(typename SS2::Halfedge_const_iterator it=ss2.halfedges_begin(); it != ss2.halfedges_end(); ++it)
    {
      if(it->id()<it->opposite()->id())
      {
        compute_edge(ss2, it, graphics_scene, gs_options, ctx);
        print_halfedge_labels(ss2, it, graphics_scene, gs_options, ctx);
      }
    }
  }

  if(gs_options.are_vertices_enabled())
  {
    for(typename SS2::Vertex_const_iterator it=ss2.vertices_begin(); it!=ss2.vertices_end(); ++it)
    {
      compute_vertex(ss2, it, graphics_scene, gs_options, ctx);
      print_vertex_label(ss2, it, graphics_scene, gs_options, ctx);
    }
  }
}

} // namespace draw_function_for_ss2

#define CGAL_SS_TYPE CGAL::Straight_skeleton_2<K>

template <class K, class GSOptions>
void add_to_graphics_scene(const CGAL_SS_TYPE &ass2,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& gs_options)
{
  draw_function_for_ss2::compute_elements(ass2, graphics_scene,
                                          gs_options);
}

template <class K>
void add_to_graphics_scene(const CGAL_SS_TYPE& ass2,
                           CGAL::Graphics_scene& graphics_scene)
{
  Graphics_scene_options<CGAL_SS_TYPE,
                  typename CGAL_SS_TYPE::Vertex_const_handle,
                  typename CGAL_SS_TYPE::Halfedge_const_handle,
                  typename CGAL_SS_TYPE::Face_const_handle>
    drawingFunctor;

  drawingFunctor.colored_edge = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Halfedge_const_handle) -> bool
  { return true; };

  drawingFunctor.colored_vertex = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Vertex_const_handle) -> bool
  { return true; };

  drawingFunctor.colored_face = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Face_const_handle) -> bool
  { return true; };

  drawingFunctor.edge_color = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Halfedge_const_handle eh) -> CGAL::IO::Color
  {
    if(eh->is_bisector())
    {
      // artificial bisectors
      if(eh->face() == eh->opposite()->face())
      { return CGAL::IO::yellow(); }

      // bisectors going to infinity (partial skeletons)
      if(eh->vertex()->has_infinite_time() ||
          eh->opposite()->vertex()->has_infinite_time())
      { return CGAL::IO::orange(); }

      return CGAL::IO::red();
    }

    // Contour edge: the weight is stored in the halfedge incident to a face
    typename CGAL_SS_TYPE::Halfedge_const_handle ch = eh->is_border() ? eh->opposite() : eh;
    if(CGAL::is_zero(ch->weight()))
    { return CGAL::IO::violet(); }

    return CGAL::IO::black();
  };

  drawingFunctor.vertex_color = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Vertex_const_handle vh) -> CGAL::IO::Color
  {
    if(vh->is_split())
      return CGAL::IO::Color(10, 10, 180);  // blue, but not flashy
    else if(vh->has_infinite_time())
      return CGAL::IO::orange();

    return CGAL::IO::Color(10, 180, 10); // green, but not flashy
  };

  drawingFunctor.face_color = []
    (const CGAL_SS_TYPE&, typename CGAL_SS_TYPE::Face_const_handle) -> CGAL::IO::Color
  { return CGAL::IO::Color(150, 170, 170); };

  add_to_graphics_scene(ass2, graphics_scene, drawingFunctor);
}

// Specialization of draw function.
template <class K, class GSOptions>
void draw(const CGAL_SS_TYPE &ass2, const GSOptions &gs_options,
          const char *title="Straight Skeleton Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(ass2, buffer, gs_options);
  draw_graphics_scene(buffer, title);
}

template <class K>
void draw(const CGAL_SS_TYPE &ass2,
          const char *title="Straight Skeleton Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(ass2, buffer);
  draw_graphics_scene(buffer, title);
}

#undef CGAL_SS_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_SS2_H