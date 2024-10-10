// Copyright(c) 2019  Foundation for Research and Technology-Hellas (Greece).
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

#ifndef CGAL_DRAW_VORONOI_DIAGRAM_2_H
#define CGAL_DRAW_VORONOI_DIAGRAM_2_H

#include <CGAL/license/Voronoi_diagram_2.h>
#include <CGAL/Qt/Basic_viewer.h>
#include <CGAL/Graphics_scene.h>
#include <CGAL/Graphics_scene_options.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/Voronoi_diagram_2/Face.h>
#include <CGAL/Voronoi_diagram_2/Handle_adaptor.h>
#include <CGAL/Voronoi_diagram_2/Vertex.h>
#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Accessor.h>

namespace CGAL {

// We need a specific graphics scene option for voronoi2 in order to allow
// to differentiate voronoi and dual vertices, and to manage rays.
template <typename DS,
          typename vertex_handle,
          typename edge_handle,
          typename face_handle>
struct Graphics_scene_options_voronoi_2 :
    public CGAL::Graphics_scene_options<DS, vertex_handle, edge_handle, face_handle>
{
  Graphics_scene_options_voronoi_2() : m_dual_vertex_color(50, 100, 180),
                                       m_ray_color(100, 0, 0),
                                       m_bisector_color(0, 100, 0),
                                       m_draw_voronoi_vertices(true),
                                       m_draw_dual_vertices(true)
  {}

  const CGAL::IO::Color& dual_vertex_color() const
  { return m_dual_vertex_color; }
  const CGAL::IO::Color& ray_color() const
  { return m_ray_color; }
  const CGAL::IO::Color& bisector_color() const
  { return m_bisector_color; }

  void dual_vertex_color(const CGAL::IO::Color& c)
  { m_dual_vertex_color=c; }
  void ray_color(const CGAL::IO::Color& c)
  { m_ray_color=c; }
  void bisector_color(const CGAL::IO::Color& c)
  { m_bisector_color=c; }

  void draw_voronoi_vertices(bool b) { m_draw_voronoi_vertices=b; }
  bool draw_voronoi_vertices() const { return m_draw_voronoi_vertices; }
  void toggle_draw_voronoi_vertices() { m_draw_voronoi_vertices=!m_draw_voronoi_vertices; }

  void draw_dual_vertices(bool b) { m_draw_dual_vertices=b; }
  bool draw_dual_vertices() const { return m_draw_dual_vertices; }
  void toggle_draw_dual_vertices() { m_draw_dual_vertices=!m_draw_dual_vertices; }

protected:
  CGAL::IO::Color m_dual_vertex_color;
  CGAL::IO::Color m_ray_color;
  CGAL::IO::Color m_bisector_color;
  bool m_draw_voronoi_vertices;
  bool m_draw_dual_vertices;
};

namespace draw_function_for_v2
{

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3  Local_point;
typedef Local_kernel::Vector_3 Local_vector;

template <class V2, class GSOptions>
void compute_vertex(const V2& v2,
                    typename V2::Vertex_iterator vh,
                    CGAL::Graphics_scene& graphics_scene,
                    const GSOptions& gs_options)
{
  if(!gs_options.draw_vertex(v2, vh))
  { return; }

  if(gs_options.colored_vertex(v2, vh))
  { graphics_scene.add_point(vh->point(), gs_options.vertex_color(v2, vh)); }
  else
  { graphics_scene.add_point(vh->point()); }
}

template <class V2, class GSOptions>
void compute_dual_vertex(const V2& /*v2*/,
                         typename V2::Delaunay_graph::Finite_vertices_iterator vi,
                         CGAL::Graphics_scene &graphics_scene,
                         const GSOptions& gs_options)
{ graphics_scene.add_point(vi->point(), gs_options.dual_vertex_color()); }

template <class V2, class GSOptions>
void add_segments_and_update_bounding_box(const V2& v2,
                                          typename V2::Halfedge_iterator he,
                                          CGAL::Graphics_scene& graphics_scene,
                                          GSOptions& gs_options)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef typename V2::Delaunay_vertex_handle Delaunay_vertex_const_handle;

  if (he->is_segment())
  {
    if(gs_options.draw_edge(v2, he))
    {
      if(gs_options.colored_edge(v2, he))
      {
        graphics_scene.add_segment(he->source()->point(), he->target()->point(),
                                   gs_options.edge_color(v2, he));
      }
      else
      {
        graphics_scene.add_segment(he->source()->point(), he->target()->point());
      }
    }
  }
  else
  {
    Delaunay_vertex_const_handle v1 = he->up();
    Delaunay_vertex_const_handle v2 = he->down();

    Kernel::Vector_2 direction(v1->point().y() - v2->point().y(),
                               v2->point().x() - v1->point().x());
    if (he->is_ray())
    {
      Kernel::Point_2 end_point;
      if (he->has_source())
      {
        end_point = he->source()->point();

        // update_bounding_box_for_ray(end_point, direction);
        Local_point lp = graphics_scene.get_local_point(end_point);
        Local_vector lv = graphics_scene.get_local_vector(direction);
        CGAL::Bbox_3 b = (lp + lv).bbox();
        graphics_scene.update_bounding_box(b);
      }
    }
    else if (he->is_bisector())
    {
      Kernel::Point_2 pointOnLine((v1->point().x() + v2->point().x()) / 2,
                                  (v1->point().y() + v2->point().y()) / 2);
      Kernel::Vector_2 perpendicularDirection(
          v2->point().x() - v1->point().x(),
          v2->point().y() - v1->point().y());

      // update_bounding_box_for_line(pointOnLine, direction,
      //                               perpendicularDirection);
      Local_point lp = graphics_scene.get_local_point(pointOnLine);
      Local_vector lv = graphics_scene.get_local_vector(direction);
      Local_vector lpv = graphics_scene.get_local_vector(perpendicularDirection);

      CGAL::Bbox_3 b = lp.bbox() + (lp + lv).bbox() + (lp + lpv).bbox();
      graphics_scene.update_bounding_box(b);
    }
  }
}

template <class V2>
Local_kernel::Point_2 get_second_point(typename V2::Halfedge_handle ray,
                                       const CGAL::Bbox_3 & m_bounding_box)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef typename V2::Delaunay_vertex_handle Delaunay_vertex_const_handle;

  Delaunay_vertex_const_handle v1 = ray->up();
  Delaunay_vertex_const_handle v2 = ray->down();

  // calculate direction of ray and its inverse
  Kernel::Vector_2 v(v1->point().y() - v2->point().y(),
                     v2->point().x() - v1->point().x());
  Local_kernel::Vector_2 inv(1 / v.x(), 1 / v.y());

  // origin of the ray
  Kernel::Point_2 p;
  if (ray->has_source())
  { p = ray->source()->point(); }
  else
  { p = ray->target()->point(); }

  // get the bounding box of the viewer
  Local_kernel::Vector_2 boundsMin(m_bounding_box.xmin(),
                                   m_bounding_box.zmin());
  Local_kernel::Vector_2 boundsMax(m_bounding_box.xmax(),
                                   m_bounding_box.zmax());
  // calculate intersection
  double txmax, txmin, tymax, tymin;

  if (inv.x() >= 0)
  {
    txmax = (boundsMax.x() - p.x()) * inv.x();
    txmin = (boundsMin.x() - p.x()) * inv.x();
  }
  else
  {
    txmax = (boundsMin.x() - p.x()) * inv.x();
    txmin = (boundsMax.x() - p.x()) * inv.x();
  }

  if (inv.y() >= 0)
  {
    tymax = (boundsMax.y() - p.y()) * inv.y();
    tymin = (boundsMin.y() - p.y()) * inv.y();
  }
  else
  {
    tymax = (boundsMin.y() - p.y()) * inv.y();
    tymin = (boundsMax.y() - p.y()) * inv.y();
  }

  if (tymin > txmin)
    txmin = tymin;
  if (tymax < txmax)
    txmax = tymax;

  Local_kernel::Point_2 p1;
  if (v.x() == 0)
  { p1 = Local_kernel::Point_2(p.x(), p.y() + tymax * v.y()); }
  else if (v.y() == 0)
  { p1 = Local_kernel::Point_2(p.x() + txmax * v.x(), p.y()); }
  else
  { p1 = Local_kernel::Point_2(p.x() + txmax * v.x(), p.y() + tymax * v.y()); }
  return p1;
}

// Halfedge_const_handle
template <class V2, class GSOptions>
void compute_rays_and_bisectors(const V2&,
                                typename V2::Halfedge_iterator he,
                                CGAL::Graphics_scene& graphics_scene,
                                const GSOptions& gs_options)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef typename V2::Delaunay_vertex_handle Delaunay_vertex_const_handle;

  Delaunay_vertex_const_handle v1 = he->up();
  Delaunay_vertex_const_handle v2 = he->down();

  Kernel::Vector_2 direction(v1->point().y() - v2->point().y(),
                             v2->point().x() - v1->point().x());
  if (he->is_ray())
  {
    if (he->has_source())
    {
      // add_ray_segment(he->source()->point(), get_second_point(he, graphics_scene.get_bounding_box()));
      graphics_scene.add_ray(he->source()->point(), direction, gs_options.ray_color());
    }
  }
  else if (he->is_bisector())
  {
    Kernel::Point_2 pointOnLine((v1->point().x() + v2->point().x()) / 2,
                                (v1->point().y() + v2->point().y()) / 2);
    graphics_scene.add_line(pointOnLine, direction, gs_options.bisector_color());
  }
}

template <class V2, class GSOptions>
void compute_face(const V2& v2,
                  typename V2::Face_iterator fh,
                  CGAL::Graphics_scene& graphics_scene,
                  const GSOptions& m_gs_options)
{
  if(fh->is_unbounded() || !m_gs_options.draw_face(v2, fh))
  { return; }

  if(m_gs_options.colored_face(v2, fh))
  { graphics_scene.face_begin(m_gs_options.face_color(v2, fh)); }
  else { graphics_scene.face_begin(); }

  typename V2::Ccb_halfedge_circulator ec_start=fh->ccb();
  typename V2::Ccb_halfedge_circulator ec=ec_start;
  do
  {
    graphics_scene.add_point_in_face(ec->source()->point());
  }
  while (++ec!=ec_start);
  graphics_scene.face_end();

  // Test: for unbounded faces (??)
  //      else {
  //            do{
  //                if( ec->has_source() ){
  //                    add_point_in_face(ec->source()->point());
  //                }
  //                else{
  //                    add_point_in_face(get_second_point(ec->twin(), graphics_scene.get_bounding_box()));
  //                }
  //            } while(++ec != ec_start);
  //        }
}

template <class V2, class GSOptions>
void compute_elements(const V2& v2,
                      CGAL::Graphics_scene& graphics_scene,
                      const GSOptions& gs_options)
{
  if(gs_options.are_vertices_enabled())
  {
    // Draw the voronoi vertices
    if (gs_options.draw_voronoi_vertices())
    {
      for (typename V2::Vertex_iterator it=v2.vertices_begin();
           it!=v2.vertices_end(); ++it)
      { compute_vertex(v2, it, graphics_scene, gs_options); }
    }

    // Draw the dual vertices
    if (gs_options.draw_dual_vertices())
    {
      for (typename V2::Delaunay_graph::Finite_vertices_iterator
             it=v2.dual().finite_vertices_begin();
           it!=v2.dual().finite_vertices_end(); ++it)
      { compute_dual_vertex(v2, it, graphics_scene, gs_options); }
    }
  }

  if(gs_options.are_edges_enabled())
  {
    // Add segments and update bounding box
    for (typename V2::Halfedge_iterator it=v2.halfedges_begin();
         it!=v2.halfedges_end(); ++it)
    { add_segments_and_update_bounding_box(v2, it,
                                           graphics_scene, gs_options); }
  }

  for (typename V2::Halfedge_iterator it=v2.halfedges_begin();
       it!=v2.halfedges_end(); ++it)
  { compute_rays_and_bisectors(v2, it, graphics_scene, gs_options); }

  if (gs_options.are_faces_enabled())
  {
    for (typename V2::Face_iterator it=v2.faces_begin(); it!=v2.faces_end(); ++it)
    { compute_face(v2, it, graphics_scene, gs_options); }
  }
}

} // namespace draw_function_for_v2

#define CGAL_VORONOI_TYPE CGAL::Voronoi_diagram_2 <DG, AT, AP>

template <class DG, class AT, class AP, class GSOptions>
void add_to_graphics_scene(const CGAL_VORONOI_TYPE &v2,
                           CGAL::Graphics_scene& graphics_scene,
                           const GSOptions& m_gs_options)
{
  draw_function_for_v2::compute_elements(v2, graphics_scene, m_gs_options);
}

template <class DG, class AT, class AP>
void add_to_graphics_scene(const CGAL_VORONOI_TYPE& v2,
                           CGAL::Graphics_scene& graphics_scene)
{
  // Default graphics view options.
  CGAL::Graphics_scene_options_voronoi_2<CGAL_VORONOI_TYPE,
                                typename CGAL_VORONOI_TYPE::Vertex_iterator,
                                typename CGAL_VORONOI_TYPE::Halfedge_iterator,
                                typename CGAL_VORONOI_TYPE::Face_iterator>
    gs_options;

  add_to_graphics_scene(v2, graphics_scene, gs_options);
}

#ifdef CGAL_USE_BASIC_VIEWER

// Specialization of draw function.
template<class DG, class AT, class AP, class GSOptions>
void draw(const CGAL_VORONOI_TYPE& av2,
          GSOptions& gs_options,
          const char *title="2D Voronoi Diagram Basic Viewer")
{
  CGAL::Graphics_scene buffer;
  add_to_graphics_scene(av2, buffer, gs_options);

  CGAL::Qt::QApplication_and_basic_viewer app(buffer, title);
  if(app)
  {
    // Here we define the std::function to capture key pressed.
    app.basic_viewer().on_key_pressed=
      [&av2, &buffer, &gs_options] (QKeyEvent* e, CGAL::Qt::Basic_viewer* basic_viewer) -> bool
      {
        const ::Qt::KeyboardModifiers modifiers = e->modifiers();
        if ((e->key() == ::Qt::Key_R) && (modifiers == ::Qt::NoButton))
        {
          basic_viewer->toggle_draw_rays();
          basic_viewer->displayMessage
            (QString("Draw rays=%1.").arg(basic_viewer->draw_rays()?"true":"false"));

          basic_viewer->redraw();
        }
        else if ((e->key() == ::Qt::Key_V) && (modifiers == ::Qt::ShiftModifier))
        {
          gs_options.toggle_draw_voronoi_vertices();
          basic_viewer->displayMessage
            (QString("Voronoi vertices=%1.").
             arg(gs_options.draw_voronoi_vertices()?"true":"false"));

          buffer.clear();
          draw_function_for_v2::compute_elements(av2, buffer, gs_options);
          basic_viewer->redraw();
        }
        else if ((e->key() == ::Qt::Key_D) && (modifiers == ::Qt::NoButton))
        {
          gs_options.toggle_draw_dual_vertices();
          basic_viewer->displayMessage(QString("Dual vertices=%1.").
                                       arg(gs_options.draw_dual_vertices()?"true":"false"));

          buffer.clear();
          draw_function_for_v2::compute_elements(av2, buffer, gs_options);
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
    app.basic_viewer().setKeyDescription(::Qt::Key_R, "Toggles rays display");
    app.basic_viewer().setKeyDescription(::Qt::Key_D, "Toggles dual vertices display");
    app.basic_viewer().setKeyDescription(::Qt::ShiftModifier, ::Qt::Key_V, "Toggles voronoi vertices display");

    // Then we run the app
    app.run();
  }
}

template<class DG, class AT, class AP>
void draw(const CGAL_VORONOI_TYPE& av2,
          const char *title="2D Voronoi Diagram Basic Viewer")
{
  CGAL::Graphics_scene_options_voronoi_2<CGAL_VORONOI_TYPE,
                                typename CGAL_VORONOI_TYPE::Vertex_iterator,
                                typename CGAL_VORONOI_TYPE::Halfedge_iterator,
                                typename CGAL_VORONOI_TYPE::Face_iterator>
    gs_options;
  draw(av2, gs_options, title);
}

#endif // CGAL_USE_BASIC_VIEWER

#undef CGAL_VORONOI_TYPE

} // End namespace CGAL

#endif // CGAL_DRAW_VORONOI_DIAGRAM_2_H
