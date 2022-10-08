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

#include <CGAL/Graphic_buffer.h>
#include <CGAL/Drawing_functor.h>
#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Voronoi_diagram_2.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Qt/init_ogl_context.h>
#include <CGAL/Random.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Voronoi_diagram_2/Face.h>
#include <CGAL/Voronoi_diagram_2/Handle_adaptor.h>
#include <CGAL/Voronoi_diagram_2/Vertex.h>
#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Accessor.h>

namespace CGAL {

namespace draw_function_for_v2
{

template <typename DS,
          typename vertex_handle,
          typename edge_handle,
          typename face_handle>
struct Drawing_functor_voronoi :
    public CGAL::Drawing_functor<DS, vertex_handle, edge_handle, face_handle>
{
  Drawing_functor_voronoi() : m_draw_rays(true),
                              m_draw_voronoi_vertices(true),
                              m_draw_dual_vertices(true),
                              m_nofaces(false)
  {

    draw_voronoi_vertices=[](const DS &, vertex_handle)->bool { return true; };
    draw_dual_vertices=[](const DS &, vertex_handle)->bool { return true; };
    draw_rays=[](const DS &, vertex_handle)->bool { return true; };

  }

  std::function<bool(const DS&, vertex_handle)> draw_voronoi_vertices;
  std::function<bool(const DS&, vertex_handle)> draw_dual_vertices;
  std::function<bool(const DS&, vertex_handle)> draw_rays;


  void disable_voronoi_vertices() { m_draw_voronoi_vertices=false; }
  void enable_voronoi_vertices() { m_draw_voronoi_vertices=true; }
  bool are_voronoi_vertices_enabled() const { return m_draw_voronoi_vertices; }

  void disable_dual_vertices() { m_draw_dual_vertices=false; }
  void enable_dual_vertices() { m_draw_dual_vertices=true; }
  bool are_dual_vertices_enabled() const { return m_draw_dual_vertices; }

  void disable_rays() { m_draw_rays=false; }
  void enable_rays() { m_draw_rays=true; }
  bool are_rays_enabled() const { return m_draw_rays; }


  void disable_nofaces() { m_nofaces=false; }
  void enable_nofaces() { m_nofaces=true; }
  bool are_nofaces_enabled() const { return m_nofaces; }

protected:
  bool m_draw_rays;
  bool m_draw_voronoi_vertices;
  bool m_draw_dual_vertices;
  bool m_nofaces;
};

typedef CGAL::Exact_predicates_inexact_constructions_kernel Local_kernel;
typedef Local_kernel::Point_3  Local_point;
typedef Local_kernel::Vector_3 Local_vector;

template <typename BufferType = float, class V2>
void compute_vertex(typename V2::Vertex_iterator vh,
                    CGAL::Graphic_buffer<BufferType> &graphic_buffer) {
  graphic_buffer.add_point(vh->point()); 
}

template <typename BufferType = float, class V2>
void compute_dual_vertex(typename V2::Delaunay_graph::Finite_vertices_iterator vi,
                        CGAL::Graphic_buffer<BufferType> &graphic_buffer) {
  graphic_buffer.add_point(vi->point(), CGAL::IO::Color(50, 100, 180));
}

template <typename BufferType = float, class V2>
void add_segments_and_update_bounding_box(typename V2::Halfedge_handle he,
                                          CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef typename V2::Delaunay_vertex_handle Delaunay_vertex_const_handle;

  if (he->is_segment()) {
    graphic_buffer.add_segment(he->source()->point(), he->target()->point());
  } else {
    Delaunay_vertex_const_handle v1 = he->up();
    Delaunay_vertex_const_handle v2 = he->down();

    Kernel::Vector_2 direction(v1->point().y() - v2->point().y(),
                                v2->point().x() - v1->point().x());
    if (he->is_ray()) {
      Kernel::Point_2 end_point;
      if (he->has_source()) {
        end_point = he->source()->point();

    // update_bounding_box_for_ray(end_point, direction);

    // update_bounding_box_for_ray
    Local_point lp = Basic_viewer_qt<>::get_local_point(end_point);
    Local_vector lv = Basic_viewer_qt<>::get_local_vector(direction);
    CGAL::Bbox_3 b = (lp + lv).bbox();
    graphic_buffer.update_bounding_box(b);
  

      }
    } else if (he->is_bisector()) {
      Kernel::Point_2 pointOnLine((v1->point().x() + v2->point().x()) / 2,
                                  (v1->point().y() + v2->point().y()) / 2);
      Kernel::Vector_2 perpendicularDirection(
          v2->point().x() - v1->point().x(),
          v2->point().y() - v1->point().y());

      // update_bounding_box_for_line(pointOnLine, direction,
      //                               perpendicularDirection);

      // update_bounding_box_for_line
      Local_point lp = Basic_viewer_qt<>::get_local_point(pointOnLine);
      Local_vector lv = Basic_viewer_qt<>::get_local_vector(direction);
      Local_vector lpv = Basic_viewer_qt<>::get_local_vector(perpendicularDirection);

      CGAL::Bbox_3 b = lp.bbox() + (lp + lv).bbox() + (lp + lpv).bbox();
      graphic_buffer.update_bounding_box(b);

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
  if (ray->has_source()) {
    p = ray->source()->point();
  } else {
    p = ray->target()->point();
  }

  // get the bounding box of the viewer
  Local_kernel::Vector_2 boundsMin(m_bounding_box.xmin(),
                                    m_bounding_box.zmin());
  Local_kernel::Vector_2 boundsMax(m_bounding_box.xmax(),
                                    m_bounding_box.zmax());
  // calculate intersection
  double txmax, txmin, tymax, tymin;

  if (inv.x() >= 0) {
    txmax = (boundsMax.x() - p.x()) * inv.x();
    txmin = (boundsMin.x() - p.x()) * inv.x();
  } else {
    txmax = (boundsMin.x() - p.x()) * inv.x();
    txmin = (boundsMax.x() - p.x()) * inv.x();
  }

  if (inv.y() >= 0) {
    tymax = (boundsMax.y() - p.y()) * inv.y();
    tymin = (boundsMin.y() - p.y()) * inv.y();
  } else {
    tymax = (boundsMin.y() - p.y()) * inv.y();
    tymin = (boundsMax.y() - p.y()) * inv.y();
  }

  if (tymin > txmin)
    txmin = tymin;
  if (tymax < txmax)
    txmax = tymax;

  Local_kernel::Point_2 p1;
  if (v.x() == 0) {
    p1 = Local_kernel::Point_2(p.x(), p.y() + tymax * v.y());
  } else if (v.y() == 0) {
    p1 = Local_kernel::Point_2(p.x() + txmax * v.x(), p.y());
  } else {
    p1 = Local_kernel::Point_2(p.x() + txmax * v.x(), p.y() + tymax * v.y());
  }
  return p1;
}

// Halfedge_const_handle
template <typename BufferType = float, class V2>
void compute_rays_and_bisectors(typename V2::Halfedge_iterator he,
                              CGAL::Graphic_buffer<BufferType> &graphic_buffer)
{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  typedef typename V2::Delaunay_vertex_handle Delaunay_vertex_const_handle;

  Delaunay_vertex_const_handle v1 = he->up();
  Delaunay_vertex_const_handle v2 = he->down();

  Kernel::Vector_2 direction(v1->point().y() - v2->point().y(),
                              v2->point().x() - v1->point().x());
  if (he->is_ray()) {
    if (he->has_source()) {
      // add_ray_segment(he->source()->point(), get_second_point(he, graphic_buffer.get_bounding_box()));
      graphic_buffer.add_ray(he->source()->point(), direction, CGAL::IO::Color(100, 0, 0));
    }
  } else if (he->is_bisector()) {
    Kernel::Point_2 pointOnLine((v1->point().x() + v2->point().x()) / 2,
                                (v1->point().y() + v2->point().y()) / 2);
    graphic_buffer.add_line(pointOnLine, direction);
  }
}

template <typename BufferType = float, class V2, class DrawingFunctor>
void compute_face(typename V2::Face_iterator fh,
                  const V2& v2,
                  CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                  const DrawingFunctor &m_drawing_functor)
{
  typedef typename V2::Ccb_halfedge_circulator  Ccb_halfedge_circulator;

  if(m_drawing_functor.colored_face(v2, fh)) {

  CGAL::IO::Color c = m_drawing_functor.face_color(v2, fh);

  Ccb_halfedge_circulator ec_start = fh->ccb();
  Ccb_halfedge_circulator ec = ec_start;

  if (!fh->is_unbounded()) {
    graphic_buffer.face_begin(c);
    do {
      graphic_buffer.add_point_in_face(ec->source()->point());
    } while (++ec != ec_start);
    graphic_buffer.face_end();
  }

  }
  // Test: for unbounded faces
  //      else {
  //            do{
  //                if( ec->has_source() ){
  //                    add_point_in_face(ec->source()->point());
  //                }
  //                else{
  //                    add_point_in_face(get_second_point(ec->twin(), graphic_buffer.get_bounding_box()));
  //                }
  //            } while(++ec != ec_start);
  //        }
}

template <typename BufferType = float, class V2, class DrawingFunctor>
void compute_elements(const V2& v2, CGAL::Graphic_buffer<BufferType> &graphic_buffer, 
                      const DrawingFunctor &m_drawing_functor)
{
  typedef typename V2::Delaunay_graph::Finite_vertices_iterator Dual_vertices_iterator;

  // Draw the voronoi vertices
  if (m_drawing_functor.are_voronoi_vertices_enabled()) {
    for (typename V2::Vertex_iterator it = v2.vertices_begin();
          it != v2.vertices_end(); ++it) {
      compute_vertex<BufferType, V2>(it, graphic_buffer);
    }
  }

  // Draw the dual vertices
  if (m_drawing_functor.are_dual_vertices_enabled()) {
    for (Dual_vertices_iterator it = v2.dual().finite_vertices_begin();
          it != v2.dual().finite_vertices_end(); ++it) {
      compute_dual_vertex<BufferType, V2>(it, graphic_buffer);
    }
  }

  // Add segments and update bounding box
  for (typename V2::Halfedge_iterator it = v2.halfedges_begin();
        it != v2.halfedges_end(); ++it) {
    add_segments_and_update_bounding_box<BufferType, V2>(it, graphic_buffer);
  }

  if(m_drawing_functor.are_rays_enabled()) {
    for (typename V2::Halfedge_iterator it = v2.halfedges_begin();
          it != v2.halfedges_end(); ++it) {
      compute_rays_and_bisectors<BufferType, V2>(it, graphic_buffer);
    }
  }

  if (!m_drawing_functor.are_nofaces_enabled()) {
    for (typename V2::Face_iterator it = v2.faces_begin();
          it != v2.faces_end(); ++it) {
      compute_face(it, v2, graphic_buffer, m_drawing_functor);
    }
  }
}

} // namespace draw_function_for_v2

template <typename BufferType = float, class V2, class DrawingFunctor>
void add_in_graphic_buffer(const V2 &v2, CGAL::Graphic_buffer<BufferType> &graphic_buffer,
                             const DrawingFunctor &m_drawing_functor)
{
  draw_function_for_v2::compute_elements(v2, graphic_buffer, m_drawing_functor);
}

template <typename BufferType = float, class V2>
void add_in_graphic_buffer(const V2 &v2, CGAL::Graphic_buffer<BufferType> &graphic_buffer, 
                           bool m_nofaces = false, bool m_draw_voronoi_vertices = true,
                           bool m_draw_dual_vertices = true ) {

  // Default functor; user can add his own functor.
  CGAL::draw_function_for_v2::Drawing_functor_voronoi<V2, typename V2::Vertex_iterator,
                typename V2::Halfedge_iterator,
                typename V2::Face_iterator>
    drawing_functor;

  drawing_functor.colored_face = [](const V2&,
                      typename V2::Face_iterator fh) -> bool
  { return true; };


  drawing_functor.face_color =  [] (const V2& alcc,
                           typename V2::Face_iterator fh) -> CGAL::IO::Color
  {
    return CGAL::IO::Color(73, 250, 117);
  };

  add_in_graphic_buffer(v2, graphic_buffer, drawing_functor);
}

// TODO: I'll add a description after idea verification
// Add custom key description (see keyPressEvent)
// setKeyDescription(::Qt::Key_R, "Toggles rays display");
// setKeyDescription(::Qt::Key_D, "Toggles dual vertices display");
// setKeyDescription(::Qt::Key_V, "Toggles voronoi vertices display");


// TODO: Expected solution, I need to add reference to drawing_functor, OR, pass m_draw_rays, m_draw_voronoi_vertices, m_draw_dual_vertices
std::function<void(QKeyEvent *, CGAL::Basic_viewer_qt<float> *)> VoronoiKeyPressEvent = [] (QKeyEvent *e, CGAL::Basic_viewer_qt<float> *_this)
{
  /// [Keypress]
  const ::Qt::KeyboardModifiers modifiers = e->modifiers();
  if ((e->key() == ::Qt::Key_R) && (modifiers == ::Qt::NoButton)) {
    // m_draw_rays = !m_draw_rays;
    std::cout << "R Pressed\n";

    _this->displayMessage(
        QString("Draw rays=%1."));//.arg(m_draw_rays ? "true" : "false"));
    // _this.update();
  } else if ((e->key() == ::Qt::Key_V) && (modifiers == ::Qt::NoButton)) {

    std::cout << "V Pressed\n";

    // m_draw_voronoi_vertices = !m_draw_voronoi_vertices;

    _this->displayMessage(
        QString("Voronoi vertices=%1."));//.arg(m_draw_voronoi_vertices? "true" : "false"));

    // TODO: What are the expected args? I think I need a kind of function wrapper to wrap drawing_functor, mesh, and buffer.
    // draw_function_for_v2::compute_elements();

    // _this->redraw();
  } else if ((e->key() == ::Qt::Key_D) && (modifiers == ::Qt::NoButton)) {

    std::cout << "D Pressed\n";

    // m_draw_dual_vertices = !m_draw_dual_vertices;

    _this->displayMessage(QString("Dual vertices=%1."));
                        //.arg(m_draw_dual_vertices ? "true" : "false"));

    // TODO: What are the expected args? I think I need a kind of function wrapper to wrap drawing_functor, mesh, and buffer.
    // draw_function_for_v2::compute_elements();

    // _this->redraw();
  } else {
    // Call the base method to process others/classicals key
    // Base::keyPressEvent(e);
  }
  /// [Keypress]
};

// Specialization of draw function.
#define CGAL_VORONOI_TYPE CGAL::Voronoi_diagram_2 <DG, AT, AP>

template<class DG,
         class AT,
         class AP, typename BufferType = float, class DrawingFunctor>
void draw(const CGAL_VORONOI_TYPE &av2,
          const DrawingFunctor &drawing_functor)
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(av2, buffer, drawing_functor);
  draw_buffer(buffer);
}

template<class DG,
         class AT,
         class AP, typename BufferType = float>
void draw(const CGAL_VORONOI_TYPE &av2,
          const char *title="2D Voronoi Diagram Basic Viewer")
{
  CGAL::Graphic_buffer<BufferType> buffer;
  add_in_graphic_buffer(av2, buffer);
  // draw_buffer(buffer);

  // Test
  draw_buffer(buffer, VoronoiKeyPressEvent);

}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_VORONOI_DIAGRAM_2_H
