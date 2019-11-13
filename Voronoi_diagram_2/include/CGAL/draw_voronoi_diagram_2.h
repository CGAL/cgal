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

#ifndef CGAL_DRAW_VORONOI_DIAGRAM_2_H
#define CGAL_DRAW_VORONOI_DIAGRAM_2_H

#include <CGAL/Qt/Basic_viewer_qt.h>
#include <CGAL/license/Voronoi_diagram_2.h>

#ifdef CGAL_USE_BASIC_VIEWER

#include <CGAL/Random.h>
#include <CGAL/Triangulation_utils_2.h>
#include <CGAL/Voronoi_diagram_2/Face.h>
#include <CGAL/Voronoi_diagram_2/Handle_adaptor.h>
#include <CGAL/Voronoi_diagram_2/Vertex.h>
#include <CGAL/Voronoi_diagram_2/basic.h>
#include <CGAL/Voronoi_diagram_2/Accessor.h>

namespace CGAL {

// Default color functor; user can change it to have its own face color
struct DefaultColorFunctorV2
{
  template <typename V2>
  static CGAL::Color run(const V2 &, const typename V2::Face_iterator /*fh*/) {
    //CGAL::Random random((unsigned int)(std::size_t)(&*fh));
    //return get_random_color(random);
    return CGAL::Color(73, 250, 117);
  }
};

// Viewer for Voronoi diagram
template <class V2, class ColorFunctor>
class SimpleVoronoiDiagram2ViewerQt : public Basic_viewer_qt
{
  typedef Basic_viewer_qt                                       Base;
  typedef typename V2::Vertex_iterator                          Vertex_const_handle;
  typedef typename V2::Delaunay_vertex_handle                   Delaunay_vertex_const_handle;
  typedef typename V2::Delaunay_graph::Finite_vertices_iterator Dual_vertices_iterator;

  typedef typename V2::Halfedge_iterator                        Halfedge_const_handle;
  typedef typename V2::Ccb_halfedge_circulator                  Ccb_halfedge_circulator;
  typedef typename V2::Halfedge_handle                          Halfedge_handle;

  typedef typename V2::Face_iterator                            Face_const_handle;

  typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

public:
  /// Construct the viewer.
  /// @param av2 the voronoi diagram to view
  /// @param title the title of the window
  /// @param anofaces if true, do not draw faces (faces are not computed; this
  /// can be useful for very big object where this time could be long)
  SimpleVoronoiDiagram2ViewerQt(QWidget *parent, const V2 &av2,
                                const char *title = "Basic Voronoi Viewer",
                                bool anofaces = false,
                                bool draw_voronoi_vertices = true,
                                bool draw_delaunay_vertices = true,
                                const ColorFunctor &fcolor = ColorFunctor())
      : // First draw: vertices; half-edges; faces; multi-color; no inverse
        // normal
        Base(parent, title, true, true, true, false, false, true, true),
        v2(av2), m_nofaces(anofaces),
        m_draw_voronoi_vertices(draw_voronoi_vertices),
        m_draw_dual_vertices(draw_delaunay_vertices), m_fcolor(fcolor) {
    // Add custom key description (see keyPressEvent)
    setKeyDescription(::Qt::Key_R, "Toggles rays display");
    setKeyDescription(::Qt::Key_D, "Toggles dual vertices display");
    setKeyDescription(::Qt::Key_V, "Toggles voronoi vertices display");

    compute_elements();
  }

protected:

  void compute_vertex(Vertex_const_handle vh) { add_point(vh->point()); }

  void compute_dual_vertex(Dual_vertices_iterator vi)
  {
    add_point(vi->point(), CGAL::Color(50, 100, 180));
  }

  void add_segments_and_update_bounding_box(Halfedge_handle he)
  {
    if (he->is_segment()) {
      add_segment(he->source()->point(), he->target()->point());
    } else {
      Delaunay_vertex_const_handle v1 = he->up();
      Delaunay_vertex_const_handle v2 = he->down();

      Kernel::Vector_2 direction(v1->point().y() - v2->point().y(),
                                 v2->point().x() - v1->point().x());
      if (he->is_ray()) {
        Kernel::Point_2 end_point;
        if (he->has_source()) {
          end_point = he->source()->point();
          update_bounding_box_for_ray(end_point, direction);
        }
      } else if (he->is_bisector()) {
        Kernel::Point_2 pointOnLine((v1->point().x() + v2->point().x()) / 2,
                                    (v1->point().y() + v2->point().y()) / 2);
        Kernel::Vector_2 perpendicularDirection(
            v2->point().x() - v1->point().x(),
            v2->point().y() - v1->point().y());
        update_bounding_box_for_line(pointOnLine, direction,
                                     perpendicularDirection);
      }
    }
  }

  Local_kernel::Point_2 get_second_point(Halfedge_handle ray)
  {
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

  void compute_rays_and_bisectors(Halfedge_const_handle he)
  {
    Delaunay_vertex_const_handle v1 = he->up();
    Delaunay_vertex_const_handle v2 = he->down();

    Kernel::Vector_2 direction(v1->point().y() - v2->point().y(),
                               v2->point().x() - v1->point().x());
    if (he->is_ray()) {
      if (he->has_source()) {
        // add_ray_segment(he->source()->point(), get_second_point(he));
        add_ray(he->source()->point(), direction, CGAL::Color(100, 0, 0));
      }
    } else if (he->is_bisector()) {
      Kernel::Point_2 pointOnLine((v1->point().x() + v2->point().x()) / 2,
                                  (v1->point().y() + v2->point().y()) / 2);
      add_line(pointOnLine, direction);
    }
  }

  void compute_face(Face_const_handle fh)
  {
    CGAL::Color c = m_fcolor.run(v2, fh);

    Ccb_halfedge_circulator ec_start = fh->ccb();
    Ccb_halfedge_circulator ec = ec_start;

    if (!fh->is_unbounded()) {
      face_begin(c);
      do {
        add_point_in_face(ec->source()->point());
      } while (++ec != ec_start);
      face_end();
    }
    // Test: for unbounded faces
    //      else {
    //            do{
    //                if( ec->has_source() ){
    //                    add_point_in_face(ec->source()->point());
    //                }
    //                else{
    //                    add_point_in_face(get_second_point(ec->twin()));
    //                }
    //            } while(++ec != ec_start);
    //        }
  }

  void compute_elements()
  {
    clear();

    // Draw the voronoi vertices
    if (m_draw_voronoi_vertices) {
      for (typename V2::Vertex_iterator it = v2.vertices_begin();
           it != v2.vertices_end(); ++it) {
        compute_vertex(it);
      }
    }

    // Draw the dual vertices
    if (m_draw_dual_vertices) {
      for (Dual_vertices_iterator it = v2.dual().finite_vertices_begin();
           it != v2.dual().finite_vertices_end(); ++it) {
        compute_dual_vertex(it);
      }
    }

    // Add segments and update bounding box
    for (typename V2::Halfedge_iterator it = v2.halfedges_begin();
         it != v2.halfedges_end(); ++it) {
      add_segments_and_update_bounding_box(it);
    }

    for (typename V2::Halfedge_iterator it = v2.halfedges_begin();
         it != v2.halfedges_end(); ++it) {
      compute_rays_and_bisectors(it);
    }

    if (!m_nofaces) {
      for (typename V2::Face_iterator it = v2.faces_begin();
           it != v2.faces_end(); ++it) {
        compute_face(it);
      }
    }
  }

  virtual void keyPressEvent(QKeyEvent *e)
  {
    /// [Keypress]
    const ::Qt::KeyboardModifiers modifiers = e->modifiers();
    if ((e->key() == ::Qt::Key_R) && (modifiers == ::Qt::NoButton)) {
      m_draw_rays = !m_draw_rays;
      displayMessage(
          QString("Draw rays=%1.").arg(m_draw_rays ? "true" : "false"));
      update();
    } else if ((e->key() == ::Qt::Key_V) && (modifiers == ::Qt::NoButton)) {
      m_draw_voronoi_vertices = !m_draw_voronoi_vertices;
      displayMessage(
          QString("Voronoi vertices=%1.").arg(m_draw_voronoi_vertices? "true" : "false"));
      compute_elements();
      redraw();
    } else if ((e->key() == ::Qt::Key_D) && (modifiers == ::Qt::NoButton)) {
      m_draw_dual_vertices = !m_draw_dual_vertices;
      displayMessage(QString("Dual vertices=%1.")
                         .arg(m_draw_dual_vertices ? "true" : "false"));
      compute_elements();
      redraw();
    } else {
      // Call the base method to process others/classicals key
      Base::keyPressEvent(e);
    }
    /// [Keypress]
  }

protected:
  const V2 &v2;
  bool m_nofaces;
  bool m_draw_voronoi_vertices;
  bool m_draw_dual_vertices;
  const ColorFunctor &m_fcolor;
};

// Specialization of draw function.
#define CGAL_VORONOI_TYPE CGAL::Voronoi_diagram_2 <DG, AT, AP>

template<class DG,
         class AT,
         class AP>
void draw(const CGAL_VORONOI_TYPE &av2,
          const char *title="2D Voronoi Diagram Basic Viewer",
          bool nofill = false,
          bool draw_voronoi_vertices = true,
          bool draw_dual_vertices = true)
{
#if defined(CGAL_TEST_SUITE)
  bool cgal_test_suite = true;
#else
  bool cgal_test_suite = qEnvironmentVariableIsSet("CGAL_TEST_SUITE");
#endif

  if (!cgal_test_suite) {
    int argc = 1;
    const char *argv[2] = {"voronoi_2_viewer", "\0"};
    QApplication app(argc, const_cast<char **>(argv));
    DefaultColorFunctorV2 fcolor;
    SimpleVoronoiDiagram2ViewerQt<CGAL_VORONOI_TYPE, DefaultColorFunctorV2>
        mainwindow(app.activeWindow(), av2, title, nofill,
                   draw_voronoi_vertices, draw_dual_vertices, fcolor);
    mainwindow.show();
    app.exec();
  }
}

} // End namespace CGAL

#endif // CGAL_USE_BASIC_VIEWER

#endif // CGAL_DRAW_VORONOI_DIAGRAM_2_H
