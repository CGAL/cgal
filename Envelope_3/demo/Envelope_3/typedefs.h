// Copyright (c) 2005  Tel-Aviv University (Israel).
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
//
//
// Author(s)     : Baruch Zukerman <baruchzu@post.tau.ac.il>
//                 Efi Fogel       <efif@post.tau.ac.il>

#ifndef CGAL_TYPEDEFS_H
#define CGAL_TYPEDEFS_H

#include <CGAL/Cartesian.h>
#include <CGAL/Lazy_exact_nt.h>
#include <CGAL/envelope_3.h>
#include <CGAL/Env_triangle_traits_3.h>
#include <CGAL/Env_surface_data_traits_3.h>
#include <CGAL/CORE_algebraic_number_traits.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Env_sphere_traits_3.h>
#include <CGAL/Env_plane_traits_3.h>
#include <CGAL/Bbox_2.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_widget_Conic_arc_2.h>
#include <CGAL/IO/Qt_widget_Linear_object_2.h>
#include <CGAL/Arr_overlay_2.h>
#include <CGAL/Sweep_line_2/Arr_overlay_traits_2.h>

#include <qpainter.h>

#ifdef CGAL_USE_GMP
  #include <CGAL/Gmpq.h>
  typedef CGAL::Gmpq                                    Base_nt;
#else
  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>
  typedef CGAL::Quotient<CGAL::MP_Float>                Base_nt;
#endif

typedef CGAL::Lazy_exact_nt<Base_nt>                    Coord_type;
// instead of
//typedef CGAL::Cartesian<Coord_type>                     Kernel;
// workaround for VC++
struct Kernel : public CGAL::Cartesian<Coord_type> {};

typedef Kernel::Segment_2                               Segment_2;

typedef Kernel::Point_2                                 Point_2;
typedef Kernel::Point_3                                 Point_3;
typedef Kernel::Triangle_2                              Triangle_2;
typedef Kernel::Triangle_3                              Triangle_3;

typedef CGAL::Env_triangle_traits_3<Kernel>             Base_traits;
typedef Base_traits::Surface_3                          Base_triangle_3;
typedef CGAL::Env_surface_data_traits_3<Base_traits, CGAL::Color>
                                                        Triangle_traits;
typedef Triangle_traits::Surface_3                      Tri_surface_3;
typedef Triangle_traits::Xy_monotone_surface_3
  Xy_monotone_tri_surface_3;
typedef CGAL::Envelope_diagram_2<Triangle_traits>       Envelope_tri_diagram_2;
typedef Triangle_traits::X_monotone_curve_2             Tri_x_monotone_curve_2;

typedef Envelope_tri_diagram_2::Ccb_halfedge_const_circulator
  Tri_ccb_halfedge_const_circulator;


typedef CGAL::CORE_algebraic_number_traits              Nt_traits;
typedef Nt_traits::Rational                             Rational;
typedef Nt_traits::Algebraic                            Algebraic;
typedef CGAL::Cartesian<Rational>                       Rat_kernel;
typedef Rat_kernel::Point_3                             Rat_point_3;
typedef CGAL::Cartesian<Algebraic>                      Alg_kernel;
typedef Alg_kernel::Point_2                             Sphere_point_2;

// instead of
//typedef CGAL::Arr_conic_traits_2<Rat_kernel, Alg_kernel, Nt_traits>
//                                                        Conic_traits_2;
// workaround for VC++
struct Conic_traits_2 : public CGAL::Arr_conic_traits_2<Rat_kernel,
			                                Alg_kernel,
                                                        Nt_traits> {} ;

typedef CGAL::Env_sphere_traits_3<Conic_traits_2>       Base_sphere_traits_3;
typedef Base_sphere_traits_3::Surface_3                 Base_sphere_3;
typedef Base_sphere_traits_3::Rat_point_3               Rat_point_3;
typedef Base_sphere_traits_3::Rational                  Rational;

typedef CGAL::Env_surface_data_traits_3<Base_sphere_traits_3, CGAL::Color>
                                                        Sphere_traits_3;
typedef Sphere_traits_3::X_monotone_curve_2
  Sphere_x_monotone_curve_2;
typedef Sphere_traits_3::Surface_3                      Sphere_3;
typedef CGAL::Envelope_diagram_2<Sphere_traits_3>
  Envelope_sphere_diagram_2;
typedef Envelope_sphere_diagram_2::Ccb_halfedge_const_circulator
  Sphere_ccb_halfedge_const_circulator;


typedef Kernel::Plane_3                                 Plane_3;
typedef CGAL::Env_plane_traits_3<Kernel>                Base_plane_traits_3;
typedef Base_plane_traits_3::Surface_3                  Base_plane_3;
typedef CGAL::Env_surface_data_traits_3<Base_plane_traits_3, CGAL::Color>
  Plane_traits_3;
typedef Plane_traits_3::Surface_3                       Plane_surface_3;
typedef Plane_traits_3::X_monotone_curve_2              Plane_x_monotone_curve_2;
typedef CGAL::Envelope_diagram_2<Plane_traits_3>        Envelope_plane_diagram_2;
typedef Envelope_plane_diagram_2::Ccb_halfedge_const_circulator
  Plane_ccb_halfedge_const_circulator;
typedef Envelope_plane_diagram_2::Halfedge_const_handle
  Plane_halfedge_const_handle;



typedef CGAL::Cartesian<double>                         Double_kernel;
typedef Double_kernel::Point_2                          Double_point_2;
typedef CGAL::Polygon_2<Double_kernel>                  Polygon_2;

template<class Arrangement, class OutputIterator>
class Faces_visitor :
  public CGAL::Arr_overlay_traits_2<typename Arrangement::Geometry_traits_2,
                                    Arrangement,
                                    Arrangement>
{
private:
  typedef typename Arrangement::Vertex_const_handle     V_const_handle;
  typedef typename Arrangement::Halfedge_const_handle   He_const_handle;
  typedef typename Arrangement::Face_const_handle       F_const_handle;

  typedef typename Arrangement::Vertex_handle           V_handle;
  typedef typename Arrangement::Halfedge_handle         He_handle;
  typedef typename Arrangement::Face_handle             F_handle;

  typedef typename Arrangement::Ccb_halfedge_const_circulator
    Ccb_halfedge_const_circulator;

  F_const_handle m_f1;
  F_const_handle m_f2;
  OutputIterator m_oi;

public:
  Faces_visitor(F_const_handle f1, F_const_handle f2, OutputIterator oi) :
    m_f1(f1),
    m_f2(f2),
    m_oi(oi)
  {}

  void create_vertex(V_const_handle, V_const_handle, V_handle) const {}
  void create_vertex(V_const_handle, He_const_handle, V_handle) const {}
  void create_vertex(V_const_handle, F_const_handle, V_handle) const {}
  void create_vertex(He_const_handle, V_const_handle, V_handle) const {}
  void create_vertex(F_const_handle, V_const_handle, V_handle) const {}
  void create_vertex(He_const_handle, He_const_handle, V_handle) const {}
  void create_edge(He_const_handle, He_const_handle, He_handle) const {}
  void create_edge(He_const_handle, F_const_handle, He_handle) const {}
  void create_edge(F_const_handle, He_const_handle, He_handle) const {}

  void create_face(F_const_handle f1, F_const_handle f2, F_handle f)
  {
    if (f1 == m_f1 && f2 == m_f2) {
      Polygon_2 pgn;
      Ccb_halfedge_const_circulator ccb = f->outer_ccb();
      Ccb_halfedge_const_circulator cc = ccb;
      do {
        pgn.push_back(Double_point_2(CGAL::to_double(cc->source()->
                                                     point().x()),
                                     CGAL::to_double(cc->source()->
                                                     point().y())));
        //created from the outer boundary of the face
      }
      while (++cc !=ccb);
      *m_oi++ = pgn;
    }
  }
};

template <class OutoutIterator>
void construct_polygon(CGAL::Qt_widget* ,
                       Tri_ccb_halfedge_const_circulator ccb,
                       OutoutIterator oi, bool is_unb, bool /* is_hole */)
{
  if (is_unb) return;
  /* running with around the outer of the face and generate from it polygon */

  Polygon_2 pgn;

  Tri_ccb_halfedge_const_circulator cc=ccb;
  do {
    pgn.push_back(Double_point_2(CGAL::to_double(cc->source()->point().x()),
                                 CGAL::to_double(cc->source()->point().y())));
    //created from the outer boundary of the face
  }
  while (++cc !=ccb);

  *oi = pgn;
}

template <class OutoutIterator>
void construct_polygon(CGAL::Qt_widget* w,
                       Sphere_ccb_halfedge_const_circulator ccb,
                       OutoutIterator oi, bool is_unb, bool /* is_hole */)
{
  if (is_unb) return;
  Polygon_2 pgn;

  Sphere_ccb_halfedge_const_circulator cc=ccb;
  do {
    const Sphere_x_monotone_curve_2& cv = cc->curve();
    // Get the co-ordinates of the curve's source and target.
    const double  sx = CGAL::to_double(cv.source().x());
    //const double  sy = CGAL::to_double(cv.source().y());
    const double  tx = CGAL::to_double(cv.target().x());
    //const double  ty = CGAL::to_double(cv.target().y());

    if (cv.orientation() == CGAL::COLLINEAR) {
      // The curve is a segment
      pgn.push_back(Double_point_2(CGAL::to_double(cc->source()->point().x()),
                                   CGAL::to_double(cc->source()->
                                                   point().y())));
    }
    else
    {

      // Draw a curves conic arc: As the arc is x-monotone, its source and its
      // target has the extreme x-coordinates.
      const bool    is_source_left = (sx < tx);
      const int     x_min = is_source_left ? w->x_pixel(sx) : w->x_pixel(tx);
      const int     x_max = is_source_left ? w->x_pixel(tx) : w->x_pixel(sx);
      const int     n = x_max - x_min + 1;

      if (n <= 0) return;

      typedef std::pair<double, double>    App_point_2;
      int                                  i;

      App_point_2  *pts = new App_point_2 [n + 1];
      cv.polyline_approximation (n, pts);

      if (cc->direction() == CGAL::ARR_RIGHT_TO_LEFT) {
        std::reverse(pts, pts+n+1);
      }

      for (i = 0; i < n; i++) {
        pgn.push_back(Double_point_2(pts[i].first, pts[i].second));
      }
      delete[] pts;
    }

    //created from the outer boundary of the face
  } while (++cc != ccb);
  *oi = pgn;
}
template <class OutoutIterator>
void construct_polygon(CGAL::Qt_widget* w,
                       Plane_ccb_halfedge_const_circulator ccb,
                       OutoutIterator oi, bool is_unb, bool is_hole)
{
  if (is_unb && !is_hole)
  {
    Plane_ccb_halfedge_const_circulator curr_ccb = ccb;
    Plane_ccb_halfedge_const_circulator non_fict;
    do {
      if (!curr_ccb->is_fictitious()) {
        non_fict = curr_ccb;
        break;
      }
      ++curr_ccb;
    }
    while(curr_ccb != ccb);

    Polygon_2 pgn;
    if (non_fict == Plane_ccb_halfedge_const_circulator()) {
      // didn't find any non-fictitous edge
      pgn.push_back(Double_point_2(w->x_min(), w->y_min()));
      pgn.push_back(Double_point_2(w->x_min(), w->y_max()));
      pgn.push_back(Double_point_2(w->x_max(), w->y_max()));
      pgn.push_back(Double_point_2(w->x_max(), w->y_min()));
      *oi = pgn;
    }
    else {
      // create two arrangements: one from the ccb of the unbounded face,
      // the second fromthe bounding box of the widget and overlay them
      // construct polygons from the intersection faces.

      Envelope_plane_diagram_2 arr_box;
      Envelope_plane_diagram_2 arr_unb_f;

      Point_2 p1(w->x_min(), w->y_min());
      Point_2 p2(w->x_min(), w->y_max());
      Point_2 p3(w->x_max(), w->y_max());
      Point_2 p4(w->x_max(), w->y_min());

      Plane_x_monotone_curve_2 c1(p1, p2);
      Plane_x_monotone_curve_2 c2(p2, p3);
      Plane_x_monotone_curve_2 c3(p3, p4);
      Plane_x_monotone_curve_2 c4(p4, p1);
      Plane_halfedge_const_handle e=
        CGAL::insert_non_intersecting_curve(arr_box, c1);
      CGAL::insert_non_intersecting_curve(arr_box, c2);
      CGAL::insert_non_intersecting_curve(arr_box, c3);
      CGAL::insert_non_intersecting_curve(arr_box, c4);

      if (e->face()->is_unbounded())
        e = e->twin();

      Plane_halfedge_const_handle he =
        insert_non_intersecting_curve(arr_unb_f, non_fict->curve());
      if (he->direction() != non_fict->direction())
        he = he->twin();

      std::list<Plane_x_monotone_curve_2> cv_list;
      Plane_halfedge_const_handle eh;
      for (eh = non_fict->next(); eh != non_fict; eh = eh->next()) {
        if (!eh->is_fictitious())
          cv_list.push_back(eh->curve());
      }
      CGAL::insert_non_intersecting_curves(arr_unb_f, cv_list.begin(),
                                           cv_list.end());

      Faces_visitor<Envelope_plane_diagram_2, OutoutIterator>
        visitor(e->face(), he->face(), oi);
      Envelope_plane_diagram_2 res;
      CGAL::overlay(arr_box, arr_unb_f, res, visitor);
    }
    return;
  }

  // its a a bounded face
  Polygon_2 pgn;

  Plane_ccb_halfedge_const_circulator cc=ccb;
  do {
    pgn.push_back(Double_point_2(CGAL::to_double(cc->source()->point().x()),
                                 CGAL::to_double(cc->source()->point().y())));
    //created from the outer boundary of the face
  }
  while (++cc !=ccb);

  *oi = pgn;
}


template<class Arr>
void draw_face(CGAL::Qt_widget* w, typename Arr::Face_const_iterator f)
{
  bool is_unb = false;
  if(f->is_unbounded())
    is_unb = true;

  Qt::RasterOp old_rasterop = w->rasterOp();
  w->get_painter().setRasterOp(Qt::XorROP);
  // make polygon from the outer ccb of the face 'f'
  std::list<Polygon_2> pgns;
  construct_polygon(w, f->outer_ccb(), std::back_inserter(pgns), is_unb,
                    false);

  std::list<Polygon_2>::iterator itr;
  for (itr = pgns.begin(); itr != pgns.end(); ++itr) {
    Polygon_2 &pgn = *itr;
    if (!pgn.is_empty()) {
      w->setFilled(true);
      CGAL::Color c = f->surface().data();
      w->setFillColor(c);
      (*w) << pgn ;  // draw the polyong
    }
  }

  typename Arr::Inner_ccb_const_iterator hit;
  for (hit = f->holes_begin(); hit != f->holes_end(); ++hit) {
    pgns.clear();
    construct_polygon(w, *hit, std::back_inserter(pgns), is_unb, true);
    CGAL_assertion(pgns.size() == 1);
    const Polygon_2& hole = pgns.front();
    if (!hole.is_empty())
      (*w) << hole;  // draw the polyong
  }
  w->get_painter().setRasterOp(old_rasterop);
  w->setFilled(false);
}

template <class Arrangement>
void draw_arr(CGAL::Qt_widget* w, const Arrangement& arr, bool draw_v,
              bool draw_e, bool draw_f)
{
  typedef typename Arrangement::Face_const_iterator     Face_const_iterator;
  typedef typename Arrangement::Edge_const_iterator     Edge_const_iterator;
  typedef typename Arrangement::Vertex_const_iterator   Vertex_const_iterator;
  if (draw_f) {
    *w <<  CGAL::BLACK;

    Face_const_iterator fit;
    for (fit = arr.faces_begin(); fit != arr.faces_end(); ++fit) {
      if (! fit->number_of_surfaces()) continue;
      draw_face<Arrangement>(w, fit);
    }
  }

  if (draw_e) {
    *w <<  CGAL::BLUE;
    Edge_const_iterator eit;
    for (eit = arr.edges_begin(); eit != arr.edges_end(); ++eit) {
      *w << eit->curve();
    }
  }

  if (draw_v) {
    *w <<  CGAL::RED;
    Vertex_const_iterator vit;
    for (vit = arr.vertices_begin(); vit != arr.vertices_end(); ++vit) {
      *w << vit->point();
    }
  }
}

#endif
