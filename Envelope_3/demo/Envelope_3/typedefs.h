// Copyright (c) 2005  Tel-Aviv University (Israel).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
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
#include <CGAL/Bbox_2.h> 
#include <CGAL/Bbox_3.h> 
#include <CGAL/IO/Color.h> 
#include <CGAL/Polygon_2.h>
#include <CGAL/IO/Qt_widget_Polygon_2.h>
#include <CGAL/IO/Qt_widget_Conic_arc_2.h>

#include <qpainter.h> 

#ifdef CGAL_USE_GMP

  #include <CGAL/Gmpq.h>

  typedef CGAL::Gmpq                                    Base_nt;

#else

  #include <CGAL/MP_Float.h>
  #include <CGAL/Quotient.h>

  typedef CGAL::Quotient<CGAL::MP_Float>                Base_nt;

#endif

typedef CGAL::Lazy_exact_nt<Base_nt>                  Coord_type;
typedef CGAL::Cartesian<Coord_type>		                Kernel;
typedef Kernel::Segment_2                             Segment_2;

typedef Kernel::Point_2				                        Point_2;
typedef Kernel::Point_3				                        Point_3;
typedef Kernel::Triangle_2				                    Triangle_2;
typedef Kernel::Triangle_3				                    Triangle_3;

typedef CGAL::Env_triangle_traits_3<Kernel>           Base_traits;
typedef Base_traits::Surface_3                        Base_triangle_3;
typedef CGAL::Env_surface_data_traits_3<Base_traits, 
                                        CGAL::Color>  Triangle_traits;
typedef Triangle_traits::Surface_3                    Tri_surface_3;
typedef Triangle_traits::Xy_monotone_surface_3        Xy_monotone_tri_surface_3;
typedef CGAL::Envelope_diagram_2<Triangle_traits>     Envelope_tri_diagram_2;
typedef Triangle_traits::X_monotone_curve_2             Tri_x_monotone_curve_2;

typedef Envelope_tri_diagram_2::Ccb_halfedge_const_circulator 
                                                      Tri_ccb_halfedge_const_circulator;


typedef CGAL::CORE_algebraic_number_traits            Nt_traits;
typedef Nt_traits::Rational                           Rational;
typedef Nt_traits::Algebraic                          Algebraic;
typedef CGAL::Cartesian<Rational>                     Rat_kernel;
typedef Rat_kernel::Point_3                           Rat_point_3;
typedef CGAL::Cartesian<Algebraic>                    Alg_kernel;
typedef Alg_kernel::Point_2                           Sphere_point_2;

typedef CGAL::Arr_conic_traits_2<Rat_kernel,
                                 Alg_kernel,
                                 Nt_traits>           Conic_traits_2;

typedef CGAL::Env_sphere_traits_3<Conic_traits_2>     Base_sphere_traits_3;
typedef Base_sphere_traits_3::Surface_3               Base_sphere_3;
typedef Base_sphere_traits_3::Rat_point_3             Rat_point_3;
typedef Base_sphere_traits_3::Rational                Rational;

typedef CGAL::Env_surface_data_traits_3<Base_sphere_traits_3, 
                                        CGAL::Color>  Sphere_traits_3;
typedef Sphere_traits_3::X_monotone_curve_2           Sphere_x_monotone_curve_2;
typedef Sphere_traits_3::Surface_3                    Sphere_3;
typedef CGAL::Envelope_diagram_2<Sphere_traits_3>     Envelope_sphere_diagram_2;
typedef Envelope_sphere_diagram_2::Ccb_halfedge_const_circulator 
                                                      Sphere_ccb_halfedge_const_circulator;

typedef CGAL::Cartesian<double>                       Double_kernel;
typedef Double_kernel::Point_2                        Double_point_2;
typedef CGAL::Polygon_2<Double_kernel>     Polygon_2;
CGAL::Bbox_2 bbox_2d(const Triangle_3& tri)
{
  const Point_3& p0 = tri.vertex(0);
  const Point_3& p1 = tri.vertex(1);
  const Point_3& p2 = tri.vertex(2);
  Triangle_2 proj_tri(Point_2(p0.x(), p0.y()),
                      Point_2(p1.x(), p1.y()),
                      Point_2(p2.x(), p2.y()));
  return proj_tri.bbox();
}

CGAL::Bbox_2 bbox_2d(const Sphere_3& s)
{
  CGAL::Bbox_3 bbox_3d = s.bbox();
  return CGAL::Bbox_2(bbox_3d.xmin(), bbox_3d.ymin(), bbox_3d.xmax(), bbox_3d.ymax());
}

void construct_polygon(CGAL::Qt_widget* w, Tri_ccb_halfedge_const_circulator ccb, Polygon_2& pgn)
{
  /* running with around the outer of the face and generate from it
    * polygon
    */
  
  Tri_ccb_halfedge_const_circulator cc=ccb;
  do 
  {
    pgn.push_back(Double_point_2(CGAL::to_double(cc->source()->point().x()),
                                 CGAL::to_double(cc->source()->point().y())));
    //created from the outer boundary of the face
  } 
  while (++cc !=ccb);
}

void construct_polygon(CGAL::Qt_widget* w, Sphere_ccb_halfedge_const_circulator ccb, Polygon_2& pgn)
{
  Sphere_ccb_halfedge_const_circulator cc=ccb;
  do
  {
    const Sphere_x_monotone_curve_2& cv = cc->curve();
    // Get the co-ordinates of the curve's source and target.
    const double  sx = CGAL::to_double(cv.source().x());
    //const double  sy = CGAL::to_double(cv.source().y());
    const double  tx = CGAL::to_double(cv.target().x());
    //const double  ty = CGAL::to_double(cv.target().y());

    if (cv.orientation() == CGAL::COLLINEAR)
    {
      // The curve is a segment
      pgn.push_back(Double_point_2(CGAL::to_double(cc->source()->point().x()),
                                   CGAL::to_double(cc->source()->point().y())));
    }
    else
    {

      // Draw a curves conic arc: As the arc is x-monotone, its source and its 
      // target has the extreme x-coordinates.
      const bool    is_source_left = (sx < tx);
      const int     x_min = is_source_left ? w->x_pixel(sx) : w->x_pixel(tx);
      const int     x_max = is_source_left ? w->x_pixel(tx) : w->x_pixel(sx);
      const int     n = x_max - x_min + 1;

      if (n <= 0)
        return ;
      
      typedef std::pair<double, double>    App_point_2;
      int                                  i;
      
      App_point_2  *pts = new App_point_2 [n + 1];
      cv.polyline_approximation (n, pts);

      if(cc->direction() == CGAL::LARGER)
      {
        std::reverse(pts, pts+n+1);
      }

      for (i=0; i < n; i++)
      {
        pgn.push_back(Double_point_2(pts[i].first, pts[i].second));
      }
      delete[] pts;
    }
    
  
    //created from the outer boundary of the face
  } while (++cc != ccb);
}


template<class Arr>
void draw_face(CGAL::Qt_widget* w, typename Arr::Face_const_iterator f)
{
  if(f->is_unbounded())
    return;
  
  Qt::RasterOp old_rasterop = w->rasterOp();
  w->get_painter().setRasterOp(Qt::XorROP);
  // make polygon from the outer ccb of the face 'f'
  Polygon_2 pgn;
  construct_polygon(w, f->outer_ccb(), pgn);


  w->setFilled(true);
  CGAL::Color c = f->surface().data();
  w->setFillColor(c);

  (*w) << pgn ;  // draw the polyong

 
  for(typename Arr::Hole_const_iterator hit = f->holes_begin();
      hit != f->holes_end();
      ++hit)
  {
    Polygon_2 hole;
    construct_polygon(w, *hit, hole);
    (*w) << hole ;  // draw the polyong
  }
  w->get_painter().setRasterOp(old_rasterop);
  w->setFilled(false);
  
}


template <class Arrangement>
void draw_arr(CGAL::Qt_widget* w, const Arrangement& arr)
{
  typedef typename Arrangement::Face_const_iterator     Face_const_iterator;
  typedef typename Arrangement::Edge_const_iterator     Edge_const_iterator;
  typedef typename Arrangement::Vertex_const_iterator   Vertex_const_iterator;

   for(Face_const_iterator fit = arr.faces_begin();
        fit != arr.faces_end();
        ++fit)
    {
      if(! fit->number_of_surfaces())
        continue;

      draw_face<Arrangement>(w, fit);
    }

    *w <<  CGAL::BLUE; 
    for(Edge_const_iterator eit = arr.edges_begin();
        eit != arr.edges_end();
        ++eit)
    {
      //draw_curve(w, eit->curve());
      *w << eit->curve();
    }
   
    *w <<  CGAL::RED; 
    for(Vertex_const_iterator vit = arr.vertices_begin();
        vit != arr.vertices_end();
        ++vit)
    {
      *w << vit->point();
    }
}

#endif
