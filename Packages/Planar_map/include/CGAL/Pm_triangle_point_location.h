// Copyright (c) 2004  Tel-Aviv University (Israel).
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
// $Source$
// $Revision$ $Date$
// $Name$
//
// Author(s)     : Idit Haran <haranidi@post.tau.ac.il>

#ifndef CGAL_PM_TRIANGLE_POINT_LOCATION_H
#define CGAL_PM_TRIANGLE_POINT_LOCATION_H

#include <CGAL/Pm_point_location_base.h>
//#include <CGAL/Planar_map_2/Pm_traits_wrap_2.h>

//#define CGAL_PM_WALK_DEBUG
//#define CGAL_PM_DEBUG



//----------------------------------------------------------
// triangulation includes
//----------------------------------------------------------
// file          : examples/Triangulation_2/constrained.C
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <iostream>
#include <stdio.h>

//----------------------------------------------------------
//Pm includes
//----------------------------------------------------------
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
//#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Pm_default_dcel.h>
//#include <CGAL/Planar_map_2.h>
//#include <CGAL/Pm_with_intersections.h>


////////////////////////////////////////////////////////
//    TRIANGLE STRATEGY
////////////////////////////////////////////////////////

CGAL_BEGIN_NAMESPACE

template <class Planar_map_> class Pm_triangle_point_location : 
  public Pm_point_location_base<Planar_map_> {
public:
  //----------------------------------------------------------
  // Pm Types
  //----------------------------------------------------------

  typedef typename Planar_map::Traits               Traits;
  typedef typename Traits::FT                       FT;
  typedef Cartesian<FT>                             Kernel;
  //typedef Quotient<CGAL::MP_Float>                  Coord_type;
  //typedef Cartesian<Coord_type>                     Kernel;
  //  typedef Arr_segment_cached_traits_2<Kernel>       Traits;

  //----------------------------------------------------------
  // Pm Types
  //----------------------------------------------------------

  typedef typename Kernel::Segment_2                        Segment;
  typedef typename Traits::Point_2                          Point_2;
  typedef typename Traits::Curve_2                          Curve_2;
  typedef typename Traits::X_monotone_curve_2               X_monotone_curve_2;
  typedef Planar_map_                                       Planar_map;
  //typedef Pm_default_dcel<Traits>                         Dcel;
  //typedef typename Planar_map_2<Dcel,Traits>              Planar_map;
  //typedef Planar_map_with_intersections_2<Planar_map>     Pmwx;
  typedef Pm_point_location_base<Planar_map>       Base;
  typedef Pm_triangle_point_location<Planar_map>   Self;
  typedef typename Planar_map::Face_iterator              Face_iterator;
  typedef typename Planar_map::Halfedge_iterator          Halfedge_iterator;
  typedef typename Planar_map::Vertex_iterator              Vertex_iterator;
  typedef typename Planar_map::Edge_iterator                Edge_iterator;
  typedef typename Planar_map::Vertex_handle                Vertex_handle;
  typedef typename Planar_map::Halfedge_handle              Halfedge_handle;
  typedef typename Planar_map::Face_handle                  Face_handle;
  typedef typename Planar_map::Halfedge_around_vertex_circulator 
                                         Halfedge_around_vertex_circulator;
  typedef typename Planar_map::Ccb_halfedge_circulator  Ccb_halfedge_circulator;
  typedef typename Base::Halfedge_handle_iterator       Halfedge_handle_iterator;
  typedef typename Base::Token                          Token;
  typedef typename Planar_map::Locate_type              Locate_type;
  typedef typename Planar_map::Traits_wrap              Traits_wrap;
  //----------------------------------------------------------

  //----------------------------------------------------------
  // Triangulation Types
  //----------------------------------------------------------
  typedef Triangulation_vertex_base_with_info_2<Vertex_handle,Kernel>  Vb;
  typedef Triangulation_face_base_with_info_2<CGAL::Color,Kernel>    Fbt;
  typedef Constrained_triangulation_face_base_2<Kernel,Fbt>          Fb;
  typedef Triangulation_data_structure_2<Vb,Fb>                      TDS;
  typedef Exact_predicates_tag                                       Itag;
  typedef Constrained_Delaunay_triangulation_2<Kernel, TDS, Itag>    CDT;
  typedef typename CDT::Point                                  CDT_Point;
  typedef typename CDT::Edge                                   CDT_Edge;
  typedef typename CDT::Face_handle                      CDT_Face_handle;
  typedef typename CDT::Vertex_handle                    CDT_Vertex_handle;
  typedef typename CDT::Finite_faces_iterator            CDT_Finite_faces_iterator;
  typedef typename CDT::Finite_vertices_iterator         CDT_Finite_vertices_iterator;
  typedef typename CDT::Finite_edges_iterator            CDT_Finite_edges_iterator;
  typedef typename CDT::Locate_type                      CDT_Locate_type;
  //----------------------------------------------------------


protected:
  typedef const Self* cPLp;

public:
  // Constructor
  Pm_triangle_point_location() : 
    Pm_point_location_base<Planar_map>(),
    pm(0),
    traits(0),
    updated_cdt(false)
  {}
  
  void init(Planar_map & pmp, Traits & tr) 
  {
#ifdef CGAL_PM_DEBUG
    std::cout << "init PL" << std::endl;
#endif

    CGAL_precondition_msg(pm == NULL,
    "Point location instance should be uninitialized "
    "(Do not use the same instance for more than one map).");

    pm = &pmp;
    traits = (Traits_wrap*)(&tr);

    triangulate_pm();
  }

  inline void insert(Halfedge_handle hh, const X_monotone_curve_2 & cv) 
  {insert_to_cdt(hh, cv); }

  Halfedge_handle locate(const typename Planar_map::Traits::Point_2 & p, Locate_type & lt) const;

  Halfedge_handle locate(const typename Planar_map::Traits::Point_2 & p, Locate_type & lt);

  Halfedge_handle vertical_ray_shoot(const typename Planar_map::Traits::Point_2& p, Locate_type& lt, bool up)
    const;
  Halfedge_handle vertical_ray_shoot(const typename Planar_map::Traits::Point_2& p, Locate_type& lt, bool up);

  inline void split_edge(const X_monotone_curve_2 &, Halfedge_handle, Halfedge_handle,
                         //additions by iddo for arrangement
                         const X_monotone_curve_2 &, const X_monotone_curve_2 &) 
  {updated_cdt = false; triangulate_pm(); }

  inline void merge_edge(const X_monotone_curve_2 &, const X_monotone_curve_2 &, Halfedge_handle, 
                         //additions by iddo for arrangement
                         const X_monotone_curve_2 &)   
  {updated_cdt = false;  triangulate_pm();}

  inline void remove_edge(Halfedge_handle) 
  {updated_cdt = false;  triangulate_pm();}

  inline void remove_edge(const Halfedge_handle_iterator &,
                          const Halfedge_handle_iterator &) 
  {updated_cdt = false; triangulate_pm(); };

  inline void clear() 
  {updated_cdt = false; triangulate_pm();}

  inline void update(const Halfedge_handle_iterator &,
                     const Halfedge_handle_iterator &,
                     const Token& token) 
  {updated_cdt = false; triangulate_pm(); }

  //function that does the triangulation
  void triangulate_pm() ;
  void insert_to_cdt(Halfedge_handle hh, const typename Planar_map::Traits::X_monotone_curve_2 &cv) ;
  
private:

#ifdef CGAL_PM_DEBUG

  void debug() {}

  void debug(const Halfedge_handle& e) const
    {
      {
        if (e!=pm->halfedges_end()) 
          std::cerr << "(" << e->source()->point() << "," 
		    << e->target()->point() << ")" << std::flush;
        else std::cerr << "(oo)";
      }
    }

#endif

public:
  inline const Traits * get_traits() const {return traits;}

protected:
  Planar_map * pm;
  Traits_wrap * traits;
  CDT cdt;
  bool updated_cdt;
};
  
CGAL_END_NAMESPACE

//#ifdef CGAL_CFG_NO_AUTOMATIC_TEMPLATE_INCLUSION
#include <CGAL/Pm_triangle_point_location.C>
//#endif

#endif //PM_TRIANGLE_POINT_LOCATION_H
