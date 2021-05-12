// Copyright (c) 2005-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org)
//
// $URL$
// $Id$
// SPDX-License-Identifier: LGPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Sebastien Loriot, Sylvain Pion

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <CGAL/CGAL_Ipelet_base.h>


namespace CGAL_triangulation{

typedef CGAL::Exact_predicates_inexact_constructions_kernel                 Kernel;
typedef CGAL::Delaunay_triangulation_2<Kernel>                              Delaunay;
typedef CGAL::Regular_triangulation_2<Kernel>                               Regular;

typedef CGAL::Triangulation_vertex_base_2<Kernel>                           Vbplus;
typedef CGAL::Constrained_triangulation_face_base_2<Kernel>                 Fbplus;
typedef CGAL::Triangulation_data_structure_2<Vbplus,Fbplus>                 TDSplus;
typedef CGAL::Exact_intersections_tag                                       Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,TDSplus,Itag>     CDTbase;
typedef CGAL::Constrained_triangulation_plus_2<CDTbase>                     CDTplus;

const std::string Slab[] = {
  "Delaunay", "Constrained Delaunay","Conforming Delaunay","Conforming Gabriel","Regular", "Help"
};

const std::string Hmsg[] = {
  "Draw a Delaunay triangulation of a set of points","Draw a Constrained Delaunay triangulation of a set of points and segments",
  "Draw a conforming Delaunay triangulation of a set of segments and points",
  "Draw a conforming Gabriel triangulation of a set of segments and points",
  "Draw a regular triangulation of a set of weighted points (circles, points)"
};

class triangulationIpelet
  : public CGAL::Ipelet_base<Kernel,6>{
public:
  triangulationIpelet()
    :CGAL::Ipelet_base<Kernel,6>("Triangulations",Slab,Hmsg){}
  void protected_run(int);
};


void triangulationIpelet::protected_run(int fn)
{

  if (fn==5) {
    show_help();
    return;
  }

  std::list<Point_2> pt_list;
  std::list<Segment_2> sg_list;
  std::list<Circle_2> cir_list;
  std::list<Polygon_2> pol_list;

  read_active_objects(
    CGAL::dispatch_or_drop_output<Point_2,Polygon_2,Circle_2,Segment_2>(
      std::back_inserter(pt_list),
      std::back_inserter(pol_list),
      std::back_inserter(cir_list),
      std::back_inserter(sg_list)
    )
  );

  Delaunay dt;
  CDTplus Cdt;
  Regular rt;


  switch(fn){
    case 0://Delaunay
      if (pt_list.empty()) {
        print_error_message("No mark selected");
        return;
      }
      dt.insert(pt_list.begin(),pt_list.end());
      draw_in_ipe(dt);
    break;
    case 2:
    case 3:
    case 1://Constraint delaunay
      if (pt_list.empty() && sg_list.empty()){
        print_error_message("No mark nor polygon selected");
      }
      //insert points
      Cdt.insert(pt_list.begin(),pt_list.end());
      //insert constraints
      for (std::list<Segment_2>::iterator it_seg=sg_list.begin();it_seg!=sg_list.end();++it_seg)
        Cdt.insert_constraint(it_seg->point(0),it_seg->point(1));
      for (std::list<Polygon_2>::iterator it_pol=pol_list.begin();it_pol!=pol_list.end();++it_pol)
        for (Polygon_2::Edge_const_iterator it_edge=it_pol->edges_begin();it_edge!=it_pol->edges_end();++it_edge)
          Cdt.insert_constraint(it_edge->point(0),it_edge->point(1));

      if (!Cdt.is_valid()){
        print_error_message("Invalid constrained triangulation");
        return;
      }

      if (!Cdt.number_of_vertices()) {
        print_error_message("No mark nor polygon selected");
        return;
      }

      if (fn==2) CGAL::make_conforming_Delaunay_2(Cdt);
      if (fn==3) CGAL::make_conforming_Gabriel_2(Cdt);
      draw_in_ipe(Cdt);
    break;
    case 4:      //regular
      if (pt_list.empty() && cir_list.empty()){
        print_error_message("No circle nor mark selected");
      }
      //insert points
      for (std::list<Point_2>::iterator it_pt=pt_list.begin();it_pt!=pt_list.end();++it_pt)
        rt.insert(Weighted_point_2(*it_pt,0));
      //insert circles
      for (std::list<Circle_2>::iterator it_cir=cir_list.begin();it_cir!=cir_list.end();++it_cir)
        rt.insert(Weighted_point_2(it_cir->center(),it_cir->squared_radius()));
      draw_in_ipe(rt);
    break;
  }
}

}

CGAL_IPELET(CGAL_triangulation::triangulationIpelet)

