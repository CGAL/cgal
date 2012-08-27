// Copyright (c) 2005-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 3 of the License,
// or (at your option) any later version.
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
// Author(s)     : Sebastien Loriot, Sylvain Pion
 
//~ #include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Cartesian.h>
#include <CGAL/CGAL_Ipelet_base.h> 
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Object.h>
#include <CGAL/intersections.h>
#include <CGAL/Apollonius_site_2.h>
#include <CGAL/Apollonius_graph_2.h>
#include <CGAL/Apollonius_graph_traits_2.h>


namespace CGAL_diagrams{

//~ typedef CGAL::Exact_predicates_inexact_constructions_kernel                   Kernel;
typedef CGAL::Cartesian<double>                                               Kernel;
typedef CGAL::Segment_Delaunay_graph_traits_2<Kernel>                         Gt;
typedef CGAL::Segment_Delaunay_graph_2<Gt>                                    SDG2;
typedef CGAL::Delaunay_triangulation_2<Kernel>                                Delaunay;
typedef CGAL::Regular_triangulation_euclidean_traits_2<Kernel,Kernel::FT>     RGt;
typedef CGAL::Regular_triangulation_2<RGt>                                    Regular;
//Apollonius
typedef CGAL::Apollonius_graph_traits_2<Kernel>                               AT;
typedef CGAL::Apollonius_graph_2<AT>                                          Apollonius;
typedef Apollonius::Site_2                                                    ASite;  
// --------------------------------------------------------------------

const std::string sublabel[] = {
  "Voronoi","Segment Voronoi skeleton", "Power Diagram", "Apollonius", "Help"
};

const std::string helpmsg[] = {
  "Draw the Voronoi diagram of a set of points and segments, circles and circle arcs",
  "Draw the segment Voronoi diagram except the bisectors between a segment and its own endpoints",
  "Draw the Power diagram of a set of weighted points (circles, points)",
  "Draw the Apollonius diagram of a set of circles"
};

class diagrammeIpelet 
  : public CGAL::Ipelet_base<Kernel,5> {
public:
  diagrammeIpelet() 
    :CGAL::Ipelet_base<Kernel,5>("Diagrams",sublabel,helpmsg){}
  void protected_run(int);
};
// --------------------------------------------------------------------

void diagrammeIpelet::protected_run(int fn)
{
  SDG2 svd;     //Voronoi for segments
  Delaunay dt;     //Voronoi of points
  Regular rt;     //power diagram
  Apollonius apo;     //apollonius
  
  bool b=false;
 
  if (fn==4) {
    show_help();
    return;
  } 
  
  std::list<Point_2> pt_list;
  std::list<Segment_2> sg_list;
  std::list<Circle_2> cir_list;
  
  Iso_rectangle_2 bbox=
  read_active_objects(
    CGAL::dispatch_or_drop_output<Point_2,Polygon_2,Circle_2,Segment_2>(
      std::back_inserter(pt_list),
      segment_grabber(std::back_inserter(sg_list)),
      std::back_inserter(cir_list),
      std::back_inserter(sg_list)
    )
  );
  
  switch(fn){
  case 1:
  case 0:
    //VORONOI
    if (pt_list.empty() && sg_list.empty()){
      print_error_message(("No mark, no segment and no polygon selected"));
      return;
    }
    
    b=!( sg_list.empty() );
    for (std::list<Segment_2>::iterator it=sg_list.begin();it!=sg_list.end();++it)
      svd.insert(it->point(0),it->point(1));

    if (b)
      svd.insert(pt_list.begin(),pt_list.end());
    else
      dt.insert(pt_list.begin(),pt_list.end());
    break;
    
    
  case 2:
    //POWER DIAGRAM
    if (pt_list.empty() && cir_list.empty()){
      print_error_message(("No mark nor circle selected"));
      return;
    }
    
    for (std::list<Circle_2>::iterator it=cir_list.begin();it!=cir_list.end();it++)
      rt.insert(Weighted_point_2(it->center(),it->squared_radius()));
    for (std::list<Point_2>::iterator it=pt_list.begin();it!=pt_list.end();it++)
      rt.insert(Weighted_point_2(*it,0));
  break;
     
     
  case 3:     
    //APOLLONIUS
    if (cir_list.empty()){
      print_error_message(("No circle selected"));
      return;
    }  
    for (std::list<Circle_2>::iterator it=cir_list.begin();it!=cir_list.end();it++)
      apo.insert(ASite(it->center(),sqrt(CGAL::to_double(it->squared_radius()))));
    for (std::list<Point_2>::iterator it=pt_list.begin();it!=pt_list.end();it++)
      apo.insert(ASite(*it,0));
  break;
  }     //end of switch

  Kernel::FT incr_len=(fn<2)?50:75;
  //slightly increase the size of the Bbox
  bbox=Iso_rectangle_2(bbox.min()+Kernel::Vector_2(-incr_len,-incr_len),
                       bbox.max()+Kernel::Vector_2(incr_len,incr_len));
  
  
  
  if(fn<2){     //recover dual objects
    if(b){
      if (fn==0) draw_dual_in_ipe(svd,bbox);
      else draw_skeleton_in_ipe(svd,bbox);
    }
    else draw_dual_in_ipe(dt,bbox);
  }
  if(fn==2) draw_dual_in_ipe(rt,bbox);
  if(fn==3) draw_dual_in_ipe(apo,bbox);
}

}

CGAL_IPELET(CGAL_diagrams::diagrammeIpelet)
