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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/CGAL_Ipelet_base.h> 
#include <boost/format.hpp>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesher_2.h>
#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>

namespace CGAL_mesh_2{
  typedef CGAL::Exact_predicates_inexact_constructions_kernel       Kernel;
  typedef CGAL::Triangulation_vertex_base_2<Kernel>                 Vb;
  typedef CGAL::Delaunay_mesh_face_base_2<Kernel>                   Fb;
  typedef CGAL::Triangulation_data_structure_2<Vb, Fb>              Tds;
  typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,Tds>    CDT;
  typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>                  Criteria;
  typedef CGAL::Delaunay_mesher_2<CDT, Criteria>                    Mesher;
    
const std::string sublabel[] ={
  "Mesh_2", "Help"
};

const std::string helpmsg[] = {
  "Mesh a polygon using CGAL::Mesh_2; Use circle centers for seeds"
};

struct IpeletMesh2
  : CGAL::Ipelet_base<Kernel,2>
{
  IpeletMesh2()
    : CGAL::Ipelet_base<Kernel,2>("Mesh_2",sublabel, helpmsg) {}

  void protected_run(int);
};


void IpeletMesh2::protected_run(int fn)
{
  if (fn==1) {
    show_help();
    return;
  }
  
  std::list<Point_2> list_of_seeds;
  
  std::list<Point_2> pt_list;
  std::list<Segment_2> sg_list;
  std::list<Circle_2> cir_list;
  std::list<Polygon_2> pol_list;
  
  Iso_rectangle_2 bbox=
    read_active_objects( 
      CGAL::dispatch_or_drop_output<Point_2,Polygon_2,Circle_2,Segment_2>(
        std::back_inserter(pt_list),
        std::back_inserter(pol_list),
        std::back_inserter(cir_list),
        std::back_inserter(sg_list)
      )
    );

  if (pt_list.empty() and sg_list.empty() and pol_list.empty()) {
    print_error_message("No mark selected");
    return;
  }

  CDT cdt;
  
  for (std::list<Point_2>::iterator it=pt_list.begin();it!=pt_list.end();++it)
    cdt.insert(*it);
  for (std::list<Segment_2>::iterator it=sg_list.begin();it!=sg_list.end();++it)
    cdt.insert_constraint(it->point(0),it->point(1));
  for (std::list<Polygon_2>::iterator it=pol_list.begin();it!=pol_list.end();++it)
    for(Polygon_2::Edge_const_iterator edge_it=it->edges_begin();edge_it!=it->edges_end();++edge_it)
      cdt.insert_constraint(edge_it->point(0),edge_it->point(1));
  for (std::list<Circle_2>::iterator it=cir_list.begin();it!=cir_list.end();++it)
    list_of_seeds.push_back(it->center());

  
  double alpha=0;
  
  int x=static_cast<int>( floor(bbox.max().x()-bbox.min().x()) );
  int y=static_cast<int>( floor(bbox.max().y()-bbox.min().y()) );
  
  int ret_val;
  boost::tie(ret_val,alpha)=request_value_from_user<double>((boost::format("Max edge length (BBox %1%x%2%)") % x % y).str() );
  if (ret_val == -1) return;  
  
  if(alpha<0){
    print_error_message("Not a good value");
    return;
  }
  
  if (list_of_seeds.empty()){
    Mesher mesher(cdt);
    mesher.set_criteria(Criteria(0.125, alpha));
    mesher.refine_mesh();
  }
  else
    CGAL::refine_Delaunay_mesh_2(cdt,list_of_seeds.begin(), list_of_seeds.end(),
      Criteria(0.125, alpha));
  
  
  for (CDT::Finite_edges_iterator it=cdt.finite_edges_begin(); it!=cdt.finite_edges_end();++it)
    if (it->first->is_in_domain() || it->first->neighbor(it->second)->is_in_domain())
      draw_in_ipe(cdt.segment(*it));
    
  group_selected_objects_();
}

}

CGAL_IPELET(CGAL_mesh_2::IpeletMesh2)
