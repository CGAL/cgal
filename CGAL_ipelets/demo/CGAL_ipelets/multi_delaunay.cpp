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

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_face_base_2.h>
#include <CGAL/Regular_triangulation_vertex_base_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include "include/CGAL_ipelets/k_delaunay.h"

namespace CGAL_multi_delaunay{

typedef  CGAL::Exact_predicates_inexact_constructions_kernel                        Kernel;
typedef Kernel::FT                                                                  FT;
typedef CGAL::Delaunay_triangulation_2<Kernel>                                      Delaunay;
typedef CGAL::Regular_triangulation_euclidean_traits_2<Kernel,FT>                   Gt;
typedef CGAL::Regular_triangulation_vertex_base_2<Gt>                               Vb;
typedef CGAL::Triangulation_vertex_base_with_info_2<std::vector<Kernel::Point_2>,Gt,Vb>   VbI;
//~ typedef CGAL::Triangulation_vertex_base_with_info_2<std::list<Kernel::Point_2>,Gt,Vb>         VbI;
typedef CGAL::Regular_triangulation_face_base_2<Kernel  >                           Fb;
typedef CGAL::Triangulation_data_structure_2<VbI,Fb>                                Tds;
typedef CGAL::Regular_triangulation_2<Gt,Tds>                                       RegularI;
typedef CGAL::Regular_triangulation_2<Gt>                                           Regular;
typedef Delaunay::Finite_vertices_iterator                                          itVert;
typedef RegularI::Finite_vertices_iterator                                          RitVert;
typedef Delaunay::Vertex_handle                                                     Vertex;
typedef RegularI::Vertex_handle                                                     VertexI;
typedef Delaunay::Vertex_circulator                                                 VCirculator;
typedef RegularI::Vertex_circulator                                                 RVCirculator;


const std::string sublabel[] = {
  "Delaunay", "Delaunay 2", "Delaunay 3","Delaunay n-1", "Delaunay k", "Voronoi", "Voronoi 2", "Voronoi 3", "Voronoi n-1", "Voronoi k", "Help"
};

const std::string hlpmsg[] = {
"Generate k-th Delaunay triangulation and k-th dual Voronoi diagram. Note : k must be smaller than the number of input points."
};

class MdelaunayIpelet 
  : public CGAL::Ipelet_base<Kernel,11> {
public:
  MdelaunayIpelet()
    : CGAL::Ipelet_base<Kernel,11>("k order Delaunay",sublabel,hlpmsg){}
  void protected_run(int);
};



void MdelaunayIpelet::protected_run(int fn)
{
  Delaunay dt;
  RegularI rti;
  Regular rt;
  //~ std::vector<Point_2> pt_list; I use instead pt_list
  
  if (fn==10){
    show_help(false);
    return;
  }
  
  std::vector<Point_2> pt_list;
  
  Iso_rectangle_2 bbox=read_active_objects( CGAL::dispatch_or_drop_output<Point_2>( std::back_inserter(pt_list) ) );  
  
  if (pt_list.empty()){
    print_error_message("No mark selected");
    return;
  }
  
  
  dt.insert(pt_list.begin(),pt_list.end());
  
  switch(fn){
    case 0://Classical Delauney
      draw_in_ipe(dt);
      break;
    case 1:
    case 6:
    case 2:
    case 7:        //Delaunay and Voronoi for 2nd-3rd order
      for (Delaunay::Finite_edges_iterator it=dt.finite_edges_begin();it!=dt.finite_edges_end();++it){
          Point_2 pt0=it->first->vertex(Delaunay::cw(it->second))->point();
          Point_2 pt1=it->first->vertex(Delaunay::ccw(it->second))->point();
          VertexI vertI_cgal = rti.insert(Weighted_point_2(CGAL::midpoint(pt0,pt1),-CGAL::to_double(CGAL::squared_distance(pt0,pt1))/4.));
          pt_list.clear();
          pt_list.push_back(pt0);
          pt_list.push_back(pt1);
          vertI_cgal -> info() = pt_list;
      }
      if(fn==1){//Delauney 2 : just regular triangulation of all midpoints of delaunay segments with weight minus the squared lenght of the edge divided by 4
        draw_in_ipe(rti);
        break;
      }
      if(fn==2 || fn==7){        //Pour l'order 3
        //CAN WE ITERATE OVER DELAUNEY TRIANGLES???
        //WE MAY COUNT SEVERAL TIME SAME TRIANGLE WITH THE FOLLOWING METHOD
        //iterate over adjacent point in the regular triangulation and compute a new wpoint for those having one commun parent from delaunay
        for (RegularI::Finite_edges_iterator it=rti.finite_edges_begin();it!=rti.finite_edges_end();++it){
          Point_2 pt0_ori0=it->first->vertex(Delaunay::cw(it->second))->info().front();
          Point_2 pt0_ori1=it->first->vertex(Delaunay::cw(it->second))->info().back();
          Point_2 pt1_ori0=it->first->vertex(Delaunay::ccw(it->second))->info().front();
          Point_2 pt1_ori1=it->first->vertex(Delaunay::ccw(it->second))->info().back();
          Point_2 pt3 = Point_2();
          if(CGAL::compare_xy(pt0_ori0,pt1_ori0)==CGAL::EQUAL || CGAL::compare_xy(pt0_ori1,pt1_ori0)==CGAL::EQUAL)
            pt3 = pt1_ori1;
          else
            if(CGAL::compare_xy(pt0_ori0,pt1_ori1)==CGAL::EQUAL || CGAL::compare_xy(pt0_ori1,pt1_ori1)==CGAL::EQUAL)
              pt3 = pt1_ori0;

          if(pt3!=Point_2()) //if adjacent wpoints comed from a delaunay triangle
            rt.insert(Weighted_point_2(CGAL::centroid(pt0_ori0,pt0_ori1,pt3),-CGAL::to_double(CGAL::squared_distance(pt0_ori0,pt0_ori1)+
              CGAL::squared_distance(pt0_ori0,pt3)+CGAL::squared_distance(pt3,pt0_ori1))/9.));
        }
        if(fn==2){//Draw 3th Delauney
          draw_in_ipe(rt);
          break;
        }
      }
    case 3://Delaunay and Voronoi of order n-1
    case 8:
      if(fn==3 ||fn==8){
        int order = pt_list.size()-1;
        double pt_x0 =0;//base to compute centroid of n-1 points
        double pt_y0 =0;//base to compute centroid of n-1 points
        double wt=0;//total weight : sum of all distances between two input points
        for(std::vector<Point_2>::iterator it_pt = pt_list.begin();it_pt!=pt_list.end();++it_pt){
          pt_x0 = pt_x0 + CGAL::to_double((*it_pt).x());
          pt_y0 = pt_y0 + CGAL::to_double((*it_pt).y());
          for(std::vector<Point_2>::iterator it_pt2 = it_pt+1;it_pt2!=pt_list.end();++it_pt2){
            wt = wt + CGAL::to_double(CGAL::squared_distance(*it_pt,*it_pt2));
          }
        }
        for(std::vector<Point_2>::iterator it_pt = pt_list.begin();it_pt!=pt_list.end();++it_pt){
          double w = wt;
          //compute centroid of the set of input points / *it_pt
          double pt_x = pt_x0 - CGAL::to_double((*it_pt).x());
          double pt_y = pt_y0 - CGAL::to_double((*it_pt).y());
          //remove from w the sum of distance of input point from *it_pt
          for(std::vector<Point_2>::iterator it_pt2 = pt_list.begin();it_pt2!=pt_list.end();++it_pt2){   //Weighted_point_2 equivalent
            w = w - CGAL::to_double(CGAL::squared_distance(*it_pt,*it_pt2));
          }
          w = - w / (double) (order*order);
          pt_x = pt_x / (double) order;
          pt_y = pt_y / (double) order;
          rt.insert(Weighted_point_2(Point_2(pt_x,pt_y),w));
        }
        if(fn==3){//Draw (n-1)th  Delaunay
          draw_in_ipe(rt);
          break;
        }
      }
    case 4:
    case 9://k-th Delauney and Voronoi
        if(fn==4 ||fn==9){
          int order;
          int ret_val;
          boost::tie(ret_val,order)=request_value_from_user<int>("Enter order");
          if (ret_val < 0){
            print_error_message("Incorrect value");
            return;  
          }
          int nb_pts = pt_list.size();
          
          if(order<1 || order>=nb_pts){
            print_error_message("Not a good order");
            return;
          }
          k_delaunay<Kernel>(rt,pt_list,order);
          if(fn==4){//Draw k-th delaunay
            draw_in_ipe(rt);
            break;
          }
        }
    case 5://Draw Voronoi diagrams
      if(fn==5) draw_dual_in_ipe(dt,bbox);
      if(fn==6) draw_dual_in_ipe(rti,bbox);
      if(fn==7 || fn==8 || fn==9) draw_dual_in_ipe(rt,bbox);
      break;
  }
}

}

CGAL_IPELET(CGAL_multi_delaunay::MdelaunayIpelet)
