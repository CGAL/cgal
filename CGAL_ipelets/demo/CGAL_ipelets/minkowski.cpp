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

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/CGAL_Ipelet_base.h>
#include <CGAL/minkowski_sum_2.h>
#include <CGAL/approximated_offset_2.h>
#include <CGAL/offset_polygon_2.h>



namespace CGAL_minkowski{

typedef CGAL::Exact_predicates_exact_constructions_kernel                       Kernel;
typedef CGAL::Polygon_with_holes_2<Kernel,std::vector<Kernel::Point_2> >        Polygon_with_holes_2;
typedef CGAL::Gps_circle_segment_traits_2<Kernel>                               Gps_traits_2;
typedef Gps_traits_2::Polygon_2                                                 Offset_polygon_2;
typedef Gps_traits_2::Polygon_with_holes_2                                      Offset_polygon_with_holes_2;


const std::string Slab[] = {
  "Minkowski Sum",
  "Polygon Offset",
  "Help"
};

const std::string Hmsg[] = {
  "Compute the Minkowski sum of two simple polygons. Origin is placed at the min point of the bounding box of the selected objects",
  "Compute the offsets of a simple polygon defined by a set of circles"
};

class SubSelectIpelet
  : public CGAL::Ipelet_base<Kernel,3>{
public:
  SubSelectIpelet()
    :CGAL::Ipelet_base<Kernel,3>("Minkowski Sum",Slab,Hmsg){}
  void protected_run(int);
};


void SubSelectIpelet::protected_run(int fn)
{
  if (fn==2) {
    show_help();
    return;
  }

  std::list<Circle_2> cir_list;
  std::list<Polygon_2> pol_list;

  Iso_rectangle_2 bbox=
    read_active_objects(
      CGAL::dispatch_or_drop_output<Polygon_2,Circle_2>(
        std::back_inserter(pol_list),
        std::back_inserter(cir_list)
      )
    );


  if (fn==0 && pol_list.size()!=2){
    print_error_message("You must select exactly two polygons");
    return;
  }


  std::list<double> r_offsets;
  for (std::list<Circle_2>::iterator it=cir_list.begin();it!=cir_list.end();++it)
    r_offsets.push_back(sqrt(CGAL::to_double(it->squared_radius())));

  IpeMatrix tfm (1,0,0,1,-CGAL::to_double((bbox.min)().x()),-CGAL::to_double((bbox.min)().y()));

  for (std::list<Polygon_2>::iterator it=pol_list.begin();it!=pol_list.end();++it)
    if(!it->is_simple()){
      print_error_message("Polygon(s) must be simple");
    }


  if (fn==0){
    Polygon_2 polygon1=*pol_list.begin();
    Polygon_2 polygon2=*++pol_list.begin();
    Polygon_with_holes_2  sum = minkowski_sum_2 (polygon1, polygon2);
    std::list<Point_2> LP;
    for (Polygon_2::iterator it=sum.outer_boundary().vertices_begin();it!= sum.outer_boundary().vertices_end();++it)
      LP.push_back(*it);
    draw_polyline_in_ipe(LP.begin(),LP.end(),true,false,false);

    for (Polygon_with_holes_2::Hole_const_iterator poly_it = sum.holes_begin(); poly_it != sum.holes_end();
          ++poly_it){
      LP.clear();
      for (Polygon_2::iterator it=poly_it->vertices_begin();it!= poly_it->vertices_end();++it)
        LP.push_back(*it);
      draw_polyline_in_ipe(LP.begin(),LP.end(),true,false,false);
    }

    create_polygon_with_holes(true);
    transform_selected_objects_(tfm);
  }
  else{
    if (r_offsets.size()==0)
      r_offsets.push_back(10);
    for (std::list<Polygon_2>::iterator it_pol=pol_list.begin();it_pol!=pol_list.end();++it_pol){
      for(std::list<double>::iterator it=r_offsets.begin();it!=r_offsets.end();++it){
        Offset_polygon_with_holes_2  offset=approximated_offset_2 (*it_pol, *it, 0.0001);
        std::list<Segment_2> LS;
        for( Offset_polygon_2::Curve_iterator itt=offset.outer_boundary().curves_begin();
          itt!=offset.outer_boundary().curves_end();++itt){
          Point_2 S=Point_2(CGAL::to_double(itt->source().x()),CGAL::to_double(itt->source().y()));
          Point_2 T=Point_2(CGAL::to_double(itt->target().x()),CGAL::to_double(itt->target().y()));
          if (itt->is_linear ())
            LS.push_back(Segment_2(S,T));
          if (itt->is_circular())
            draw_in_ipe(Circular_arc_2(itt->supporting_circle(),S,T,itt->supporting_circle().orientation()));
        }
        draw_in_ipe(LS.begin(),LS.end());
      }
    }
  }
}

}


CGAL_IPELET(CGAL_minkowski::SubSelectIpelet)

