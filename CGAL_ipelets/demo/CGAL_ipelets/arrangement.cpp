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



#define S_D 2.5f
#include <CGAL/config.h> // to include before NDEBUG is defined, to
                         // workaround the check in the testsuite
#ifndef NDEBUG
#define NDEBUG //points are not on circular arcs
#endif
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Arr_circle_segment_traits_2.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <CGAL/Object.h>
#include <CGAL/CGAL_Ipelet_base.h> 


namespace CGAL_argt{

typedef CGAL::Exact_predicates_exact_constructions_kernel           Kernel;
typedef CGAL::Arr_circle_segment_traits_2<Kernel>                   Traits;
typedef Traits::X_monotone_curve_2                                  X_monotone_curve_2;
typedef Traits::Curve_2                                             Curve_2;
typedef std::list<Curve_2>                                          Curve_list;
typedef std::list<X_monotone_curve_2>                               X_monotone_list;

const std::string sublabel[] = {
  "Segmentation","Help" 
};

const std::string helpmsg[] = {
  "Segmentation of a set of segments, circles and circle arcs"
};


class ArrPolyIpelet 
  : public CGAL::Ipelet_base<Kernel,2>{
public:
  ArrPolyIpelet()
    :CGAL::Ipelet_base<Kernel,2>("Arrangement",sublabel,helpmsg){}
  void protected_run(int);
};

void ArrPolyIpelet::protected_run(int fn){
  if (fn==1) {
    show_help();
    return;
  }

  X_monotone_list output_curves;
  Curve_list input_curves;
  //Argt 
  std::list<Segment_2> sg_list;
  std::list<Circle_2> cir_list;
  std::list<Polygon_2> pol_list;
  std::list<Circular_arc_2> arc_list;
  
  read_active_objects(
    CGAL::dispatch_or_drop_output<Polygon_2,Circle_2,Segment_2,Circular_arc_2>(
      std::back_inserter(pol_list),
      std::back_inserter(cir_list),
      std::back_inserter(sg_list),
      std::back_inserter(arc_list)
    ),
    true,true
  );

  for (std::list<Polygon_2>::iterator it=pol_list.begin();it!=pol_list.end();++it)
    for(Polygon_2::Edge_const_iterator edge_it=it->edges_begin();edge_it!=it->edges_end();++edge_it)
      input_curves.push_back(Curve_2(edge_it->point(0),edge_it->point(1)));
  
  for (std::list<Segment_2>::iterator it=sg_list.begin();it!=sg_list.end();++it)
    input_curves.push_back(Curve_2(it->point(0),it->point(1)));
  
  for (std::list<Circle_2>::iterator it=cir_list.begin();it!=cir_list.end();++it)
    input_curves.push_back(Curve_2(it->center(),sqrt(CGAL::to_double(it->squared_radius()))));
  
  for (std::list<Circular_arc_2>::iterator it=arc_list.begin();it!=arc_list.end();++it)
    input_curves.push_back(
      Curve_2( CGAL::cpp11::get<0>(*it).center(),
               sqrt(CGAL::to_double(CGAL::cpp11::get<0>(*it).squared_radius())),
               CGAL::cpp11::get<3>(*it),
               Traits::Point_2(CGAL::cpp11::get<1>(*it).x(),CGAL::cpp11::get<1>(*it).y()),
               Traits::Point_2(CGAL::cpp11::get<2>(*it).x(),CGAL::cpp11::get<2>(*it).y())
             )
      );

  Traits T;
  CGAL::compute_subcurves(input_curves.begin(),input_curves.end(),std::back_inserter(output_curves),false,T);
  


  
  for (X_monotone_list::iterator it=output_curves.begin();it!=output_curves.end();++it){
    Point_2 S(CGAL::to_double(it->source().x()),CGAL::to_double(it->source().y()));
    Point_2 T(CGAL::to_double(it->target().x()),CGAL::to_double(it->target().y()));
    if (it->is_linear ())
      draw_in_ipe(Segment_2(S,T));
    if (it->is_circular())
      draw_in_ipe(Circular_arc_2(it->supporting_circle(),S,T,it->supporting_circle().orientation()));
  }
  return;
}

}

CGAL_IPELET(CGAL_argt::ArrPolyIpelet)
