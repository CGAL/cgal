// Copyright (c) 2009  INRIA Sophia-Antipolis (France).
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
// Author(s)     : Sebastien Loriot

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>



#include <CGAL/CGAL_Ipelet_base.h> 

// --------------------------------------------------------------------


namespace CGAL_bbox_restriction{

typedef CGAL::Exact_predicates_inexact_constructions_kernel               Kernel;  
  
const std::string sublabel[] ={
  "Bounding box restriction", "Help"
};

const std::string helpmsg[] = {
  "Restrict a set of objects to the bounding box of a set of points."
};

struct hilbertsortIpelet
  : CGAL::Ipelet_base<Kernel,2>
{
  hilbertsortIpelet()
    : CGAL::Ipelet_base<Kernel,2>("Bounding box restriction",sublabel, helpmsg) {}

  void protected_run(int);
};


void hilbertsortIpelet::protected_run(int fn)
{
  if (fn==1) {
    show_help();
    return;
  }


  
  std::vector<Circle_2> cir_list;
  std::vector<Circular_arc_2> arc_list;
  std::vector<Polygon_2> poly_list;
  std::vector<Segment_2> seg_list;
  std::vector<Point_2>  pt_list;

  read_active_objects( 
      CGAL::dispatch_or_drop_output<Point_2,Circle_2,Polygon_2,Circular_arc_2,Segment_2>(
                       std::back_inserter(pt_list),std::back_inserter(cir_list),
                       std::back_inserter(poly_list),std::back_inserter(arc_list),
                       std::back_inserter(seg_list)
      ),true,true
  );

  
  if (pt_list.size()<2) {
    print_error_message("No point selected to define a bounding box");
    return;
  }
  
  CGAL::Bbox_2 bbox_2=pt_list.begin()->bbox();
  for (std::vector<Point_2>::iterator it=pt_list.begin();it!=pt_list.end();++it)
    bbox_2=bbox_2+it->bbox();
  Iso_rectangle_2 bbox(bbox_2.xmin(),bbox_2.ymin(),bbox_2.xmax(),bbox_2.ymax());
  
  
  draw_in_ipe(bbox);
  draw_in_ipe(cir_list.begin(),cir_list.end(),bbox,false);
  draw_in_ipe(seg_list.begin(),seg_list.end(),bbox,false);
  draw_in_ipe(arc_list.begin(),arc_list.end(),bbox,false);
  draw_in_ipe(poly_list.begin(),poly_list.end(),bbox,false);
  
}

}

CGAL_IPELET(CGAL_bbox_restriction::hilbertsortIpelet)
