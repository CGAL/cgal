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
#include <CGAL/partition_2.h>
#include <CGAL/CGAL_Ipelet_base.h> 


namespace CGAL_convex_part{

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  
const std::string Slab[] = {
  "Y monotone partition", "Greene's approx Convex Partition","Approx Convex Partition","Optimal Convex Partition", "Help"
};

const std::string Hmsg[] = {
  "Y monotonic decomposition of a polygon",
  "Approximation of convex decomposition of a polygon using Greene's algorithm",
  "Approximation of convex decomposition of a polygon using Hertel and Mehlhorn's algorithm",
  "Optimal convex decomposition of a polygon"
};

class ConvexpartitionIpelet 
  : public CGAL::Ipelet_base<Kernel,5>{
public:
  ConvexpartitionIpelet()
    :CGAL::Ipelet_base<Kernel,5>("Polygon Partition",Slab,Hmsg){}
  void protected_run(int);
};


void ConvexpartitionIpelet::protected_run(int fn)
{
  
  if (fn==4) {
    show_help();
    return;
  }

  std::list<Polygon_2> pol_list;
  read_active_objects( CGAL::dispatch_or_drop_output<Polygon_2>( std::back_inserter(pol_list) ) );

  
  
  if (pol_list.size ()==0){
    print_error_message("No polygon selected");
    return;
  }
  
  for (std::list<Polygon_2>::iterator itp=pol_list.begin();itp!=pol_list.end();++itp){
    //~ Polygon_2 polygon=*itp;
    //~ std::list<Polygon_2> partition_polys;
    CGAL::Polygon_2<Kernel,std::list<Kernel::Point_2> > polygon(itp->vertices_begin(),itp->vertices_end());
    std::list<CGAL::Polygon_2<Kernel,std::list<Kernel::Point_2> > > partition_polys;
    
    if (!polygon.is_simple()){
      print_error_message("Polygon must be simple");
      continue;   
    }
    
    if (polygon.orientation()!=CGAL::COUNTERCLOCKWISE)
      polygon.reverse_orientation();
    
    switch(fn){
    case 0:
    CGAL::y_monotone_partition_2(polygon.vertices_begin(), 
                                  polygon.vertices_end(),
                                  std::back_inserter(partition_polys));
    break;
    
    case 1:
    CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(), 
                                  polygon.vertices_end(),
                                  std::back_inserter(partition_polys));
    break;

    case 2:
    CGAL::approx_convex_partition_2(polygon.vertices_begin(), 
                                  polygon.vertices_end(),
                                  std::back_inserter(partition_polys));
    break;
    
    case 3:
    CGAL::optimal_convex_partition_2(polygon.vertices_begin(), 
                                  polygon.vertices_end(),
                                  std::back_inserter(partition_polys));
    break;
    }
    
    draw_in_ipe(partition_polys.begin(),partition_polys.end());
  }
}

}







CGAL_IPELET(CGAL_convex_part::ConvexpartitionIpelet)

