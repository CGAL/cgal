// Copyright (c) 2005-2009  INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; version 2.1 of the License.
// See the file LICENSE.LGPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:  $
// $Id:  $
// 
//
// Author(s)     : Sebastien Loriot, Sylvain Pion

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include<boost/shared_ptr.hpp>
#include <CGAL/CGAL_Ipelet_base.h> 
#include<CGAL/create_offset_polygons_2.h>



namespace CGAL_skeleton{

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
  
const std::string Slab[] = {
  "Skeleton", "Polygon offset", "Help"
};

const std::string Hmsg[] = {
  "TODO",
  "TODO",
};

class SkeletonIpelet 
  : public CGAL::Ipelet_base<Kernel,3>{
    typedef CGAL::Polygon_2<Kernel>           Polygon_2_with_vector;
    
  typedef boost::shared_ptr<Polygon_2_with_vector>        PolygonPtr ;
  typedef std::vector<PolygonPtr>             PolygonPtrVector ;    
public:
  SkeletonIpelet()
    :CGAL::Ipelet_base<Kernel,3>("Skeleton and offset",Slab,Hmsg){}
  void protected_run(int);
};


void SkeletonIpelet::protected_run(int fn)
{
  
  if (fn==2) {
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
    Polygon_2_with_vector polygon(itp->vertices_begin(),itp->vertices_end());
      
    if (!polygon.is_simple()){
      print_error_message("Polygon must be simple");
      continue;   
    }
    
    if (polygon.orientation()!=CGAL::COUNTERCLOCKWISE)
      polygon.reverse_orientation();
    
    Kernel::FT lOffset=1;
    
    PolygonPtrVector inner_offset_polygons = CGAL::create_interior_skeleton_and_offset_polygons_2(lOffset,polygon);
    PolygonPtrVector outer_offset_polygons = CGAL::create_exterior_skeleton_and_offset_polygons_2(lOffset,polygon);
    
    

    typedef std::vector< boost::shared_ptr< Polygon_2_with_vector > > PolygonVector ;
    for( PolygonVector::const_iterator pi = inner_offset_polygons.begin() ; pi != inner_offset_polygons.end() ; ++ pi ){
      Polygon_2_with_vector poly=**pi;
      draw_in_ipe(poly);
    }

    
    //~ switch(fn){
    //~ case 0:
    //~ CGAL::y_monotone_partition_2(polygon.vertices_begin(), 
                                  //~ polygon.vertices_end(),
                                  //~ std::back_inserter(partition_polys));
    //~ break;
    
    //~ case 1:
    //~ CGAL::greene_approx_convex_partition_2(polygon.vertices_begin(), 
                                  //~ polygon.vertices_end(),
                                  //~ std::back_inserter(partition_polys));
    //~ break;
    //~ }
    
    //~ std::list<Polygon_2>::const_iterator   poly_it;
    //~ draw_in_ipe(partition_polys.begin(),partition_polys.end());
  }
}

}







CGAL_IPELET(CGAL_skeleton::SkeletonIpelet)

