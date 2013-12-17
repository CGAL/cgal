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
#include <boost/format.hpp>
#include <CGAL/CGAL_Ipelet_base.h> 
#include <algorithm>
#include <CGAL/point_generators_2.h>
#include <CGAL/copy_n.h>
#include <CGAL/random_selection.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/random_polygon_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Join_input_iterator.h>
#include <CGAL/function_objects.h>
#include <CGAL/copy_n.h>


namespace CGAL_generator{
  
  
const std::string sublabel[] ={
  "Points in a disk","Points on a grid","Points in a square","Points on a convex hull","Polygon","Segments in a square", "Circles (center in a square)","Help"
};

const std::string hlpmsg[] ={
"Generate random inputs. You have to specify the size of the bounding box and the number of elements"};

typedef CGAL::Exact_predicates_inexact_constructions_kernel               Kernel;

struct generator
  : CGAL::Ipelet_base<Kernel,8>
{
  typedef CGAL::Creator_uniform_2<Kernel::FT,Point_2>                     Creator;
  typedef CGAL::Random_points_in_square_2<Point_2,Creator>                Point_generator;
  typedef CGAL::Creator_uniform_2<Kernel::FT,Point_2>                     Pt_creator;
  
  generator()
    : CGAL::Ipelet_base<Kernel,8>("Generators",sublabel, hlpmsg){};
  void protected_run(int);
};


void generator::protected_run(int fn)
{
  if (fn==7) {
    show_help(false);
    return;
  }

  std::list<Point_2> pt_list;
  std::list<Circle_2> cir_list;

  Iso_rectangle_2 bbox=
  read_active_objects(
		      CGAL::dispatch_or_drop_output<Point_2,Circle_2>(
      std::back_inserter(pt_list),
      std::back_inserter(cir_list)
    )
  );

  Kernel::Vector_2 origin;
  double size=200;
  int ret_val;


  if (fn==0){      
    if (cir_list.size()==0) { print_error_message(("Selection must be a circle")); return;}
    Circle_2  circ=*cir_list.begin();
    size =  sqrt(circ.squared_radius());
    origin= circ.center()-CGAL::ORIGIN;
  }else{
    size = (bbox.xmax()-bbox.xmin())/2;
    origin= Kernel::Vector_2((bbox.xmin()+bbox.xmax())/2,(bbox.ymin()+bbox.ymax())/2);
    if (size<1){
      size=200;
      //boost::tie(ret_val,size)=request_value_from_user<int>((boost::format("Size (default : %1%)") % size).str());
      //if (ret_val == -1) return;
      //if (ret_val == 0) size=200;
      origin =  Kernel::Vector_2(200,200);
    }
  }

  int nbelements=30;
  
  boost::tie(ret_val,nbelements)=request_value_from_user<int>((boost::format("Number of elements (default : %1%)") % nbelements).str() );
  if (ret_val == -1) return;
  if (ret_val == 0) nbelements=30;
  

  if(nbelements < 3){
    print_error_message("Not a good value");
    return;
  }
  
  std::vector<Point_2> points;
  std::vector<Segment_2> segments;
  
  if (fn==5)
    points.reserve(nbelements);
  else
    segments.reserve(nbelements);
  
  #ifdef CGAL_USE_IPE_7
  get_IpePage()->deselectAll();
  #else
  get_IpePage()->DeselectAll();
  #endif
  
  switch(fn){
    case 0:{//random point in a circle
      CGAL::Random_points_in_disc_2<Point_2,Creator> gs( size);
      CGAL::cpp11::copy_n( gs, nbelements, std::back_inserter(points));
      }
    break;
    
    case 1://random point on a grid
    points_on_square_grid_2( size, nbelements, std::back_inserter(points),Creator());
    break;
    
    case 6:
    case 2://points in a square : side =   
    {CGAL::Random_points_in_square_2<Point_2, Creator> gc (size);
    CGAL::cpp11::copy_n( gc, nbelements, std::back_inserter(points));
    }
    break;
    
    case 3:{//draw random set of point on a convex hull 
       CGAL::random_convex_set_2(nbelements, std::back_inserter(points),
        Point_generator( size));
    }
    break;
    
    
    case 4:
      // create k-gon and write it into a window: 
      CGAL::random_polygon_2(nbelements, std::back_inserter(points),Point_generator(size));
      for ( std::vector<Point_2>::iterator it=points.begin(); it!=points.end(); ++it) *it = *it + origin; 
      draw_polyline_in_ipe(points.begin(),points.end(),true);
      return;
    
    case 5://Random segments  
    typedef CGAL::Random_points_in_square_2<Point_2, Creator> P1;
    typedef CGAL::Random_points_in_square_2<Point_2, Creator> P2;
    
    P1 p1 (size);
    P2 p2 (size);
    typedef CGAL::Creator_uniform_2< Point_2, Segment_2> Seg_creator;
    typedef CGAL::Join_input_iterator_2< P1, P2, Seg_creator> Seg_iterator;
    Seg_iterator g( p1, p2);
    CGAL::cpp11::copy_n( g, nbelements, std::back_inserter(segments) );
    break;
  };
  
  if (fn==6){
    CGAL::Random random;
    for (std::vector<Point_2>::iterator it_pt=points.begin();it_pt!=points.end();++it_pt)
      draw_in_ipe(Circle_2(*it_pt+origin,pow(random.get_double(size/20.,size/2.),2) ));
    group_selected_objects_();
  }
  else
    if (!points.empty()){// Translate and draw points
      for ( std::vector<Point_2>::iterator it=points.begin(); it!=points.end(); ++it) *it = *it + origin; 
      draw_in_ipe(points.begin(),points.end());
    }
    else
      if (!segments.empty()){// Translate and draw segments
	for ( std::vector<Segment_2>::iterator it=segments.begin(); it!=segments.end(); ++it) 
	  *it = Segment_2( it->source() + origin, it->target() + origin); 
        draw_in_ipe(segments.begin(),segments.end());
      }
}

}

CGAL_IPELET(CGAL_generator::generator)
