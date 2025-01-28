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
#include <CGAL/CGAL_Ipelet_base.h>
#include <CGAL/Regular_triangulation_2.h>

#include "include/CGAL_ipelets/k_delaunay.h"


namespace CGAL_multi_regular{


typedef  CGAL::Exact_predicates_inexact_constructions_kernel           Kernel;
typedef CGAL::Regular_triangulation_2<Kernel>                          Regular;
typedef Regular::Finite_edges_iterator                                 itEdge;

// --------------------------------------------------------------------

const std::string sublabel[] = {
  "Regular", "Regular 2", "Regular 3","Regular n-1", "Regular k", "Power Diagram", "Power Diagram 2", "Power Diagram 3", "Power Diagram n-1", "Power Diagram k", "Help"
};

const std::string hlpmsg[] = {
  "Generate k-th regular triangulation and k-th dual Power diagram. Note : k must be smaller than the number of input circles."
};

class MregularIpelet
  : public CGAL::Ipelet_base<Kernel,11> {
public:
  MregularIpelet()
    : CGAL::Ipelet_base<Kernel,11>("k-order Regular",sublabel,hlpmsg){}
  void protected_run(int);
};

// --------------------------------------------------------------------


// --------------------------------------------------------------------


void MregularIpelet::protected_run(int fn)
{
  Regular rt;
  std::vector<Weighted_point_2> input_wpt;

  if (fn==10) {
    show_help(false);
    return;
  }

  Iso_rectangle_2 bbox=
    read_active_objects(
      CGAL::dispatch_or_drop_output<Point_2,Circle_2>(
        wpoint_grabber(std::back_inserter(input_wpt)),
        wpoint_grabber(std::back_inserter(input_wpt))
      )
    );

  if (!input_wpt.size()) {
    print_error_message("No circle selected");
    return;
  }


  int order = 0;
  if(fn==0 || fn==5) order = 1;
  if(fn==1 || fn==6) order = 2;
  if(fn==2 || fn==7) order = 3;
  if(fn==3 || fn==8) order = input_wpt.size()-1;;
  if(fn==4 || fn==9){

    int ret_val;
    boost::tie(ret_val,order)=request_value_from_user<int>("Enter order");
    if (ret_val < 0){
      print_error_message("Incorrect value");
      return;
    }
    if(order<1 || order>=(int) input_wpt.size()){
      print_error_message("Not a good order");
      return;
    }
  }
  k_delaunay<Kernel>(rt,input_wpt,order);
  if(fn<5)//Draw k-th regular triangulation
    draw_in_ipe(rt);
  else{//Draw kth Power diagram
    double incr_len=75;
    bbox=Iso_rectangle_2((bbox.min)()+Kernel::Vector_2(-incr_len,-incr_len),
                         (bbox.max)()+Kernel::Vector_2(incr_len,incr_len));
    draw_dual_in_ipe(rt,bbox);        //draw Voronoi Diagram
  }
}

}

CGAL_IPELET(CGAL_multi_regular::MregularIpelet)
