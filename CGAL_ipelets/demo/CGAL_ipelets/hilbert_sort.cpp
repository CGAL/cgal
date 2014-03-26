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
#include <CGAL/hilbert_sort.h>
#include <CGAL/CGAL_Ipelet_base.h> 

// --------------------------------------------------------------------


namespace CGAL_hilbert_sort{

typedef CGAL::Exact_predicates_inexact_constructions_kernel               Kernel;  

const std::string sublabel[] ={
  "Hilbert sorting curve", "Help"
};

const std::string helpmsg[] = {
  "Sort the points along a Hilbert curve"
};

struct hilbertsortIpelet
  : CGAL::Ipelet_base<Kernel,2>
{
  hilbertsortIpelet()
    : CGAL::Ipelet_base<Kernel,2>("Hilbert sort",sublabel, helpmsg) {}

  void protected_run(int);
};


void hilbertsortIpelet::protected_run(int fn)
{
  if (fn==2) {
    show_help();
    return;
  }


  
  std::vector<Point_2> pt_list;

  read_active_objects( CGAL::dispatch_or_drop_output<Point_2>( std::back_inserter(pt_list) ) );

  if (pt_list.empty()) {
    print_error_message("No mark selected");
    return;
  }
  
  if (fn==0) CGAL::hilbert_sort(pt_list.begin(), pt_list.end());
  else       CGAL::hilbert_sort(pt_list.begin(), pt_list.end(), CGAL::Hilbert_sort_middle_policy());

  draw_polyline_in_ipe(pt_list.begin(), pt_list.end());
}

}

CGAL_IPELET(CGAL_hilbert_sort::hilbertsortIpelet)
