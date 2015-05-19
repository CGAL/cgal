// Copyright (c) 2013  INRIA Sophia Antipolis -  Mediterranee,  (France).
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
// Author(s)     : Olivier Devillers


#include <CGAL/Exact_circular_kernel_2.h>

#include <CGAL/CGAL_Ipelet_base.h> 
#include <CGAL/Object.h>


 
#include <CGAL/Cartesian.h>
namespace CGAL_distance_ipelet{


typedef CGAL::Exact_circular_kernel_2 Kernel;

// --------------------------------------------------------------------

const std::string sublabel[] = {
  "2 marks",
  "2 marks (cm)"
  "2 marks (inch)"
  "Help"
};

const std::string helpmsg[] = {
  "Distance between two marks in ipe screen pts",
  "Distance between two marks in centimeters when printed",
  "Distance between two marks in inches when printed",
};

class distanceIpelet 
  : public CGAL::Ipelet_base<Kernel,4> {
public:
  distanceIpelet() 
    :CGAL::Ipelet_base<Kernel,4>("Distance",sublabel,helpmsg){}
  void protected_run(int);
};
// --------------------------------------------------------------------

void distanceIpelet::protected_run(int fn)
{
  if (fn==3) {
    show_help();
    return;
  } 
  
  std::list<Point_2> pt_list;

  int i=get_IpePage()->primarySelection(); 

  if (i<0) {
    print_error_message(("Nothing selected"));
    return;
  }

  Iso_rectangle_2 bbox=
  read_active_objects(
		      CGAL::dispatch_or_drop_output<Point_2>(
      std::back_inserter(pt_list)
    )
  );

  if (pt_list.empty()) {print_error_message(("No mark selected")); return;}
  std::list<Point_2>::iterator it=pt_list.begin();
  Point_2 p1=*it; ++it;
  if (pt_list.end()==it) {
    print_error_message(("Only one mark selected")); return;}
  Point_2 p2=*it; ++it;
  if (pt_list.end()!=it) {
    print_error_message(("More than two marks selected")); return;}

  double length = sqrt( CGAL::to_double(CGAL::squared_distance(p1,p2)) );
  char message[50];
  if (fn==0)
    sprintf(message,"Distance between marks is %f in ipe pts",length);
  else if (fn==1) 
    sprintf(message,"Distance between marks is %f cm",0.0353*length);
  else if (fn==2) 
    sprintf(message,"Distance between marks is %f inches",0.0139*length);
  print_error_message(message);
  return;
}
}

CGAL_IPELET(CGAL_distance_ipelet::distanceIpelet)
