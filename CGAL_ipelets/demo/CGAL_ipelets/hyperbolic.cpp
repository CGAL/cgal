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

#include "include/CGAL_ipelets/pencils.h"


 
#include <CGAL/Cartesian.h>
namespace CGAL_hyperbolic{


typedef CGAL::Exact_circular_kernel_2 Kernel;

// --------------------------------------------------------------------

const std::string sublabel[] = {
  "Line through two points",
  "Segment through two points",
  "Bisector of two points",
  "Circle by center and point",
  "Circle center", 
  "Help"
};

const std::string helpmsg[] = {
  "Draw the hyperbolic line trough two points in Poincare disk",
  "Draw the hyperbolic segment trough two points in Poincare disk",
  "Draw the hyperbolic bisector of two points in Poincare disk",
  "Draw the hyperbolic circle given the center (primary selection) and a point in Poincare disk",
  "Draw the hyperbolic center given a circle (primary selection) in Poincare disk",
};

class hyperbolicIpelet 
  : public CGAL::Ipelet_base<Kernel,6> {
public:
  hyperbolicIpelet() 
    :CGAL::Ipelet_base<Kernel,6>("Hyperbolic",sublabel,helpmsg){}
  void protected_run(int);
};
// --------------------------------------------------------------------

void hyperbolicIpelet::protected_run(int fn)
{
  Circle_2 circ;     //constructed circle:
  Circle_2 p1,p2;
  Circle_2  poincare,selected;
  
  if (fn==5) {
    show_help();
    return;
  } 
  
  std::list<Point_2> pt_list,pt_list1;
  std::list<Circle_2> cir_list,cir_list1;

  int i=get_IpePage()->primarySelection(); 

  if (i<0) {
    print_error_message(("No mark or circle selected"));
    return;
  }

  read_one_active_object(get_IpePage()->object(i),CGAL::dispatch_or_drop_output<Point_2,Circle_2>(
       std::back_inserter(pt_list1),
       std::back_inserter(cir_list1))); 

  Iso_rectangle_2 bbox=
  read_active_objects(
		      CGAL::dispatch_or_drop_output<Point_2,Circle_2>(
      std::back_inserter(pt_list),
      std::back_inserter(cir_list)
    )
  );

  
  std::list<Point_2>::iterator it1=pt_list1.begin();
  std::list<Circle_2>::iterator cit1=cir_list1.begin();
  std::list<Point_2>::iterator it=pt_list.begin();
  std::list<Circle_2>::iterator cit=cir_list.begin();

  if (fn!=4){
    if (pt_list.empty() || cir_list.empty()){
      print_error_message(("Two marks and a circle have to be selected"));
      return;
    }
  }else{
    if (cir_list.empty()){
      print_error_message(("Two circles have to be selected"));
      return;
    }
  }
  
  poincare=*cit;++cit;
  if(fn==4){
    if( (cit==cir_list.end()) || (cit1==cir_list1.end())){
      print_error_message(("Two circles have to be selected"));
      return;
    }
    if (*cit1==poincare) poincare=*cit;
    selected=*cit1;
  }else{
    p1=Circle_2(*it,0);
    ++it;
    if (it!=pt_list.end())  {
      p2=Circle_2(*it,0);
      ++it;
    }else{ 
      print_error_message(("Two marks and a circle have to be selected")); 
      return;
    }
    if( (it!=pt_list.end())||(cit!=cir_list.end())){
      print_error_message(("Only two marks and a circle have to be selected")); 
      return;
    }
  }

  if (fn==3){//primary selection must be a point (p1)
    if (pt_list1.empty()){
      print_error_message(("Primary selection must be a mark (center)"));
      return;
    }  
    if (*it1 != p1.center()) {
      //swap
      circ = p1;
      p1 = p2;
      p2 = circ;
    }
    if (*it1 != p1.center()) {
      print_error_message(("Primary selection must be a mark (center)"));
      return;
    }  
  }

  switch(fn){
  case 0:
    // Circle orthogonal to p1, p2, and poincare
    circ = compute_circle_orthogonal<Kernel>(p1,p2,poincare);
    break; //goto clip
  case 1:
    // Circle orthogonal to p1, p2, and poincare
    circ = compute_circle_orthogonal<Kernel>(p1,p2,poincare);
    if (orientation(poincare.center(),p1.center(),p2.center())>0)
      draw_in_ipe(Circular_arc_2(circ,p2.center(),p1.center(),circ.orientation()));
    else if (orientation(poincare.center(),p1.center(),p2.center())==0){
      print_error_message(
	"degenerate case, hyperbolic line is on a diameter of Poincare disk");
      draw_in_ipe(Segment_2(p1.center(),p2.center()));
    }else
      draw_in_ipe(Circular_arc_2(circ,p1.center(),p2.center(),circ.orientation()));
    return;
  case 2:
    // Circle of pencil generated by p1 p2 orthogonal to poincare
    circ = compute_circle_in_pencil<Kernel>(poincare,p1,p2);
    break; //goto clip
  case 3:
    // Circle of pencil p1 poincare through p2
    circ = compute_circle_in_pencil<Kernel>(p2,poincare,p1);
    draw_in_ipe(circ);
    return;
  case 4:
    // Zere radius circle of pencil selected poincare inside
    // translate so that Poincare : x^2+y^2=A
    // and selected : x^2+y^2 -2ax -2by +a^2+ b^2=C
    // look for l  so that l.Poincare + selected has zero radius
    double a=CGAL::to_double(selected.center().x())-CGAL::to_double(poincare.center().x());
    double b=CGAL::to_double(selected.center().y())-CGAL::to_double(poincare.center().y());
    double C=CGAL::to_double(selected.squared_radius());
    double A=CGAL::to_double(poincare.squared_radius());
    double B=A+C-a*a-b*b;
    double delta=B*B-4*A*C;
    double l=(-B+sqrt(delta))/2/A;
    l = 1/(1+l);
    Point_2 center=poincare.center()+ (l*(selected.center()-poincare.center()));
    draw_in_ipe(center);
    return;
  }     //end of switch

  // detect degenerate case
  if (circ==Circle_2()){ 
    Kernel::Vector_2 v;
    if (fn==2) v= Kernel::Vector_2
       (p2.center().y()-p1.center().y(),p2.center().x()-p1.center().x());
    else v=p2.center()-p1.center();
    Kernel::FT sqr_length=poincare.squared_radius() / v.squared_length();
    double length = sqrt( CGAL::to_double(sqr_length) );
    v = Kernel::FT(length)*v;
    Point_2 q1=poincare.center()+ v;
    Point_2 q2=poincare.center()- v;
    print_error_message(
	"degenerate case, hyperbolic line is a diameter of Poincare disk");
    Kernel::Segment_2 s(q1,q2);
    draw_in_ipe(s);
    return;
    }

  // clip circ by poincare 
  std::vector< CGAL::Object > result;
  Kernel::Circular_arc_point_2 L,R;
  std::pair<Kernel::Circular_arc_point_2, unsigned > the_pair;

  CGAL::intersection(circ, poincare, std::back_inserter(result));
  assert (result.size()==2);
  assign(the_pair, result[0]);
  L = the_pair.first;
  assign(the_pair, result[1]);
  R = the_pair.first;
  Point_2 LL(CGAL::to_double(L.x()),CGAL::to_double(L.y()));
  Point_2 RR(CGAL::to_double(R.x()),CGAL::to_double(R.y()));
  assert( LL.x() <= RR.x());
  Circular_arc_2 arc;
  if ( orientation(poincare.center(),circ.center(),LL) >0)
    arc = Circular_arc_2(circ,LL,RR,circ.orientation());
  else arc = Circular_arc_2(circ,RR,LL,circ.orientation());
  draw_in_ipe( arc );
}
}

CGAL_IPELET(CGAL_hyperbolic::hyperbolicIpelet)
