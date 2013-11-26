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
// Author(s)     : Olivier Devillers


#include <CGAL/Exact_circular_kernel_2.h>

#include <CGAL/CGAL_Ipelet_base.h> 
#include <CGAL/Object.h>

#include "include/CGAL_ipelets/pencils.h"


/*
template <typename T>
double prob_2() {
CGAL::Random_points_in_square_2<Point_2> g(1.0);
double prob = 0.0;
for (int i = 0; i < 10000; i++) {
Point_2 p1, p2, p3, p4, p5, p6;
p1 = *g++; p2 = *g++; p3 = *g++;
p4 = *g++; p5 = *g++; p6 = *g++;
// the pi's are points inherited from the Cartesian kernel Point_2, so,
// the orientation predicate can be called on them
if(CGAL::orientation(p1, p2, p3) != CGAL::COUNTERCLOCKWISE) std::swap(p1, p3);
T o1 = T(p1, p2, p3);
if(CGAL::orientation(p4, p5, p6) != CGAL::COUNTERCLOCKWISE) std::swap(p4, p6);
T o2 = T(p4, p5, p6);
typedef typename CGAL::CK2_Intersection_traits<Kernel, T, T>::type
Intersection_result;
std::vector<Intersection_result> res;
CGAL::intersection(o1, o2, std::back_inserter(res));
prob += (res.size() != 0) ? 1.0 : 0.0;
}
return prob/10000.0;
}
int main()
{
std::cout << "What is the probability that two arcs formed by" << std::endl;
std::cout << "three random counterclockwise-oriented points on" << std::endl;
std::cout << "an unit square intersect? (wait a second please)" << std::endl;
std::cout << "The probability is: " << prob_2<Circular_arc_2>() <<
std::endl << std::endl;
std::cout << "And what about the probability that two circles formed by"
<< std::endl;
std::cout << "three random counterclockwise-oriented points on" << std::endl;
std::cout << "an unit square intersect? (wait a second please)" << std::endl;
std::cout << "The probability is: " << prob_2<Circle_2>() << std::endl;
return 0;
}
*/
 
#include <CGAL/Cartesian.h>
namespace CGAL_hyperbolic{


typedef CGAL::Exact_circular_kernel_2 Kernel;

// --------------------------------------------------------------------

const std::string sublabel[] = {
  "Line through two points","Circle by center and point", "Help"
};

const std::string helpmsg[] = {
  "Draw the hyperbolic line trough two points",
  "Draw the hyperbolic bisector of two points",
  "Draw the hyperbolic circle given the center (primary selection) and a point",
};

class hyperbolicIpelet 
  : public CGAL::Ipelet_base<Kernel,4> {
public:
  hyperbolicIpelet() 
    :CGAL::Ipelet_base<Kernel,4>("Hyperbolic",sublabel,helpmsg){}
  void protected_run(int);
};
// --------------------------------------------------------------------

void hyperbolicIpelet::protected_run(int fn)
{
  Circle_2 circ;     //constructed circle:
  Circle_2 p1,p2;
  Circle_2  poincare;
  
  if (fn==3) {
    show_help();
    return;
  } 
  
  std::list<Point_2> pt_list,pt_list1;
  std::list<Circle_2> cir_list,cir_list1;

  int i=get_IpePage()->primarySelection(); 
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

  if (pt_list.empty() || cir_list.empty()){
    print_error_message(("two marks and a circle have to be selected"));
    return;
  }
  
  poincare=*cit;++cit;
  p1=Circle_2(*it,0);
  ++it;
  if (it!=pt_list.end())  {
    p2=Circle_2(*it,0);
    ++it;
  }else{ print_error_message(("two marks and a circle have to be selected")); return;}
  if( (it!=pt_list.end())||(cit!=cir_list.end()))
    { print_error_message(("only two marks and a circle have to be selected")); return;}

  if (fn==2){//primary selection must be a point (p1)
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
    // Circle of pencil generated by p1 p2 orthogonal to poincare
    circ = compute_circle_in_pencil<Kernel>(poincare,p1,p2);
    break; //goto clip
  case 2:
    // Circle of pencil p1 poincare through p2
    circ = compute_circle_in_pencil<Kernel>(p2,poincare,p1);
    draw_in_ipe(circ);
    return;
  }     //end of switch

  // clip circ by poincare 
  std::vector< CGAL::Object > result;
  Kernel::Circular_arc_point_2 S,T;
  std::pair<Kernel::Circular_arc_point_2, unsigned > the_pair;

  CGAL::intersection(circ, poincare, std::back_inserter(result));
  assert (result.size()==2);
  assign(the_pair, result[0]);
  assign(the_pair, result[1]);
  S = the_pair.first;
  T = the_pair.first;
  Point_2 SS(CGAL::to_double(S.x()),CGAL::to_double(S.y()));
  Point_2 TT(CGAL::to_double(T.x()),CGAL::to_double(T.y()));
  Circular_arc_2 arc(circ,SS,TT,circ.orientation());
  draw_in_ipe( arc );

}
}

CGAL_IPELET(CGAL_hyperbolic::hyperbolicIpelet)
