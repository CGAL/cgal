#include <iostream>
#include <list>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Weighted_point.h>


#include <CGAL/grabbers.h>

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Circle_2 Circle_2;
typedef CGAL::Weighted_point<Point_2,Kernel::FT> Weighted_point;
typedef CGAL::Polygon_2<Kernel> Polygon_2;


int main(int, char*[])
{
  std::list<Point_2> pt_list;

  std::list<Segment_2> sg_list;
  sg_list.push_back(Segment_2(Point_2(0,0),Point_2(1,1)));
  sg_list.push_back(Segment_2(Point_2(0,1),Point_2(1,2)));
  
  std::copy(sg_list.begin(), sg_list.end(),CGAL::internal::point_grabber<Kernel>(std::back_inserter(pt_list)));  

  assert (pt_list.size()==4);
  
  std::list<Polygon_2> pol_list;
  pol_list.push_back(Polygon_2(pt_list.begin(),pt_list.end()));
  pol_list.push_back(Polygon_2(pt_list.begin(),pt_list.end()));

  std::copy(pol_list.begin(), pol_list.end(),CGAL::internal::point_grabber<Kernel>(std::back_inserter(pt_list)));  
  
  assert (pt_list.size()==12);
  
  std::copy(pol_list.begin(), pol_list.end(),CGAL::internal::segment_grabber<Kernel>(std::back_inserter(sg_list)));
            
  assert (sg_list.size()==10);

  std::list<Circle_2> l_cir;
  l_cir.push_back( Circle_2(Point_2(0,0),4) );
  
  std::list<Weighted_point> l_wp;
  
  std::copy(l_cir.begin(), l_cir.end(),CGAL::internal::wpoint_grabber<Kernel>(std::back_inserter(l_wp)));

  std::copy(pt_list.begin(), pt_list.end(),CGAL::internal::wpoint_grabber<Kernel>(std::back_inserter(l_wp)));

  assert (l_wp.size()==13);

  return 0;
}
