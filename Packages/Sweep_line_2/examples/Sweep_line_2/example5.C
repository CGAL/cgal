// examples/Sweep_line/example5.C
// ------------------------------

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/sweep_to_produce_points_2.h>
#include <CGAL/IO/Pm_iostream.h>
#include <iostream>
#include <vector>

#include <CGAL/IO/Segment_circle_Window_stream.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/IO/cgal_window.h>

typedef CGAL::Quotient<CGAL::MP_Float>          NT;
typedef CGAL::Cartesian<NT>                     Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel>  Traits;

typedef Traits::Point                           Point;
typedef Traits::Curve                           Curve;

typedef CGAL::Pm_default_dcel<Traits>           Dcel;   
typedef CGAL::Planar_map_2<Dcel,Traits>         PM;

int main()
{
  PM                  pm;
  int                 num_segments;
  std::vector<Curve>  segments;
  
  std::cout << " * * * Demonstrating a trivial usage of the sweep line ";
  std::cout << "algorithm" << std::endl << std::endl;
  
  // Read input
  std::cin >> num_segments;
  
  while (num_segments--) 
  {
    Curve cv;
    std::cin >> cv;
    segments.push_back(cv);
  }    

  // Use a sweep to find the points induced in the arrangement  
  std::vector<Point> points;
  Traits traits;
  CGAL::sweep_to_produce_points_2(segments.begin(),
                                  segments.end(), 
                                  traits, 
                                  std::back_inserter(points));

  // Write output
  std::cout << " * * * Printing list of all points in the arrangement ";
  std::cout << "induced by the input segments" << std::endl;
  
  for (std::vector<Point>::iterator p_iter = points.begin();
       p_iter != points.end(); ++p_iter)
    std::cout<< *p_iter << std::endl;
  
  return 0;
}
