// examples/Sweep_line/example1.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>

#include <CGAL/leda_real.h>
#include <CGAL/Arr_segment_exact_traits.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/sweep_to_produce_points_2.h>

#include <CGAL/IO/Pm_iostream.h>

#include <iostream>
#include <vector>

//uncomment if you have LEDA installed.
//#include <CGAL/IO/Segment_circle_Window_stream.h>
//#include <CGAL/IO/Pm_Window_stream.h>
//#include <CGAL/IO/leda_window.h>  //used for visualization.


typedef double                              NT;
typedef CGAL::Cartesian<NT>                 R;
typedef CGAL::Arr_segment_exact_traits<R>   Traits;


typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Pm_default_dcel<Traits>                Dcel;   
typedef CGAL::Planar_map_2<Dcel, Traits>             PM;

using std::cout;
using std::cin;
using std::endl;
using std::vector;


int main()
{
  PM                 pm;
  int                num_segments;
  vector<Curve>      segments;
  
  cout << " * * * Demonstrating a trivial usage of the sweep line algorithm" << endl <<endl;

  cin >> num_segments;
  
  while (num_segments--) {
    
    Curve cv;
    
    cin >> cv;

    segments.push_back(cv);
  }    
  
  vector<Point> points;
  Traits traits;
  CGAL::sweep_to_produce_points_2(segments.begin(),
                                  segments.end(), 
                                  traits, 
                                  std::back_inserter(points));

  cout << " * * * Printing list of all points in the arrangemen induced by the input segments" << endl;
  
  for (vector<Point>::iterator p_iter = points.begin();
       p_iter != points.end(); ++p_iter)
    cout<< *p_iter <<endl;
  
}





