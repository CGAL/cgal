//example1

#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <iostream.h>
#include <vector>


#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/Arr_segment_exact_traits.h>

#include <CGAL/sweep_to_construct_planar_map.h>

#include <CGAL/IO/Pm_iostream.h>

//uncomment if you have LEDA installed.
//#include <CGAL/IO/Pm_Window_stream.h>
//#include <CGAL/IO/leda_window.h>  //used for visualization.

typedef CGAL::Quotient<int>              NT;
typedef CGAL::Cartesian<NT>              R;

typedef CGAL::Arr_segment_exact_traits<R>          Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Pm_default_dcel<Traits>                Dcel;   
typedef CGAL::Planar_map_2<Dcel, Traits>             PM;

using namespace std;


int main(int argc, char* argv[])
{
  PM                 pm;
  int                num_segments;
  vector<Curve>     segments;
  
  cin >> num_segments;
  
  NT        x1, y1, x2, y2;
  
  while (num_segments--) {
    cin >> x1 >> y1 >> x2 >> y2;
    
    segments.push_back(Curve(Point(x1,y1), Point(x2,y2)));
  }    
  
  CGAL::sweep_to_construct_planar_map(segments.begin(),segments.end(), pm);
  
  cout<<pm;
  
  //CGAL::Window_stream W(700, 700);
  //W.init(-10, 10, -10);
  //W.set_mode(leda_src_mode);
  //W.set_node_width(3);
  //W.display();
  //W << pm; 
}
