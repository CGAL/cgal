//example2

#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <iostream.h>
#include <vector>

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h> 


#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/sweep_to_construct_planar_map.h>

#include <CGAL/Arr_polyline_traits.h>

#include <CGAL/IO/Pm_iostream.h>
#include <CGAL/IO/Arr_polyline_traits_iostream.h>


//uncomment if you have LEDA installed.
//#ifndef CGAL_IO_PM_WINDOW_STREAM_H
//#include <CGAL/IO/Pm_Window_stream.h>
//#endif
//#ifndef CGAL_ARR_POLYLINE_TRAITS_WINDOW_STREAM_H  
//#include <CGAL/IO/Arr_polyline_traits_Window_stream.h>
//#endif
//#include <CGAL/IO/leda_window.h>  //used for visualization.

using namespace std;

typedef CGAL::Quotient<int>                  NT;
typedef CGAL::Cartesian<NT>                  R;
typedef CGAL::Arr_polyline_traits<R>         Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

typedef CGAL::Pm_default_dcel<Traits>                Dcel;   
typedef CGAL::Planar_map_2<Dcel, Traits>             PM;

template <class Container>
void read_polylines(Container& curves)
{
  int      num_polylines = 0;
  NT       x,y; 

  std::cin >> num_polylines;
  std::cout << "number of polylines is : " << num_polylines << std::endl;

  while (num_polylines--) {
    Curve      polyline;
    int        num_points;
    
    std::cin >> num_points;

    while (num_points--) {
      std::cin >> x >> y;
      Point s(x, y);
      
      polyline.push_back(s);
    }
    
    curves.push_back(polyline);

    polyline.clear();
  }
}

int main(int argc, char* argv[])
{
  PM                 pm;
  vector<Curve>      polylines;
  
  read_polylines(polylines);

  CGAL::sweep_to_construct_planar_map(polylines.begin(),polylines.end(), pm);
  
  std::cout << pm;
  
  //CGAL::Window_stream W(700, 700);
  //W.init(-10, 10, -10);
  //W.set_mode(leda_src_mode);
  //W.set_node_width(3);
  //W.button("finish",2);
  //W.display();
  //W << pm;
}







