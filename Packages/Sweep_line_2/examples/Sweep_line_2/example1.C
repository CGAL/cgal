// examples/Sweep_line/example1.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/IO/Pm_iostream.h>
#include <iostream>
#include <vector>


#include <CGAL/IO/Pm_Window_stream.h>
#include <CGAL/IO/cgal_window.h>  //used for visualization.

typedef CGAL::Quotient<CGAL::MP_Float>           NT;
typedef CGAL::Cartesian<NT>                      Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel>   Traits;

typedef Traits::Point                            Point;
typedef Traits::X_curve                          X_curve;
typedef Traits::Curve                            Curve;

typedef CGAL::Pm_default_dcel<Traits>            Dcel;   
typedef CGAL::Planar_map_2<Dcel, Traits>         PM;
typedef CGAL::Pm_file_writer<PM>                 Pm_writer;

int main()
{
  PM                 pm;
  int                num_segments;
  std::vector<Curve> segments;
  
  std::cout << " * * * Demonstrating a trivial use of the sweep line algorithm"
	    << std::endl << std::endl;


  // Read input

  std::cin >> num_segments;
  
  NT x1, y1, x2, y2;
  
  while (num_segments--) 
  {
    std::cin >> x1 >> y1 >> x2 >> y2;
    
    segments.push_back(Curve(Point(x1,y1), Point(x2,y2)));
  }    

  // Construct the planar map  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(segments.begin(),segments.end(),
					traits, pm);

  // Write output 
  
  std::cout << " * * * Printing list of all halfedges of the resulting Planar";
  std::cout << " map" << std::endl;

  Pm_writer verbose_writer(std::cout, pm, true);
  verbose_writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());

  // Use a window visualization
  
  CGAL::Window_stream W(700, 700);
  W.init(-10, 10, -10);
  W.set_mode(CGAL::src_mode);
  W.set_node_width(3);
  W.display();
  W << pm;
}
