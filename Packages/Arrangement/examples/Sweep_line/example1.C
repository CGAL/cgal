// examples/Sweep_line/example1.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>

#include <CGAL/Arr_segment_exact_traits.h>

#include <CGAL/sweep_to_construct_planar_map.h>

#include <CGAL/IO/Pm_iostream.h>

#include <iostream>
#include <vector>

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
typedef CGAL::Pm_file_writer<PM>                     Pm_writer;

using std::cout;
using std::cin;
using std::endl;
using std::vector;

int main()
{
  PM                 pm;
  int                num_segments;
  vector<Curve>      segments;
  
  cout << " * * * Demonstrating a trivial use of the sweep line algorithm" << endl <<endl;

  cin >> num_segments;
  
  NT        x1, y1, x2, y2;
  
  while (num_segments--) {
    cin >> x1 >> y1 >> x2 >> y2;
    
    segments.push_back(Curve(Point(x1,y1), Point(x2,y2)));
  }    
  
  CGAL::sweep_to_construct_planar_map(segments.begin(),segments.end(), pm);
  
  cout << " * * * Printing list of all halfedges of the resulting Planar map" << endl;
  Pm_writer verbose_writer(cout, pm, true);
  verbose_writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());
  
  //CGAL::Window_stream W(700, 700);
  //W.init(-10, 10, -10);
  //W.set_mode(leda_src_mode);
  //W.set_node_width(3);
  //W.display();
  //W << pm; 
}

