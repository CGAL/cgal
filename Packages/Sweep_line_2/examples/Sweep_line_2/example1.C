// examples/Sweep_line/example1.C
// ------------------------------

#if defined(CGAL_CFG_NO_LONGNAME_PROBLEM) || defined(_MSC_VER)
// Define shorter names to please linker (g++/egcs)
#define Arrangement_2                           _Ar
#define Cartesian                               _Cr
#define Quotient                                _Qt
#define Planar_map_2                            _PM
#define Point_2                                 _Pt
#define allocator                               _All
#define Pm_default_dcel                         _PDD
#define Arr_segment_exact_traits                _AST
#define Segment_2                               _Sg
#define Pm_change_notification                  _PCN
#define Sweep_curves_base_2                     _SCB
#define _Rb_tree                                _RT
#define Sweep_curves_to_planar_map_utils        _SCPMU
#define X_curve_plus_id                         _XCPI
#define sweep_to_construct_planar_map_2         _SCPM
#define vector_iterator                         _VI
#define Point_plus_handle                       _PPH
#define Intersection_point_node                 _IPN
#endif

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

// #include <CGAL/IO/cgal_window.h>  //used for visualization.
// #include <CGAL/IO/Pm_Window_stream.h>

typedef CGAL::Quotient<CGAL::MP_Float>           NT;
typedef CGAL::Cartesian<NT>                      Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel>   Traits;

typedef Traits::Point_2                          Point_2;
typedef Traits::X_curve_2                        X_curve_2;
typedef Traits::Curve_2                          Curve_2;

typedef CGAL::Pm_default_dcel<Traits>            Dcel;   
typedef CGAL::Planar_map_2<Dcel, Traits>         PM;
typedef CGAL::Pm_file_writer<PM>                 Pm_writer;
 
int main()
{
  PM                   pm;
  int                  num_segments;
  std::vector<Curve_2> segments;
  
  std::cout << " * * * Demonstrating a trivial use of the sweep line algorithm"
	    << std::endl << std::endl;

  // Read input
  std::cin >> num_segments;
  
  NT x1, y1, x2, y2;
  
  while (num_segments--) 
  {
    std::cin >> x1 >> y1 >> x2 >> y2;
    segments.push_back(Curve_2(Point_2(x1, y1), Point_2(x2, y2)));
  }    
  // Construct the planar map  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(segments.begin(), segments.end(),
					traits, pm);

  // Write output 
  std::cout << " * * * Printing list of all halfedges of the resulting Planar";
  std::cout << " map" << std::endl;

  Pm_writer verbose_writer(std::cout, pm, true);
  verbose_writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());
  
  // Use a window visualization
  // CGAL::Window_stream W(700, 700, "CGAL Sweep Line Example");
  // W.init(-10, 10, -10);
  // W.set_node_width(3);
  // W.display();
  // W << pm;
  // W.read_mouse();

  return 0;
}
