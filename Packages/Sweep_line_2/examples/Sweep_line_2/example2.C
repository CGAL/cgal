// examples/Sweep_line/example2.C
// ------------------------------
#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h> 
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>
#include <CGAL/Arr_polyline_traits.h>
#include <CGAL/IO/Pm_iostream.h>
#include <iostream>
#include <vector>

// #include <CGAL/IO/cgal_window.h>
// #include <CGAL/IO/Pm_Window_stream.h>
// #include <CGAL/IO/Arr_polyline_traits_Window_stream.h>

typedef CGAL::Quotient<CGAL::MP_Float>       NT;
typedef CGAL::Cartesian<NT>                  Kernel;
typedef CGAL::Arr_polyline_traits<Kernel>    Traits;

typedef Traits::Point_2                      Point_2;
typedef Traits::Curve_2                      Curve_2;

typedef CGAL::Pm_default_dcel<Traits>        Dcel;   
typedef CGAL::Planar_map_2<Dcel, Traits>     PM;
typedef CGAL::Pm_file_writer<PM>             Pm_writer;

CGAL_BEGIN_NAMESPACE

std::ostream & operator<<(std::ostream & os, const ::Curve_2 & cv)
{
  typedef ::Curve_2::const_iterator Points_iterator;
  
  os << cv.size() << std::endl;
  for (Points_iterator points_iter = cv.begin(); 
       points_iter != cv.end(); points_iter++)
    os << " " << *points_iter;

  return os;
}

std::istream & operator>>(std::istream & in, ::Curve_2 & cv)
{
  std::size_t size;
  in >> size;

  for (unsigned int i = 0; i < size; i++){
    ::Point_2 p;
    in >> p;
    cv.push_back(p);  
  }
  
  return in;
}

CGAL_END_NAMESPACE

// Read polylines from the input

template <class Container>
void read_polylines(Container& curves)
{
  int num_polylines = 0;

  std::cin >> num_polylines;
  std::cout << "number of polylines is : " << num_polylines << std::endl;

  while (num_polylines--) {
    Curve_2 polyline;
    
    std::cin >> polyline;
    curves.push_back(polyline);
  }
}

int main()
{
  PM                   pm;
  std::vector<Curve_2> polylines;
  
  // Read input 
  read_polylines(polylines);

  // Construct the planar map  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(polylines.begin(), polylines.end(), 
                                        traits, pm);

  // Write output 
  std::cout << " * * * Printing list of all halfedges of the resulting ";
  std::cout << "Planar map" << std::endl;
  
  Pm_writer verbose_writer(std::cout, pm, true);
  
  verbose_writer.write_halfedges(pm.halfedges_begin(), pm.halfedges_end());

  // Use a window visualization
  // CGAL::Window_stream W(700, 700);
  // W.init(-10, 150, -5);
  // W.set_node_width(3);
  // W.display();
  // W << pm;
  // W.read_mouse();

  return 0;
}
