// file: examples/Sweep_line_2/example2.C

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h> 
#include <CGAL/Sweep_line_2.h> 
#include <CGAL/Arr_segment_cached_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <iostream>
#include <vector>
#include <list>

typedef CGAL::Quotient<CGAL::MP_Float>                NT;
typedef CGAL::Cartesian<NT>                           Kernel;
typedef CGAL::Arr_segment_cached_traits_2<Kernel>     Seg_traits;
typedef CGAL::Arr_polyline_traits_2<Seg_traits>       Traits;

typedef Traits::Point_2                               Point_2;
typedef Traits::Curve_2                               Curve_2;
typedef Traits::X_monotone_curve_2                    X_monotone_curve_2;

typedef std::list<Curve_2>                            Curves_list;
typedef Curves_list::iterator                         Curves_iter;
typedef CGAL::Sweep_line_2<Curves_iter, Traits>       Sweep_line;

typedef std::list<X_monotone_curve_2>                 X_curves_list;
typedef X_curves_list::iterator                       X_curves_iter;

// Read one polyline:
Curve_2 read_polyline ()
{
  std::size_t size;
  std::cout << 
    "enter number of points and then the (x,y) values for each point: ";
  std::cin >> size;

  std::list<Point_2>  pts;
  for (unsigned int i = 0; i < size; i++)
  {
    Traits::Point_2 p;
    std::cin >> p;
    pts.push_back(p);  
  }
  std::cout << std::endl;

  return (Curve_2 (pts.begin(), pts.end()));
}

// Read a list of polylines from the input:
template <class Container>
void read_polylines(Container & curves)
{
  int  num_polylines = 0;

  std::cin >> num_polylines;
  std::cout << "number of polylines is : " << num_polylines << std::endl;

  while (num_polylines--) 
    curves.push_back(read_polyline());

  return;
}

int main()
{
  // Read the input polylines.
  Curves_list polylines;
  read_polylines(polylines);
  
  // Use a sweep to create the sub-curves.  
  Traits traits;
  X_curves_list subcurves;
  Sweep_line  sl(&traits);

  sl.get_subcurves (polylines.begin(), polylines.end(), 
		    std::back_inserter(subcurves));

  // Write the output sub-curves.
  X_curves_iter scv_iter;
  for (scv_iter = subcurves.begin(); scv_iter != subcurves.end(); scv_iter++)  
    std::cout << (*scv_iter) << std::endl;

  return (0);
}
