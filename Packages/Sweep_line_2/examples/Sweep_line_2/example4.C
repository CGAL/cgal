// examples/Sweep_line/example4.C
// ------------------------------

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h> 
#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/sweep_to_produce_subcurves_2.h>
#include <CGAL/Arr_polyline_traits.h>
#include <iostream>
#include <vector>

typedef CGAL::Quotient<CGAL::MP_Float>       NT;
typedef CGAL::Cartesian<NT>                  Kernel;

typedef CGAL::Arr_polyline_traits<Kernel>    Traits;

typedef Traits::Point_2                      Point_2;
typedef Traits::Curve_2                      Curve_2;

CGAL_BEGIN_NAMESPACE

std::ostream & operator<<(std::ostream & os, const Curve_2 & cv)
{
  typedef Curve_2::const_iterator Points_iterator;
  
  os << cv.size() << std::endl;
  for (Points_iterator points_iter = cv.begin(); 
       points_iter != cv.end(); points_iter++)
    os << " " << *points_iter;

  return os;
}

std::istream & operator>>(std::istream & in, Curve_2 & cv)
{
  std::size_t size;
  in >> size;

  for (unsigned int i = 0; i < size; i++){
    Traits::Point_2 p;
    in >> p;
    cv.push_back(p);  
  }
  return in;
}

CGAL_END_NAMESPACE

// Read polylines from the input

template <class Container>
void read_polylines(Container & curves)
{
  int  num_polylines = 0;

  std::cin >> num_polylines;
  std::cout << "number of polylines is : " << num_polylines << std::endl;

  while (num_polylines--) {
    Curve_2 polyline;
    std::cin>>polyline;
    curves.push_back(polyline);
    polyline.clear();
  }
}
int main()
{
  // Read input
  std::list<Curve_2> polylines;
  read_polylines(polylines);
  
  // Use a sweep to create the sub curves  
  Traits traits;
  std::list<Curve_2> subcurves;
  CGAL::sweep_to_produce_subcurves_2(polylines.begin(),polylines.end(), 
                                     traits, std::back_inserter(subcurves));

  // Write output
  for (std::list<Curve_2>::iterator scv_iter = subcurves.begin(); 
       scv_iter != subcurves.end(); scv_iter++)    
    std::cout << *scv_iter << std::endl;

  return 0;
}
