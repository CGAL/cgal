// examples/Sweep_line/example3.C
// ------------------------------

#include <CGAL/Cartesian.h>
#include <CGAL/MP_Float.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/sweep_to_produce_subcurves_2.h>
#include <iostream>
#include <vector>

typedef CGAL::Quotient<CGAL::MP_Float>         NT;
typedef CGAL::Cartesian<NT>                    Kernel;
typedef CGAL::Arr_segment_exact_traits<Kernel> Traits;

typedef Traits::Point_2                        Point_2;
typedef Traits::Curve_2                        Curve_2;


int main()
{
  // Read input
  int                 num_segments;
  std::list<Curve_2>  segments;
  std::cin >> num_segments;
  
  NT  x1, y1, x2, y2;
  
  while (num_segments--) 
  {
    std::cin >> x1 >> y1 >> x2 >> y2;
    segments.push_back(Curve_2(Point_2(x1, y1), Point_2(x2, y2)));
  }    

  // Use a sweep to create the sub curves  
  Traits traits;
  std::list<Curve_2> subcurves;
  CGAL::sweep_to_produce_subcurves_2(segments.begin(), 
                                     segments.end(), 
                                     traits, 
                                     std::back_inserter(subcurves));
  
  // Write output

  std::cout <<"The number of disjoint interior sub segments is "
            << subcurves.size();
  std::cout << std::endl;
  
  for (std::list<Curve_2>::iterator scv_iter = subcurves.begin(); 
       scv_iter != subcurves.end(); scv_iter++)
    std::cout<< *scv_iter<< std::endl;

  return 0;
}
