// examples/Sweep_line/example3.C
// ------------------------------
#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/sweep_to_produce_subcurves_2.h>
#include <iostream>
#include <vector>

typedef CGAL::Quotient<int>                 NT;
typedef CGAL::Cartesian<NT>                 R;
typedef CGAL::Arr_segment_exact_traits<R>   Traits;
typedef Traits::Point                       Point;
typedef Traits::X_curve                     X_curve;
typedef Traits::Curve                       Curve;


int main()
{
  int                num_segments;
  std::list<Curve>        segments;
  
  std::cin >> num_segments;
  
  NT        x1, y1, x2, y2;
  
  while (num_segments--) {
    std::cin >> x1 >> y1 >> x2 >> y2;
    
    segments.push_back(Curve(Point(x1,y1), Point(x2,y2)));
  }    
  
  Traits traits;
  std::list<Curve>  subcurves;
  CGAL::sweep_to_produce_subcurves_2(segments.begin(), 
                                     segments.end(), 
                                     traits, 
                                     std::back_inserter(subcurves));
  
  std::cout<<"The number of disjoint interior sub segments is "<< subcurves.size();
  std::cout<< std::endl;
  
  for (std::list<Curve>::iterator scv_iter = subcurves.begin(); 
       scv_iter != subcurves.end(); 
       scv_iter++)
    std::cout<< *scv_iter<< std::endl;
}
