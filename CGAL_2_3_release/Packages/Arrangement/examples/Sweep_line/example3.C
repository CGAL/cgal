// examples/Sweep_line/example3.C
// ------------------------------
#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/sweep_to_produce_planar_map_subcurves.h>

#include <iostream>
#include <vector>

typedef CGAL::Quotient<int>              NT;
typedef CGAL::Cartesian<NT>              R;

typedef CGAL::Arr_segment_exact_traits<R>          Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;


using namespace std;


int main()
{
  int                num_segments;
  list<Curve>        segments;
  
  cin >> num_segments;
  
  NT        x1, y1, x2, y2;
  
  while (num_segments--) {
    cin >> x1 >> y1 >> x2 >> y2;
    
    segments.push_back(Curve(Point(x1,y1), Point(x2,y2)));
  }    
  
  Traits traits;
  list<Curve>  subcurves;
  CGAL::sweep_to_produce_planar_map_subcurves(segments.begin(), 
					      segments.end(), 
					      traits, 
					      subcurves);
  
  cout<<"The number of disjoint interior sub segments is "<< subcurves.size();
  cout<<endl;

  for (list<Curve>::iterator scv_iter = subcurves.begin(); 
       scv_iter != subcurves.end(); 
       scv_iter++)
    cout<<*scv_iter<<endl;
}
