//example3

#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <iostream.h>
//#include <fstream.h>
#include <vector>

#include <CGAL/Quotient.h>
#include <CGAL/Cartesian.h>

#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/sweep_to_produce_planar_map_subcurves.h>

//#include <CGAL/IO/Pm_iostream.h>

typedef CGAL::Quotient<int>              NT;
typedef CGAL::Cartesian<NT>              R;

typedef CGAL::Arr_segment_exact_traits<R>          Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;


using namespace std;


int main(int argc, char* argv[])
{
  int                num_segments;
  list<Curve>        segments;
  
  std::ifstream f_curves(argv[1]);
  cin >> num_segments;
  
  double        x1, y1, x2, y2;
  
  while (num_segments--) {
    cin >> x1 >> y1 >> x2 >> y2;
    
    segments.push_back(Curve(Point(x1,y1), Point(x2,y2)));
  }    
  
  Traits traits;
  list<Curve>  subcurves;
  CGAL::sweep_to_produce_planar_map_subcurves(segments.begin(),segments.end(), traits, subcurves);
  

  for (list<Curve>::iterator scv_iter = subcurves.begin(); scv_iter != subcurves.end(); scv_iter++)
    cout<<*scv_iter<<endl;
}









