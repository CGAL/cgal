//example4

#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <iostream.h>
//#include <fstream.h>
#include <vector>

#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h> 

#include <CGAL/Pm_default_dcel.h>
#include <CGAL/Planar_map_2.h>


#include <CGAL/sweep_to_produce_planar_map_subcurves.h>

#include <CGAL/Arr_polyline_traits.h>
#include <CGAL/IO/Arr_polyline_traits_iostream.h>


using namespace std;

typedef CGAL::Quotient<int>                  NT;
typedef CGAL::Cartesian<NT>                  R;
typedef CGAL::Arr_polyline_traits<R>         Traits;

typedef Traits::Point                                 Point;
typedef Traits::X_curve                               X_curve;
typedef Traits::Curve                                 Curve;

using namespace std;

template <class Container>
void read_polylines(Container& curves)
{
  int      num_polylines = 0;
  double   x,y;

  cin >> num_polylines;
  std::cout<<"number of polylines is : " << num_polylines<<std::endl;

  while (num_polylines--) {
    Curve      polyline;
    int        num_points;
    
    cin >> num_points;

    while (num_points--) {
      cin >> x >> y;
      Point s(x, y);
      
      polyline.push_back(s);
    }
    
    curves.push_back(polyline);

    polyline.clear();
  }
}

int main(int argc, char* argv[])
{
  list<Curve>      polylines;
    
  read_polylines(polylines);
  
  Traits traits;
  list<Curve> subcurves;
  CGAL::sweep_to_produce_planar_map_subcurves(polylines.begin(),polylines.end(), traits, subcurves);  
 
  for (list<Curve>::iterator scv_iter = subcurves.begin(); scv_iter != subcurves.end(); scv_iter++)
    cout<<*scv_iter<<endl;
}




