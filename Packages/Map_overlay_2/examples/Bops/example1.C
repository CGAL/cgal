#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>
//#include <CGAL/IO/Pm_Window_stream.h>

//#include <CGAL/leda_rational.h>
#include <CGAL/Bops_polygon_2.h>

#include <iostream.h>
#include <vector>
#include <list>

//#include <CGAL/IO/cgal_window.h>  //used for visualization.


typedef CGAL::Quotient<double>       NT;
//typedef leda_rational              NT;
typedef CGAL::Cartesian<NT>         K;
typedef K::Point_2                  Point;
typedef CGAL::Polygon_2<K>          Polygon;
typedef Polygon::Segment_2          Segment_2;
typedef CGAL::Arr_segment_exact_traits<K>  Traits;


using std::cin; 
using std::cout; 
using std::endl;

void read_polygon(const char* filename, Polygon& polygon)
{  
  std::ifstream file(filename);
  unsigned int n;  

  file >> n;
  
  int    x, y;
  while (n--) {
    file >> x >> y;
    polygon.push_back(Point(x,y));
  }
  
  file.close();
}


int  main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cout << "usage: Segment_sweep_ovl_from_file filename\n";
    exit(1);
  }

  Polygon  poly1, poly2;
  
  read_polygon(argv[1],poly1);
  read_polygon(argv[2],poly2);
  
  Traits traits;
  std::list<Polygon>    polygons;
  std::list<Segment_2>  curves;
  std::list<Point>      points;
  
  CGAL::Union(poly1,poly2,traits,
              std::back_inserter(polygons),
              std::back_inserter(curves),
              std::back_inserter(points));
 
  cout<<"The number of resulting polygons: "<<polygons.size()<<endl;
 
 return 0;
}
