#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
//#include <CGAL/Polygon_2.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>
#include <CGAL/IO/Pm_Window_stream.h>

#include <CGAL/leda_rational.h>

#include <CGAL/Bops/Bops_traits.h>
#include <CGAL/Bops_polygon_2.h>

#include <iostream.h>
#include <vector>
#include <list>

#include <CGAL/IO/cgal_window.h>  //used for visualization.


//typedef CGAL::Quotient<int>         NT;
typedef leda_rational               NT;
typedef CGAL::Cartesian<NT>         K;
typedef K::Point_2                  Point_2;
typedef CGAL::Bops_traits_2<K>      BopsTraits;
typedef BopsTraits::Polygon_2       Polygon_2;
typedef Polygon_2::Segment_2        Segment_2;
typedef CGAL::Arr_segment_exact_traits<K>  Traits;


using std::cin; 
using std::cout; 
using std::endl;

void read_polygon(const char* filename, Polygon_2& polygon)
{  
  std::ifstream file(filename);
  unsigned int n;  

  file >> n;
  
  double    x, y;
  while (n--) {
    file >> x >> y;
    polygon.push_back(Point_2(x,y));
  }
  
  file.close();
}

int  main(int argc, char* argv[])
{
  if (argc != 3) {
    std::cout << "usage: Segment_sweep_ovl_from_file filename\n";
    exit(1);
  }

  Polygon_2  poly1, poly2;
  
  read_polygon(argv[1],poly1);
  read_polygon(argv[2],poly2);
  
  Traits traits;
  BopsTraits bops_traits;
  std::list<Polygon_2>    polygons;
  std::list<Segment_2>    curves;
  std::list<Point_2>      points;
  
  CGAL::intersection(poly1,poly2,traits,bops_traits,
                     std::back_inserter(polygons),
                     std::back_inserter(curves),
                     std::back_inserter(points));
 
  cout<<"The number of resulting polygons: "<<polygons.size()<<endl;
  
  for (std::list<Polygon_2>::iterator iter = polygons.begin();
       iter != polygons.end(); ++iter)
    cout << *iter;
  
 return 0;
}
