#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Arr_leda_segment_exact_traits.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>
#include <CGAL/Bops/Polygons_bops_traits.h>
#include <CGAL/Polygons_bops_2.h>

#include <iostream.h>
#include <vector>
#include <list>

// uncomment if LEDA is installed.
//#include <CGAL/IO/cgal_window.h>  //used for visualization.
//#include <CGAL/IO/Pm_Window_stream.h>

typedef CGAL::Quotient<int>        NT;
typedef CGAL::Cartesian<NT>        K;
typedef K::Point_2                 Point_2;
typedef CGAL::Polygon_2<K>         Polygon;
typedef Polygon::Segment_2         Segment_2;

typedef CGAL::Polygons_bops_traits_2<K>   Bops_traits;

//typedef CGAL::Arr_segment_exact_traits<K>  Traits;

using std::cin; 
using std::cout; 
using std::endl;

void read_polygon(Polygon& polygon)
{  
  unsigned int n;  

  cin >> n;
  
  int    x, y;
  while (n--) {
    cin >> x >> y;
    polygon.push_back(Point_2(x,y));
  }
}

int  main()
{
  Polygon  poly1, poly2;
  
  read_polygon(poly1);
  read_polygon(poly2);
  
  cout<<"Polygons are"<<endl;
  cout<<poly1<<endl;
  cout<<poly2<<endl;
  
  //Bops_traits traits;
  //BopsTraits bops_traits;
  std::list<Polygon>      polygons;
  std::list<Segment_2>    curves;
  std::list<Point_2>      points;
  
  CGAL::Polygons_union_2<Bops_traits> poly_union;
  
  poly_union(poly1,poly2, 
             std::back_inserter(polygons),
             std::back_inserter(curves),
             std::back_inserter(points));
  
  //CGAL::Union(poly1,poly2,traits,bops_traits,
  //            std::back_inserter(polygons),
  //            std::back_inserter(curves),
  //            std::back_inserter(points));
 
  cout<<"The number of resulting polygons: "<<polygons.size()<<endl;
  
  for (std::list<Polygon>::iterator iter = polygons.begin();
       iter != polygons.end(); ++iter)
    cout << *iter;
  
  cout<<endl;
  
 return 0;
}
