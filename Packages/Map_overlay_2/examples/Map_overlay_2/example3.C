#include <CGAL/Cartesian.h> //CGAL definitions that need to be before anything else
#include <CGAL/Quotient.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Map_overlay_default_dcel.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Arr_polyline_traits.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>

#include <iostream>
#include <vector>
#include <list>

typedef CGAL::Quotient<int>                 NT;
typedef CGAL::Cartesian<NT>                 K;
typedef CGAL::Arr_polyline_traits<K>        Traits;
typedef Traits::Point_2                     Point_2;
typedef Traits::X_curve_2                   X_curve_2;
typedef Traits::Curve_2                     Curve_2;

typedef CGAL::Map_overlay_default_dcel<Traits> Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>        Planar_map;
typedef CGAL::Map_overlay_2<Planar_map>        MapOverlay;

typedef CGAL::Pm_walk_along_line_point_location<Planar_map>  PmWalkPL;

using std::cin;
using std::cout;
using std::endl;

// This is provided dute to the fact that the extractor (is >> ) 
// is not defined per each kind of curve. 
class Read_polyline {
public:
  Curve_2 operator()(std::istream& is) {
    Curve_2 cv;
    
    std::size_t  size;
    
    is >> size;

    for (unsigned int i = 0; i < size; ++i){
      int x,y;
      
      is>>x>>y;
      
      cv.push_back(Point_2(x,y));  
    }
    return cv;
  }
};

int  main()
{
  std::vector<Curve_2> polylines;
  Planar_map   pm1, pm2;
  int   num_curves1, num_curves2;
  
  std::cin >> num_curves1;
  while (num_curves1--) {
    Curve_2 cv = Read_polyline()(std::cin);
    polylines.push_back(cv);
  } 
  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(polylines.begin(),polylines.end(), 
                                        traits, pm1);

  polylines.clear();
  
  std::cin >> num_curves2;
  while (num_curves2--) {
    Curve_2 cv = Read_polyline()(std::cin);
    polylines.push_back(cv); 
  }
  
  CGAL::sweep_to_construct_planar_map_2(polylines.begin(),polylines.end(), 
                                        traits, pm2);
  
  MapOverlay map1(pm1), map2(pm2);
  cout<<"Calling map overlay\n";
  MapOverlay map_overlay(map1, map2);
  
  MapOverlay::Change_notification  notifier;
  for (Planar_map::Vertex_const_iterator v_iter=map_overlay.subdivision().vertices_begin();
         v_iter!=map_overlay.subdivision().vertices_end(); ++v_iter){
    if (notifier.get_first_halfedge_above(v_iter) != v_iter->incident_halfedges() && 
        notifier.get_second_halfedge_above(v_iter) != v_iter->incident_halfedges())
      // means both halfedge above are not null.
      std::cout<<v_iter->point()<<" is an intersection point" << std::endl;
    else
      std::cout<<v_iter->point()<<" is an endpoint" << std::endl;
  }
  
  return 0;
}




