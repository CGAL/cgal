#include <CGAL/basic.h> //CGAL definitions that need to be before anything else
#include <CGAL/Cartesian.h>
#include <CGAL/Quotient.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Bop_default_dcel.h>
#include <CGAL/Map_overlay_default_notifier.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Boolean_operations_2.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/IO/Pm_Window_stream.h>
#include <iostream.h>
#include <vector>
#include <list>
#include <CGAL/IO/cgal_window.h>  //used for visualization -


typedef CGAL::Quotient<int>                 NT;
typedef CGAL::Cartesian<NT>                 R;
typedef CGAL::Arr_segment_exact_traits<R>   Traits;
typedef Traits::Point                       Point;
typedef Traits::X_curve                     X_curve;
typedef Traits::Curve                       Curve;

typedef CGAL::Bop_default_dcel<Traits>         Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>        Planar_map;
typedef CGAL::Map_overlay<Planar_map>          MapOverlay;
typedef CGAL::Boolean_operations<MapOverlay>   Bops;
typedef CGAL::Pm_walk_along_line_point_location<Planar_map>  PmWalkPL;

typedef Bops::Faces_container                  Faces_container;
typedef Bops::Halfedges_container              Halfedges_container;
typedef Bops::Vertices_container               Vertices_container;


using namespace std;


int  main()
{
  Planar_map  pm1, pm2;
  int   num_curves1, num_curves2;
  
  NT   x1, y1, x2, y2;
  
  std::cin >> num_curves1;
  std::vector<Curve>  curves;
    
  while (num_curves1--) {
    std::cin >> x1 >> y1 >> x2 >> y2;
    curves.push_back(X_curve(Point(x1,y1), Point(x2,y2)));
  }
  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(curves.begin(), curves.end(), traits, pm1);
  
  curves.clear();
  
  std::cin >> num_curves2;
  while (num_curves2--) {
    std::cin >> x1 >> y1 >> x2 >> y2;
    curves.push_back(X_curve(Point(x1,y1), Point(x2,y2)));
  }
   
  CGAL::sweep_to_construct_planar_map_2(curves.begin(), curves.end(), traits, pm2);
  
  // ignoring unbounded face in boolean operations.
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);
  
  std::cout<<"Creating a boolean-operation object"<<std::endl;
  Bops bop(pm1, pm2);
  
  Halfedges_container halfedges_result;
  bop.intersection(halfedges_result);
  
  CGAL::Window_stream W(700, 700);
  W.init(-10,10,-10);
  W.set_node_width(3);
  W.display();
  
  W<<CGAL::BLUE;
  W<<pm1;
  W<<CGAL::RED;
  W<<pm2;

  // drawing all edges in the intersection of pm1 and pm2.
  for (Halfedges_container::iterator h_iter=halfedges_result.begin();
       h_iter!=halfedges_result.end(); ++h_iter){
    W << CGAL::PURPLE;
    W << (*h_iter)->curve();
  }
      
  return 0;
}


