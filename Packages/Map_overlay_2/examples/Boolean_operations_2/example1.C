#include <CGAL/Cartesian.h> //CGAL definitions that need to be before anything else
#include <CGAL/Quotient.h>
#include <CGAL/Pm_walk_along_line_point_location.h>
#include <CGAL/Bop_default_dcel.h>
#include <CGAL/Map_overlay_default_notifier.h>
#include <CGAL/Map_overlay.h>
#include <CGAL/Boolean_operations_2.h>
#include <CGAL/Arr_segment_exact_traits.h>
#include <CGAL/Planar_map_2.h>
#include <CGAL/sweep_to_construct_planar_map_2.h>
#include <iostream>
#include <vector>
#include <list>

// uncomment if LEDA is installed.
//#include <CGAL/IO/cgal_window.h>  //used for visualization -
//#include <CGAL/IO/Pm_Window_stream.h>

typedef CGAL::Quotient<int>                 NT;
typedef CGAL::Cartesian<NT>                 K;
typedef CGAL::Arr_segment_exact_traits<K>   Traits;
typedef Traits::Point_2                       Point_2;
typedef Traits::X_curve_2                     X_curve_2;
typedef Traits::Curve_2                       Curve_2;

typedef CGAL::Bop_default_dcel<Traits>           Dcel;
typedef CGAL::Planar_map_2<Dcel,Traits>          Planar_map;
typedef CGAL::Map_overlay_2<Planar_map>          MapOverlay;
typedef CGAL::Boolean_operations_2<MapOverlay>   Bops;
typedef CGAL::Pm_walk_along_line_point_location<Planar_map>  PmWalkPL;

typedef Bops::Faces_container                  Faces_container;
typedef Bops::Halfedges_container              Halfedges_container;
typedef Bops::Vertices_container               Vertices_container;


using std::cin;
using std::cout;
using std::endl;

int  main()
{
  Planar_map  pm1, pm2;
  int   num_curves1, num_curves2;
  
  NT   x1, y1, x2, y2;
  
  cin >> num_curves1;
  std::vector<Curve_2>  curves;
  
  while (num_curves1--) {
    cin >> x1 >> y1 >> x2 >> y2;
    
    Point_2 p1(x1,y2), p2(x2,y2);
    curves.push_back(X_curve_2(p1, p2));
  }
  
  Traits traits;
  CGAL::sweep_to_construct_planar_map_2(curves.begin(), 
                                        curves.end(), traits, pm1);
  
  curves.clear();
  
  cin >> num_curves2;
  while (num_curves2--) {
    cin >> x1 >> y1 >> x2 >> y2;
    
    Point_2 p1(x1,y2), p2(x2,y2);
    curves.push_back(X_curve_2(p1, p2));
  }
   
  CGAL::sweep_to_construct_planar_map_2(curves.begin(), 
                                        curves.end(), traits, pm2);
  
  // ignoring unbounded face in boolean operations.
  pm1.unbounded_face()->set_ignore_bop(false); 
  pm2.unbounded_face()->set_ignore_bop(false);
  
  cout<<"Creating a boolean-operation object"<<endl;
  Bops bop(pm1, pm2);
  
  Halfedges_container halfedges_result;
  bop.intersection(halfedges_result);
  
  for (Halfedges_container::iterator h_iter=halfedges_result.begin();
       h_iter!=halfedges_result.end(); ++h_iter)
    cout << (*h_iter)->curve() << endl;

  // uncomment if LEDA is installed.
  //CGAL::Window_stream W(700, 700);
  //W.init(-10,10,-10);
  //W.set_node_width(3);
  //W.display();
  
  //W<<CGAL::BLUE;
  //W<<pm1;
  //W<<CGAL::RED;
  //W<<pm2;

  // drawing all edges in the intersection of pm1 and pm2.
  //for (Halfedges_container::iterator h_iter=halfedges_result.begin();
  //     h_iter!=halfedges_result.end(); ++h_iter){
  //  W << CGAL::PURPLE;
  //  W << (*h_iter)->curve();
  //}
      
  return 0;
}
