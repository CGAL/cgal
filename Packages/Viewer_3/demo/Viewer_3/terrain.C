#include <CGAL/basic.h>
#include <CGAL/Gmpz.h>
#include <sys/time.h>
#include <fstream>
#include <stack>
#include <set>
#include <iostream>
#include <vector>
#include <string>


#include <CGAL/Homogeneous.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Point_3.h>


#include <CGAL/Viewer_stream.h>
//#include "demo_win.h"
#ifndef V_UTILS
#include <CGAL/v_utils.h>
#endif 

#include <CGAL/draw_CGAL_Objects.h>

typedef double db;
typedef CGAL::Cartesian<db> Rpp;
typedef CGAL::Triangulation_euclidean_traits_xy_3<Rpp>  Gtt;
typedef Gtt::Point    Point_t ;
typedef Gtt::Segment  Segment_t;
typedef Gtt::Triangle Triangle_t;
typedef CGAL::Triangulation_vertex_base_2<Gtt> Vbt; 
typedef CGAL::Constrained_triangulation_face_base_2<Gtt>   Fbt; 
typedef CGAL::Triangulation_default_data_structure_2<Gtt,Vbt,Fbt> Tdst;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gtt,Tdst> Triangulation;
typedef Triangulation::Face_handle Face_handle;
typedef Triangulation::Vertex_handle Vertex_handle;
typedef Triangulation::Line_face_circulator  Line_face_circulator;
typedef Triangulation::Face_iterator  Face_iterator;

typedef CGAL::Point_3<Rpp> point3;
typedef CGAL::Triangle_3<Rpp> triangle3;


std::list<point3> lpt;


std::list<point3> set_z_zero(std::list<point3> lp)
{

  std::list<point3>::iterator it;
  std::list<point3> res;
  for (it=lp.begin(); it!=lp.end(); it++) 
    res.push_back(point3(it->hx(),it->hy(),-200,1));
  return res;
}
 


int main()
{

  Triangulation dtp;

  std::ifstream is("./terrain_2.dat");
  CGAL::set_ascii_mode(is); 
  std::istream_iterator<Point_t> it(is);
  std::istream_iterator<Point_t> end;
  
  for( ; it != end; it++) {
    dtp.insert( *it);
  }
  
  CGAL::Drawable_triangulation_3<Triangulation> ddtp(dtp,CGAL::PURPLE,CGAL::WHITE);

  CGAL::Viewer_3 W(500);
  W.init_window_thread(); 
  //  W.set_custom_panel(demo_panel);
  W.add_drawable(&ddtp,1);
  W.display();
  stop();
 pthread_join(W.get_window_thread(), NULL);
}
