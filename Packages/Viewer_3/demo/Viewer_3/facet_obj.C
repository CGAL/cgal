
#include "get_data.h"
#include <CGAL/Viewer_stream.h>

typedef CGAL::Point_3<rep_t> point_t;
typedef CGAL::Tetrahedron_3<rep_t> tetra;



int main(int argc, char *argv[]) 
{
  CGAL::Viewer_3 W(500);
  W.init_window_thread();
  stop();
  Delaunay_3 tr;
  
  tr.insert(point_t(100,100,100));
  tr.insert(point_t(-100,450,-100));
  tr.insert(point_t(0,0,0));
  tr.insert(point_t(200,-130,170));
  tr.insert(point_t(100,-100,-100));
  tr.insert(point_t(400,100,300));
  tr.insert(point_t(210,140,0));
  tr.insert(point_t(-200,50,310));
  tr.insert(point_t(-100,100,-100));
  tr.insert(point_t(100,-300,0));
  tr.insert(point_t(500,-100,-200));
  tr.insert(point_t(200,140,-200));
  tr.insert(point_t(-140,400,-200));
  tr.insert(point_t(100,-350,0));
  tr.insert(point_t(500,-300,-250));
  tr.insert(point_t(300,-100,250));

  CGAL::Drawable_facets_object_3<Delaunay_3>
    df(tr,CGAL::ORANGE,CGAL::BLACK,"triangulation3");
  W.add_drawable(&df);
  W.display();
  pthread_join(W.get_window_thread(), NULL);
  //  W.main_loop();
}
