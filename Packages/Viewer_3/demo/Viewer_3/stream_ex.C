#include <CGAL/Cartesian.h>
#include <CGAL/Viewer_stream.h>

typedef CGAL::Cartesian<double> rep_t;
typedef CGAL::Point_3<rep_t> point_t;

int main(int argc, char *argv[]) 
{

  CGAL::Viewer_3 W(500);
  W.init_window_thread();
  stop();
  point_t p1(100,50,0);
  point_t p2(200,50,0);
  point_t p3(300,50,0);  
  point_t p4(100,200,100);
  point_t p5(200,200,100);
  point_t p6(300,200,100);
  point_t p7(100,300,-100);
  point_t p8(200,300,-100);
  point_t p9(300,300,-100);
  W << CGAL::set_precision(100);
  W << CGAL::set_color_1(CGAL::RED) << p1 << CGAL::set_precision(10) << p2 
      << CGAL::set_precision(100) << p3;
  W.display();
  stop();
  W << CGAL::set_size(50);
  W << CGAL::set_color_1(CGAL::ORANGE) << p4 << p5 << p6;
  W.display();
  stop();
  W << CGAL::set_style(CGAL::WIRE)<< CGAL::set_precision(10);
  W << CGAL::set_color_1(CGAL::BLUE) << p7 << p8 << p9;
  W.display();
  pthread_join(W.get_window_thread(), NULL);
}
