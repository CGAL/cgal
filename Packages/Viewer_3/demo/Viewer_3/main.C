#include <CGAL/Cartesian.h>
#include <CGAL/Viewer_stream.h>

#include "custom_win.h"
#include "myhandler.h"


typedef CGAL::Point_3<rep_t> point_t;
typedef CGAL::Point_2<rep_t> point2;
typedef CGAL::Segment_3<rep_t> segment;
typedef CGAL::Segment_2<rep_t> segment2;
typedef CGAL::Triangle_3<rep_t> triangle;
typedef CGAL::Triangle_2<rep_t> triangle2;
typedef CGAL::Tetrahedron_3<rep_t> tetra;
typedef CGAL::Line_3<rep_t> Line;
typedef CGAL::Line_2<rep_t> Line2;
typedef CGAL::Ray_3<rep_t> Ray;
typedef CGAL::Ray_2<rep_t> Ray2;
typedef CGAL::Circle_2<rep_t> circle2;
std::list<point_t> all_points;



int main(int argc, char *argv[]) 
{

  CGAL::Viewer_3 W(500);
  W.init_window_thread();
  W.set_custom_panel(my_panel);
  W.set_mouse_push_handler(myhandler);
   stop();


  point_t p1(100,50,0);
  point_t p2(200,50,0);
  point_t p3(300,50,0);

  CGAL::Drawable_point_3<point_t> dp1(p1,CGAL::RED,CGAL::FILL,25,50);

  CGAL::Drawable_point_3<point_t> dp2(p2,CGAL::RED,CGAL::FILL,25,50);
  CGAL::Drawable_point_3<point_t> dp3(p3,CGAL::RED,CGAL::FILL,25,50);


  point_t p4(100,200,100);
  point_t p5(200,200,100);
  point_t p6(300,200,100);

  CGAL::Drawable_point_3<point_t> dp4(p4,CGAL::VIOLET,CGAL::FILL,5,50);
  CGAL::Drawable_point_3<point_t> dp5(p5,CGAL::GRAY,CGAL::FILL,10,50);
  CGAL::Drawable_point_3<point_t> dp6(p6,CGAL::DEEPBLUE,CGAL::FILL,15,50);

  point_t p7(100,300,-100);
  point_t p8(200,300,-100);
  point_t p9(300,300,-100);

  CGAL::Drawable_point_3<point_t> dp7(p7,CGAL::ORANGE,CGAL::FILL,10,15);
  CGAL::Drawable_point_3<point_t> dp8(p8,CGAL::GREEN,CGAL::WIRE,10,15);
  CGAL::Drawable_point_3<point_t> dp9(p9,CGAL::GREEN,CGAL::RAW,10,15);

  point_t p10(100,400,500);
  point_t p11(200,400,500);
  point_t p12(300,400,500);

  CGAL::Drawable_point_3<point_t> dp10(p10,CGAL::WHITE,CGAL::RAW,5,8);
  CGAL::Drawable_point_3<point_t> dp11(p11,CGAL::PURPLE,CGAL::FILL,5,15);
  CGAL::Drawable_point_3<point_t> dp12(p12,CGAL::YELLOW,CGAL::WIRE,50,20);

  Line ln(point_t(300,300,1),point_t(300,300,-1));
  Ray ry1(point_t(300,300,0), point_t(400,300,0));
  Ray ry2(point_t(300,300,0), point_t(400,350,0));

  CGAL::Drawable_line_3<Line> dln(ln,CGAL::RED,CGAL::FILL,5,20);
  CGAL::Drawable_ray_3<Ray> dry1(ry1,CGAL::BLUE,CGAL::FILL,5,10);
  CGAL::Drawable_ray_3<Ray> dry2(ry2,CGAL::BLUE,CGAL::RAW,5,10);


  triangle tr(point_t(10,10,-200),point_t(300,10,-200),point_t(200,200,-100));
  
  tetra
    tet(point_t(100,100,-200),point_t(400,100,-300),point_t(400,300,-100),point_t(250,250, 100));


  CGAL::Drawable_triangle_3<triangle> dtr(tr,CGAL::GRAY,CGAL::FILL);
  CGAL::Drawable_tetrahedron_3<tetra> dtet(tet,CGAL::ORANGE,CGAL::FILL);

  segment s(p10,p11);
  CGAL::Drawable_segment_3<segment> ds1(s,CGAL::VIOLET,CGAL::FILL,5,10);


 std::list<point_t> lp;
 lp.push_back(point_t(10,10,0));
 lp.push_back(point_t(20,20,0));
 lp.push_back(point_t(30,30,0));
 lp.push_back(point_t(40,40,0));
 std::list<point_t>::iterator first=lp.begin(), last=lp.end();
 CGAL::Drawable_points_set_3<std::list<point_t>::iterator,point_t>
   dlp(first,last,CGAL::RED,CGAL::FILL,5,30);


 W.add_drawable(&dp1);
 W.add_drawable(&dp2);
 W.add_drawable(&dp3);
 W.display();
 stop(); 
 W.add_drawable(&dp4);
 W.add_drawable(&dp5);
 W.add_drawable(&dp6);

 W.add_drawable(&dp7);
 W.add_drawable(&dp8);
 W.add_drawable(&dp9);

 W.add_drawable(&dp10,2);
 W.add_drawable(&dp11,2);
 W.add_drawable(&dp12,2);

 W.add_drawable(&dln,3);
 W.add_drawable(&dry1,3);
 W.add_drawable(&dry2,3);


 W.add_drawable(&dtr,4);
 W.add_drawable(&dtet,5);

 W.add_drawable(&ds1,6);

 W.add_drawable(&dlp,7);
 W.display();
 pthread_join(W.get_window_thread(), NULL);
 // W.main_loop();
}

