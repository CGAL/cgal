#include <CGAL/Cartesian.h>

#include <sys/time.h>
#include <iostream>

#include <CGAL/IO/Color.h>
#include <CGAL/Bbox_3.h>

#include "PS_Stream_3.C"

typedef CGAL::Cartesian<double> D;

typedef CGAL::Bbox_3 PS_BBox3;

typedef D::Direction_3      Direction;
typedef D::Point_3          Point3;
typedef D::Point_2          Point2;
typedef D::Plane_3          Plane3;
typedef D::Line_3           Line3;
typedef D::Segment_3        Segment3;
typedef D::Triangle_3       Triangle3;
typedef D::Tetrahedron_3    Tetrahedron;

int main()
{
  PS_BBox3 bb3(-500,-500,-500,500,500,500);
  double x,y,z,lx,ly,lz;
  int filling;
  
  std::cerr << "Filling 0,1,2,3 : ";std::cin >> filling;
  std::string filename;
  std::cerr << "Enter file name: ";std::cin >> filename;
  std::cerr << "Enter view x: ";std::cin >> x;
  std::cerr << "Enter view y: ";std::cin >> y;
  std::cerr << "Enter view z: ";std::cin >> z;
  
  Direction dir(x,y,z);
  std::cerr << "Enter light x: ";std::cin >> lx;
  std::cerr << "Enter light y: ";std::cin >> ly;
  std::cerr << "Enter light z: ";std::cin >> lz;
  
  Direction light(lx,ly,lz);
  CGAL::PS_Stream_3 ps(bb3,dir,light,300,filename.c_str(),CGAL::PS_Stream::READABLE_EPS);
  ps.set_current_filling((CGAL::FILLING)filling);
  
  Point3 p(-400,300,0);
  Point3 q(400,300,0);
  Point3 r(-200,0,100);Point3 s(200,100,0);Line3 line(r,s);
  
  Point3 t(0,100,200);Point3 u(300,100,200);Point3 v(150,300,200);Triangle3 triangle(t,u,v);
  
  Point3 w(-100,100,-100);Point3 m(0,300,-200);Segment3 segment(w,m);
  
  Point3 a(-200,-300,200);Point3 b(-100,-400,400);
  Point3 c(300,-200,200);Point3 d(100,0,300);
  Tetrahedron tetrahedre(a,b,c,d);
  
  ps.set_current_filling(CGAL::UNIFORM_FILL);
  ps.set_fill_color(CGAL::BLACK);
  ps.set_border_color(CGAL::RED);
  ps << p;
  
  ps.set_current_filling(CGAL::UNIFORM_FILL);
  ps.set_border_color(CGAL::BLUE);
  ps.set_fill_color(CGAL::RED);
  ps << q;
  
  ps.set_current_filling(CGAL::NO_FILL);
  ps.set_border_color(CGAL::VIOLET);
  ps.set_fill_color(CGAL::BLUE);
  ps << line;
 
  ps.set_current_filling(CGAL::UNIFORM_FILL);
  ps.set_border_color(CGAL::BLACK);
  ps.set_fill_color(CGAL::ORANGE);
  ps << triangle;
  
  ps.set_current_filling(CGAL::WIRED_CULLBACK_FACING);
  ps.set_border_color(CGAL::BLACK);
  ps.set_fill_color(CGAL::ORANGE);
  ps << segment;
  
  ps.set_current_filling(CGAL::NORMAL_FILL);
  ps.set_fill_color(CGAL::GRAY);
  ps.set_border_color(CGAL::GREEN);
  ps << tetrahedre;

  struct timeval * t1 = (struct timeval *)malloc(sizeof(struct timeval));
  struct timezone *t1z = (struct timezone *)malloc(sizeof(struct
							timezone));
  struct timeval * t2 = (struct timeval *)malloc(sizeof(struct timeval));
  struct timezone *t2z = (struct timezone *)malloc(sizeof(struct
							timezone));
  gettimeofday(t1,t1z);
  ps.display();
  gettimeofday(t2,t2z);

  printf("Temps de calcul en secondes pour display() : %ld sec.\n", t2->tv_sec - t1->tv_sec);
  
  std::cerr << "Creation du fichier : " << filename  << std::endl;
  
  return 0;
}
