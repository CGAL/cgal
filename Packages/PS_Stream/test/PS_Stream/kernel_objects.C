
#include <sys/time.h>
#include <iostream>

#include <CGAL/Cartesian.h>
#include <CGAL/IO/Color.h>
#include <CGAL/Point_3.h>
#include <CGAL/Line_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Bbox_3.h>
#include <CGAL/Segment_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangle_3.h>

#include "PS_Stream_3.C"


using namespace std;

typedef CGAL::Cartesian<double> D;
//typedef CGAL::Cartesian<leda_real> R;

typedef CGAL::Bbox_3 PS_BBox3;

typedef CGAL::Direction_3< D > Direction;
typedef CGAL::Point_3< D > Point3;
typedef CGAL::Point_2< D > Point2;
typedef CGAL::Plane_3< D > Plane3;
typedef CGAL::Line_3< D > Line3;
typedef CGAL::Segment_3< D > Segment3;
typedef CGAL::Triangle_3< D > Triangle3;
typedef CGAL::Tetrahedron_3< D > Tetrahedron;

int main(void) {
  
  PS_BBox3 bb3(-500,-500,-500,500,500,500);
  double x,y,z,lx,ly,lz;
  int filling;
  
  cerr << "Filling 0,1,2,3 : ";cin >> filling;
  std::string filename;
  cerr << "Enter file name: ";cin >> filename;
  cerr << "Enter view x: ";cin >> x;
  cerr << "Enter view y: ";cin >> y;
  cerr << "Enter view z: ";cin >> z;
  
  Direction dir(x,y,z);
  cerr << "Enter light x: ";cin >> lx;
  cerr << "Enter light y: ";cin >> ly;
  cerr << "Enter light z: ";cin >> lz;
  
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
  
  cerr << "Creation du fichier : " << filename  << endl;
  
return 0;
}
