#include <sys/time.h>
#include <fstream>
#include <iostream>

#include <CGAL/Cartesian.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Point_3.h>
#include <CGAL/Direction_3.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>
"
#include "PS_Stream_3.C"

using namespace std;

typedef CGAL::Cartesian<double> D;
//typedef CGAL::Cartesian<leda_real> R;

typedef CGAL::Bbox_3 PS_BBox3;
typedef CGAL::Point_3< D > Point3;
typedef CGAL::Direction_3< D > Direction;

typedef CGAL::Halfedge_data_structure_polyhedron_default_3<D> HDS;
typedef  CGAL::Polyhedron_default_traits_3<D> Traits;
typedef  CGAL::Polyhedron_3<Traits,HDS> Polyhedron;

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
  ps.set_border_color(CGAL::RED);ps.set_fill_color(CGAL::BLUE);
  Point3 p1( 100.0, 50.0, 0.0);
  Point3 q1( 50.0, 100.0, 0.0);
  Point3 r1( 10.0, 0.0, 100.0);
  Point3 s1( 0.0, 200.0, 0.0);
  Polyhedron P;  
  P.make_tetrahedron(p1,q1,r1,s1);  
  ps << P;

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
