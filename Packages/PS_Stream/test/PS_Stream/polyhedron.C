#include <CGAL/Cartesian.h>

#include <sys/time.h>
#include <fstream>
#include <iostream>

#include <CGAL/IO/Color.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/Halfedge_data_structure_polyhedron_default_3.h>
#include <CGAL/Polyhedron_default_traits_3.h>
#include <CGAL/Polyhedron_3.h>

#include "PS_Stream_3.C"

typedef CGAL::Cartesian<double> D;

typedef CGAL::Bbox_3     PS_BBox3;
typedef D::Point_3       Point3;
typedef D::Direction_3   Direction;

typedef CGAL::Halfedge_data_structure_polyhedron_default_3<D> HDS;
typedef CGAL::Polyhedron_default_traits_3<D> Traits;
typedef CGAL::Polyhedron_3<Traits,HDS> Polyhedron;

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
  
  std::cerr << "Creation du fichier : " << filename  << std::endl;
  
  return 0;
}
