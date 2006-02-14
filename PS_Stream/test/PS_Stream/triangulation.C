#include <CGAL/Cartesian.h>

#include <sys/time.h>
#include <fstream>
#include <iostream>

#include <CGAL/Triangulation_2.h> 
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <CGAL/Triangulation_euclidean_traits_xy_3.h>

#include <CGAL/IO/Color.h>

#include <CGAL/Bbox_3.h>

#include "PS_Stream_3.C"

typedef CGAL::Cartesian<double> D;

typedef CGAL::Bbox_3    PS_BBox3;
typedef D::Direction_3  Direction;
typedef D::Point_3      Point3;

typedef double Nt;
typedef CGAL::Cartesian<Nt> Rp;
typedef CGAL::Triangulation_euclidean_traits_xy_3<Rp>  Gt;
typedef Gt::Point    Point ;
typedef Gt::Segment  Segment;
typedef Gt::Triangle Triangle;
typedef CGAL::Triangulation_vertex_base_2<Gt> Vb; 
typedef CGAL::Triangulation_face_base_2<Gt>   Fb; 
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds> Triangulation;

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
  
  Triangulation tr;
  std::cout << "Lecture du fichier ./data/terrain_2.dat" << std::endl;
  std::ifstream is("./data/terrain_2.dat");
  CGAL::set_ascii_mode(is); 
  std::istream_iterator<Point> it(is);
  std::istream_iterator<Point> end;
  
  for( ; it != end; it++) {
    tr.insert( *it);
  }
  
  ps << tr;

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
