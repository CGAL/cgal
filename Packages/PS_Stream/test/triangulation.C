#include <sys/time.h>
#include <fstream>
#include <iostream>

#include <CGAL/Triangulation_2.h> 
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Triangulation_euclidean_traits_xy_3.h>

#include <CGAL/IO/Color.h>
#include <CGAL/Point_3.h>
#include <CGAL/Direction_3.h>

#include <CGAL/Bbox_3.h>

#include "PS_Stream_3.C"


using namespace std;

typedef CGAL::Cartesian<double> D;
//typedef CGAL::Cartesian<leda_real> R;

typedef CGAL::Bbox_3 PS_BBox3;
typedef CGAL::Direction_3< D > Direction;
typedef CGAL::Point_3< D > Point3;

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
  
  Triangulation tr;
  cout << "Lecture du fichier ./data/terrain_2.dat" << endl;
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
  
  cerr << "Creation du fichier : " << filename  << endl;
  
return 0;
}
