
#include <CGAL/Simple_cartesian.h>
#include <CGAL/bounding_box.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Random.h>

#include <iostream>
#include <fstream>

typedef CGAL::Simple_cartesian<double>::Point_3 Point_3;

int main(int argc, char* argv[])
{
#if 0
  std::cerr << "ascii to binary\n";
  int d;
  int N;
  Point_3 p;
  
  std::ifstream ascii(argv[1]);
  ascii >> d >> N; 
  std::ofstream binary(argv[2], std::ios::out | std::ios::binary);
  CGAL::set_binary_mode(binary);
  CGAL::write(binary, d);
  CGAL::write(binary, N);
  for(int i=0; i < N; i++){
    ascii >> p;
    binary << p;
  }

  #elif 1
  std::cerr << " binary bbox\n";
  int d=0;
  int N=0;
  Point_3 p;
  
  std::ifstream binary(argv[1], std::ios::in | std::ios::binary);
  std::ofstream bbox(argv[2], std::ios::out | std::ios::binary);
  CGAL::set_binary_mode(binary);
  CGAL::set_binary_mode(bbox);
  CGAL::read(binary,d);
  CGAL::read(binary,N);
  CGAL::write(bbox, d);
  CGAL::write(bbox, N);

  std::vector<Point_3> points;

  while(binary >> p){
    points.push_back(p);
  }

  CGAL::Bbox_3 bb = CGAL::bounding_box(points.begin(), points.end()).bbox();
  double dx = bb.xmax() - bb.xmin();
  double dy = bb.ymax() - bb.ymin();
  double dz = bb.zmax() - bb.zmin();

  for(int i=0; i < N; i++){
    double  rx  = CGAL::default_random.get_double();
    double  ry  = CGAL::default_random.get_double();
    double  rz  = CGAL::default_random.get_double();
    bbox << Point_3(bb.xmin() + dx * rx , bb.ymin() + dy * ry , bb.zmin() + dz * rz);
    /*bbox << bb.xmin() + dx * rx << " ";
    bbox << bb.ymin() + dy * ry << " ";
    bbox << bb.zmin() + dz * rz << std::endl;*/
  }
  bbox.close();
  binary.close();

#else
  std::cerr << " reading\n";
  int d=0;
  int N=0;
  Point_3 p;
  
  std::ifstream binary(argv[1], std::ios::in | std::ios::binary);
  std::ofstream ascii(argv[2]);
  CGAL::set_binary_mode(binary);
  CGAL::read(binary,d);
  CGAL::read(binary,N);
  ascii << d << std::endl << N << std::endl;
  for(int i=0; i < N; i++){
    binary >> p;
    ascii << p << std::endl;
  }
  ascii.close();
  binary.close();
  
#endif
  return 0;
}
