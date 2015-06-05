#include <CGAL/IO/STL_reader.h>
#include <CGAL/array.h>

#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char* argv[])
{
  std::ifstream input((argc>1)?argv[1]:"data/triangle.stl",
                      std::ios::in | std::ios::binary);

  std::vector< CGAL::cpp11::array<double,3> > points;
  std::vector< CGAL::cpp11::array<int,3> > faces;

  CGAL::read_STL( input,
                  points,
                  faces,
                  true);
  
  std::cout.precision(17);
  std::cout << "OFF\n" << points.size() << " " << faces.size()  << " 0" << std::endl;
  for(int i=0; i < points.size(); i++){
    std::cout << points[i][0] << " " << points[i][1] << " " << points[i][2]<< std::endl;
  }

  for(int i=0; i < faces.size(); i++){
    std::cout << "3 " << faces[i][0] << " " << faces[i][1] << " " << faces[i][2] << std::endl;
  }

  return 0;
}
