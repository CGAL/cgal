#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/box_intersection_d.h>
#include <iostream>
#include <fstream>
#include <iterator>
#include <sstream>

#include <boost/lexical_cast.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polygon_2<Kernel> Polygon_2;
typedef Kernel::Point_2 Point_2;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Segment_2 Segment_2;
typedef Kernel::Direction_2 Direction_2;

typedef std::pair<Point_2,std::string> Point;
typedef std::list<Point>::iterator iterator;


void contour(int count, int n, int m, int dim)
{
  if(count == m){
    std::cout << n << std::endl;
  }
 
  std::cin.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
  std::list<std::pair<Point_2,std::string> > points;

  for(int i = 0; i <n; i++){
    std::string line;
    std::getline(std::cin, line);
    if(count == m){
      std::cout << line << std::endl;
    }
  }
}

int main(int argc, char* argv[])
{
  int dim = 2; // the dimension where the coordinates are constant
  int m = 0;
  if(argc > 1){
    dim = boost::lexical_cast<int>(argv[1]);
  } 
  if(argc > 2){
    m = boost::lexical_cast<int>(argv[2]);
  }
  int n;
  int count = 0;
  while(std::cin >> n){
    contour(count, n, m, dim);
    count++;
  }
  std::cerr << "done" << std::endl;
  return 0;
}
