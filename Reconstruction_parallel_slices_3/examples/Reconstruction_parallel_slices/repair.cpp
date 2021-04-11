#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/box_intersection_d.h>
#include <CGAL/property_map.h>
#include <CGAL/Reconstruction_from_parallel_slices_3/check_and_fix_input.h>
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

typedef std::pair<iterator,iterator> Segment;
typedef std::vector<Segment>::iterator segment_iterator;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,segment_iterator> Box;

struct Get_point_from_pair{
  //classical typedefs
  typedef const std::pair<Point_2, std::string>& key_type;
  typedef Point_2 value_type;
  typedef const value_type& reference;
  typedef boost::readable_property_map_tag category;
  friend reference get(const Get_point_from_pair&, key_type p) {
    return p.first;
  }
};

void contour(int count, int n, int dim)
{
  std::cerr << "\ncontour  " << count << " with " << n << " vertices" << std::endl;

  std::cin.ignore((std::numeric_limits<std::streamsize>::max)(), '\n');
  std::list<std::pair<Point_2,std::string> > points;

  for(int i = 0; i <n; i++){
    std::string line;
    std::getline(std::cin, line);
    std::istringstream iss(line);

    double x,y,z;
    iss >> x >> y >> z;

    Point_2 p = (dim==0) ? Point_2(y,z) : (dim==1) ? Point_2(x,z) : Point_2(x,y);
    if(! points.empty()){
      if(p != points.back().first){
      points.push_back(std::make_pair(p, line));
      } else {
        std::cerr << "identical consecutive points " << p << std::endl;
      }
    } else {
     points.push_back(std::make_pair(p, line));
    }
  }
  if(points.size() <= 3){
    std::cerr << "ignore segment" <<  std::endl;
    std::cerr << points.front().first << " -- " << (++points.begin())->first << std::endl;
    return;
  }

  Get_point_from_pair ppmap;

  if( CGAL::Reconstruction_from_parallel_slices::
        find_safe_start_to_fix_consecutive_overlapping_segments(points, ppmap) )
  {
    CGAL::Reconstruction_from_parallel_slices::
      fix_consecutive_overlapping_segments(points, ppmap);

    if(points.size() <= 3){
      std::cerr << "ignore segment" <<  std::endl;
      std::cerr << points.front().first << " -- " << get(ppmap, *(++points.begin())) << std::endl;
      return;
    }

    if ( !CGAL::Reconstruction_from_parallel_slices::
           fix_self_intersections(points, ppmap) )
      std::cerr << "Could not fix the self-intersection\n";
  } else {
    std::cerr << "No safe start point found" << std::endl;
  }

  std::cout << points.size() <<std::endl;
  for(iterator it = points.begin(); it!= points.end(); ++it){
    std::cout << it->second << std::endl;
  }
}

int main(int argc, char** argv)
{
  if (argc!=2)
  {
    std::cerr << "Please provide the constant coordinate index (0,1,2)\n";
    return 1;
  }
  int constant_coordinate = boost::lexical_cast<int>(argv[1]);
  int n;
  int count = 0;
  while(std::cin >> n){
    contour(count, n, constant_coordinate);
    count++;
  }
  std::cerr << "done" << std::endl;
  return 0;
}
