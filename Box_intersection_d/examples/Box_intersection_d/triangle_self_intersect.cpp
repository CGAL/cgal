#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/box_intersection_d.h>
#include <vector>
#include <fstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef Kernel::Point_3                                       Point_3;
typedef Kernel::Triangle_3                                    Triangle_3;
typedef std::vector<Triangle_3>                               Triangles;
typedef Triangles::iterator                                   Iterator;
typedef CGAL::Box_intersection_d::Box_with_handle_d<double,3,Iterator> Box;

struct Report {
  Triangles* triangles;

  Report(Triangles& triangles)
    : triangles(&triangles)
  {}

  // callback functor that reports all truly intersecting triangles
  void operator()(const Box& a, const Box& b) const 
  {
    std::cout << "Box " << (a.handle() - triangles->begin()) << " and "
              << (b.handle() - triangles->begin()) << " intersect";
    if ( ! a.handle()->is_degenerate() && ! b.handle()->is_degenerate()
         && CGAL::do_intersect( *(a.handle()), *(b.handle()))) {
    std::cout << ", and the triangles intersect also";
    }
    std::cout << '.' << std::endl;
  }
};


int main(int argc, char*argv[])
{
  std::ifstream in((argc>1)?argv[1]:"data/triangles.xyz");
  Triangles triangles;
  Triangle_3 t;
  while(in >> t){
    triangles.push_back(t);
  }
  
  // Create the corresponding vector of bounding boxes
  std::vector<Box> boxes;
  for ( Iterator i = triangles.begin(); i != triangles.end(); ++i)
    boxes.push_back( Box( i->bbox(), i));
  
  // Run the self intersection algorithm with all defaults
  CGAL::box_self_intersection_d( boxes.begin(), boxes.end(), Report(triangles));
  return 0;
}
