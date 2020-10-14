#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Octree.h>
#include <CGAL/Octree/IO.h>
#include <CGAL/Random.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef std::vector<Point_2> Point_vector;

typedef CGAL::Octree::Octree<Point_vector> Quadtree;

int main(int argc, char **argv)
{
  CGAL::Random r;

  Point_vector points_2d;
  for (std::size_t i = 0; i < 5; ++ i)
    points_2d.emplace_back(r.get_double(-1., 1.),
                           r.get_double(-1., 1.));

  Quadtree quadtree(points_2d);
  quadtree.refine(10, 1);
  std::cerr << "Quadtree = " << std::endl
            << quadtree << std::endl;

  return EXIT_SUCCESS;
}
