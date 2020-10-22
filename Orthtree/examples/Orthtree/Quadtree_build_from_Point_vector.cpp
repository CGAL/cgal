#include <fstream>
#include <iostream>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Quadtree.h>
#include <CGAL/Random.h>

// Type Declarations
typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::Point_2 Point_2;
typedef std::vector<Point_2> Point_vector;

typedef CGAL::Quadtree<Kernel, Point_vector> Quadtree;

int main()
{
  CGAL::Random r;

  Point_vector points_2d;
  for (std::size_t i = 0; i < 5; ++ i)
    points_2d.emplace_back(r.get_double(-1., 1.),
                           r.get_double(-1., 1.));

  Quadtree quadtree(points_2d);
  quadtree.refine(10, 5);


  std::ofstream opoints ("points.xyz");
  opoints.precision(18);
  for (const auto& p : points_2d)
    opoints << p << " 0" << std::endl;


  {
    std::ofstream ofile ("quadtree.polylines.txt");
    ofile.precision(18);
    quadtree.dump_to_polylines(ofile);
  }

  quadtree.grade();

  {
    std::ofstream ofile ("quadtree_graded.polylines.txt");
    ofile.precision(18);
    quadtree.dump_to_polylines(ofile);
  }

  return EXIT_SUCCESS;
}
