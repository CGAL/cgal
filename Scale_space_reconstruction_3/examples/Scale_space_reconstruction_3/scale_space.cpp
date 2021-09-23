#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/Timer.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;
typedef Kernel::Point_3                                         Point;

typedef CGAL::Scale_space_surface_reconstruction_3<Kernel>      Reconstruction;

typedef Reconstruction::Facet_const_iterator                    Facet_iterator;

int main(int argc, char** argv)
{
  if (argc != 2)
  {
    std::cerr << "Error, no input file provided\n";
    return 1;
  }

  std::cerr << "Reading " << std::flush;
  std::vector<Point> points;
  if(!CGAL::IO::read_points(argv[1], std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }
  std::cerr << "done: " << points.size() << " points." << std::endl;

  std::cerr << "Reconstruction ";
  CGAL::Timer t;
  t.start();

  // Construct the mesh in a scale space.
  Reconstruction reconstruct (points.begin(), points.end());
  reconstruct.increase_scale(4);
  reconstruct.reconstruct_surface();
  std::cerr << "done in " << t.time() << " sec." << std::endl;
  t.reset();
  std::ofstream out ("out.off");
  out << reconstruct;
  std::cerr << "Writing result in " << t.time() << " sec." << std::endl;
  std::cerr << "Done." << std::endl;
  return EXIT_SUCCESS;
}
