
#include <fstream>
#include <iostream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Point_set_3.h>
#include <CGAL/Point_set_3/IO.h>
#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Scale_space_reconstruction_3/Advancing_front_mesher.h>
#include <CGAL/Scale_space_reconstruction_3/Jet_smoother.h>

#include <CGAL/IO/read_off_points.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3<Kernel>                    Reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Advancing_front_mesher<Kernel>    Mesher;
typedef CGAL::Scale_space_reconstruction_3::Jet_smoother<Kernel>              Smoother;

typedef Kernel::Point_3 Point;
typedef CGAL::Point_set_3<Point> Point_set;

typedef Reconstruction::Facet_const_iterator                   Facet_iterator;

int main(int argc, char** argv)
{
  if (argc!=2){
    std::cerr << "Error, no input file provided\n";
    return 1;
  }
  // Read the data.
  Point_set points;
  std::ifstream in(argv[1]);
  std::cerr << "Reading " << std::flush;
  in >> points;

  if (points.empty())
  {
    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }
  std::cerr << "done: " << points.size() << " points." << std::endl;

  std::cerr << "Reconstruction ";

  CGAL::Timer t;
  t.start();

  // Construct the mesh in a scale space.
  Reconstruction reconstruct (points.points().begin(), points.points().end());
  reconstruct.increase_scale<Smoother> (4);
  reconstruct.reconstruct_surface (Mesher (0.5));

  std::cerr << "done in " << t.time() << " sec." << std::endl;

  t.reset();

  std::ofstream out ("out.off");
  out << "OFF" << std::endl << points.size() << " " << reconstruct.number_of_facets() << " 0" << std::endl;

  for (Point_set::iterator it = points.begin(); it != points.end(); ++ it)
    out << points.point(*it) << std::endl;

  for (Reconstruction::Facet_iterator it = reconstruct.facets_begin();
       it != reconstruct.facets_end(); ++ it)
    out << "3 " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << std::endl;

  std::cerr << "Writing result in " << t.time() << " sec." << std::endl;

  std::cerr << "Done." << std::endl;

  return EXIT_SUCCESS;
}
