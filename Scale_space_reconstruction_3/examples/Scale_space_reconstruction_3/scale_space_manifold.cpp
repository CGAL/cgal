#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Scale_space_surface_reconstruction_3.h>

#include <CGAL/IO/read_points.h>
#include <CGAL/Timer.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >                Reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother< Kernel > Smoother;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher< Kernel >    Mesher;

typedef Reconstruction::Point                                   Point;

typedef Reconstruction::Facet_const_iterator                    Facet_iterator;
typedef Mesher::Facet_const_iterator                            Mesher_iterator;
typedef CGAL::Timer Timer;

int main(int argc, char* argv[]) {
  // Read the data
  std::string fname = argc==1?CGAL::data_file_path("points_3/kitten.off"):argv[1];

  std::cerr << "Reading " << std::flush;
  std::vector<Point> points;
  if(!CGAL::IO::read_points(fname, std::back_inserter(points)))
  {
    std::cerr << "Error: cannot read file" << std::endl;
    return EXIT_FAILURE;
  }

  std::cerr << "done: " << points.size() << " points." << std::endl;

  Timer t;
  t.start();

  // Construct the mesh in a scale space.
  Reconstruction reconstruct(points.begin(), points.end() );
  Smoother smoother(10, 200 );
  reconstruct.increase_scale(4, smoother);

  Mesher mesher(smoother.squared_radius(),
                false, // Do not separate shells
                true // Force manifold output
               );
  reconstruct.reconstruct_surface(mesher);

  std::cerr << "Reconstruction done in " << t.time() << " sec." << std::endl;

  std::ofstream out("out.off");
  // Write the reconstruction.
  for(Facet_iterator it = reconstruct.facets_begin(); it != reconstruct.facets_end(); ++it )
    // We write a '3' in front so that it can be assembled into an OFF file
    out << "3 " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << '\n';
  out.close();

  std::ofstream garbage("garbage.off");
  // Write facets that were removed to force manifold output
  for(Mesher_iterator it = mesher.garbage_begin(); it != mesher.garbage_end(); ++it )
    garbage << "3 " << (*it)[0] << " " << (*it)[1] << " " << (*it)[2] << '\n';

  std::cerr << "Done." << std::endl;

  return EXIT_SUCCESS;
}
