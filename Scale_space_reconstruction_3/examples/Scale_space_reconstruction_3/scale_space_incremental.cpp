#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/IO/OFF.h>

#include <algorithm>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >                Reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother< Kernel > Smoother;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher< Kernel >    Mesher;

typedef Reconstruction::Point                                   Point;

typedef Reconstruction::Facet_const_iterator                   Facet_iterator;


int main(int argc, char* argv[])
{
    // Read the data.
    std::string fname = argc==1?CGAL::data_file_path("points_3/kitten.off"):argv[1];
    std::cout << "Reading " << std::flush;
    std::vector<Point> points;
    if(!CGAL::IO::read_points(fname, std::back_inserter(points)))
    {
      std::cerr << "Error: cannot read file" << std::endl;
      return EXIT_FAILURE;
    }

    std::cout << "done: " << points.size() << " points." << std::endl;

    // Construct the reconstruction
    Reconstruction reconstruct;

    // Add the points.
    reconstruct.insert( points.begin(), points.end() );

    // Two passes
    for (std::size_t i = 0; i < 2; ++ i)
    {
      // Construct the smoother with parameters for
      // the neighborhood squared radius estimation.
      Smoother smoother( 10, 100 );

      // Advance the scale-space several steps.
      // This automatically estimates the scale-space.
      reconstruct.increase_scale( 2, smoother );

      // Reconstruct the surface from the current scale-space.
      std::cout << "Neighborhood squared radius is "
                << smoother.squared_radius() << std::endl;

      Mesher mesher (smoother.squared_radius());
      reconstruct.reconstruct_surface(mesher);
      if (i == 0)
      {
        std::cout << "First reconstruction done." << std::endl;
        CGAL::IO::write_OFF("reconstruction1.off",
                            reconstruct.points(),
                            reconstruct.facets(),
                            CGAL::parameters::stream_precision(17));
      }
      else
      {
        std::cout << "Second reconstruction done." << std::endl;
        CGAL::IO::write_OFF("reconstruction2.off",
                            reconstruct.points(),
                            reconstruct.facets(),
                            CGAL::parameters::stream_precision(17));
      }
    }

    std::cout << "Reconstructions are ready to be examined in your favorite viewer" << std::endl;
    return EXIT_SUCCESS;
}
