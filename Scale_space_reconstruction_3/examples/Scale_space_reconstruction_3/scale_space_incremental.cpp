#include <fstream>
#include <iostream>
#include <algorithm>

#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >                Reconstruction;
typedef CGAL::Scale_space_reconstruction_3::Weighted_PCA_smoother< Kernel > Smoother;
typedef CGAL::Scale_space_reconstruction_3::Alpha_shape_mesher< Kernel >    Mesher;

typedef Reconstruction::Point                                   Point;

typedef Reconstruction::Facet_const_iterator                   Facet_iterator;

// function for writing the reconstruction output in the off format
void dump_reconstruction(const Reconstruction& reconstruct, std::string name)
{
  std::ofstream output(name.c_str());
  output << "OFF " << reconstruct.number_of_points() << " "
         << reconstruct.number_of_facets() << " 0\n";

  std::copy(reconstruct.points_begin(),
            reconstruct.points_end(),
            std::ostream_iterator<Point>(output,"\n"));
  for( Facet_iterator it = reconstruct.facets_begin(); it != reconstruct.facets_end(); ++it )
      output << "3 " << *it << std::endl;
}

int main(int argc, char* argv[]) {
    // Read the data.
    std::vector<Point> points;
    if (argc!=2){
      std::cerr << "Error, no input file provided\n";
      return 1;
    }
    std::ifstream in(argv[1]);
    std::cout << "Reading " << std::flush;
    if( !in || !CGAL::read_off_points( in, std::back_inserter( points ) ) ) {
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
        // Write the reconstruction.
        dump_reconstruction(reconstruct, "reconstruction1.off");
      }
      else
      {
        std::cout << "Second reconstruction done." << std::endl;
        // Write the reconstruction.
        dump_reconstruction(reconstruct, "reconstruction2.off");
      }
    }

    std::cout << "Reconstructions are ready to be examinated in your favorite viewer" << std::endl;
    return EXIT_SUCCESS;
}
