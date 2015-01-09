#include <fstream>
#include <iostream>
#include <algorithm>

#include <CGAL/Scale_space_surface_reconstruction_3.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/IO/read_off_points.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >    Reconstruction;

typedef Reconstruction::Point                                   Point;
typedef std::vector< Point >                                    Point_collection;

typedef Reconstruction::Triple_const_iterator                   Triple_iterator;

// function for writing the reconstruction output in the off format
void dump_reconstruction(const Reconstruction& reconstruct, std::string name)
{
  std::ofstream output(name.c_str());
  output << "OFF " << reconstruct.number_of_points() << " "
         << reconstruct.number_of_triangles() << " 0\n";

  std::copy(reconstruct.points_begin(),
            reconstruct.points_end(),
            std::ostream_iterator<Point>(output,"\n"));
  for( Triple_iterator it = reconstruct.surface_begin(); it != reconstruct.surface_end(); ++it )
      output << "3 " << *it << std::endl;
}

int main(int argc, char* argv[]) {
    // Read the data.
    Point_collection points;
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

    // Construct the reconstruction with parameters for
    // the neighborhood squared radius estimation.
    Reconstruction reconstruct( 10, 100 );

    // Add the points.
    reconstruct.insert( points.begin(), points.end() );

    // Advance the scale-space several steps.
    // This automatically estimates the scale-space.
    reconstruct.increase_scale( 2 );

    // Reconstruct the surface from the current scale-space.
    std::cout << "Neighborhood squared radius is "
              << reconstruct.neighborhood_squared_radius() << std::endl;
    reconstruct.reconstruct_surface();
    std::cout << "First reconstruction done." << std::endl;

    // Write the reconstruction.
    dump_reconstruction(reconstruct, "reconstruction1.off");

    // Advancing the scale-space further and visually compare the reconstruction result
    reconstruct.increase_scale( 2 );

    // Reconstruct the surface from the current scale-space.
    std::cout << "Neighborhood squared radius is "
              << reconstruct.neighborhood_squared_radius() << std::endl;
    reconstruct.reconstruct_surface();
    std::cout << "Second reconstruction done." << std::endl;

    // Write the reconstruction.
    dump_reconstruction(reconstruct, "reconstruction2.off");

    std::cout << "Reconstructions are ready to be examinated in your favorite viewer" << std::endl;
    return EXIT_SUCCESS;
}
