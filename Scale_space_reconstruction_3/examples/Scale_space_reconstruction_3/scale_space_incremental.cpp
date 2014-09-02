#include <fstream>
#include <iostream>

#include <CGAL/Scale_space_surface_reconstruction_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/read_off_points.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >    Reconstruction;

typedef Reconstruction::Point                                   Point;
typedef std::vector< Point >                                    Pointset;

typedef Reconstruction::Const_triple_iterator                   Triple_iterator;

int main(void) {
    // Read the data.
	Pointset points;
	std::ifstream in("kitten.off");
    std::cout << "Reading " << std::flush;
    if( !in || !CGAL::read_off_points( in, std::back_inserter( points ) ) ) {
        std::cerr << "Error: cannot read file" << std::endl;
        return EXIT_FAILURE;
    }
	std::cout << "done: " << points.size() << " points." << std::endl;
    
	// Construct the reconstruction to estimate the neighborhood radius.
    // The -1 indicates the neighborhood radius is not yet set.
	Reconstruction reconstruct( -1 );

    // Add the points.
    reconstruct.insert( points.begin(), points.end() );

    // Advance the scale-space several steps.
    // This automatically estimates the scale-space.
    reconstruct.advance_scale_space( 2 );

    // Re-estimate the neighborhood radius.
    reconstruct.estimate_neighborhood_radius( 10, 100 );

    // Advance the scale-space further.
    reconstruct.advance_scale_space( 2 );

    // Manually set the neighborhood radius.
    reconstruct.set_neighborhood_squared_radius( 0.05 );

    // Reconstruct the surface from the current scale-space.
    reconstruct.reconstruct_surface();
    std::cout << "Reconstruction done:" << std::endl;

    // Write the reconstruction.
    std::cout << "Neighborhood radius^2 = " << reconstruct.neighborhood_squared_radius() << std::endl;
    for( Triple_iterator it = reconstruct.surface_begin(); it != reconstruct.surface_end(); ++it )
        std::cout << *it << std::endl;

	std::cout << "Done." << std::endl;
    return EXIT_SUCCESS;
}