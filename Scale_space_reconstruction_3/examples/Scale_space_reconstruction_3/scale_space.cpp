
#include <fstream>
#include <iostream>

#include <CGAL/Scale_space_surface_reconstruction_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/IO/read_off_points.h>
#include <CGAL/Timer.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel     Kernel;

typedef CGAL::Scale_space_surface_reconstruction_3< Kernel >    Reconstruction;

typedef Reconstruction::Point                                   Point;
typedef std::vector< Point >                                    Point_collection;

typedef Reconstruction::Triple_const_iterator                   Triple_iterator;
typedef CGAL::Timer Timer;

int main(int argc, char* argv[]) {
    if (argc!=2){
      std::cerr << "Error, no input file provided\n";
      return 1;
    }
    // Read the data.
    Point_collection points;
    std::ifstream in(argv[1]);
    std::cerr << "Reading " << std::flush;
    if( !in || !CGAL::read_off_points( in, std::back_inserter( points ) ) ) {
        std::cerr << "Error: cannot read file" << std::endl;
        return EXIT_FAILURE;
    }
    std::cerr << "done: " << points.size() << " points." << std::endl;

    Timer t;
    t.start();
	// Construct the mesh in a scale space.
	Reconstruction reconstruct( 10, 200 );
	reconstruct.reconstruct_surface( points.begin(), points.end(), 4 );
        std::cerr << "Reconstruction done in " << t.time() << " sec." << std::endl;
        t.reset();
    std::ofstream out ("out.off");
    // Write the reconstruction.
    std::cerr << "Neighborhood radius^2 = " << reconstruct.neighborhood_squared_radius() << std::endl;
    for( std::size_t shell = 0; shell < reconstruct.number_of_shells(); ++shell ) {
        std::cerr << "Shell " << shell << std::endl;
        for( Triple_iterator it = reconstruct.shell_begin( shell ); it != reconstruct.shell_end( shell ); ++it )
          out << "3 "<< *it << '\n'; // We write a '3' in front so that it can be assembled into an OFF file
    }
        std::cerr << "Writing result in " << t.time() << " sec." << std::endl;
    std::cerr << "Done." << std::endl;
    return EXIT_SUCCESS;
}
