// polygonal_surface_reconstruction_test.cpp
//
// Test the Polygonal Surface Reconstruction method with different 
//    - kernels(Simple_cartesian, EPICK)
//    - solvers(GLPK, SCIP)
//    - use/ignore provided planar segmentation
//    - input file formats (pwn, ply). For ply format, a property "segment_index"
//      must be present storing the plane index for each point(-1 if the point is
//      not assigned to a plane).

#include <iostream>

#ifdef CGAL_USE_GLPK
#include <CGAL/GLPK_mixed_integer_program_traits.h>
typedef CGAL::GLPK_mixed_integer_program_traits<double>				GLPK_Solver;
#endif

#ifdef CGAL_USE_SCIP
#include <CGAL/SCIP_mixed_integer_program_traits.h>
typedef CGAL::SCIP_mixed_integer_program_traits<double>				SCIP_Solver;
#endif

#if defined(CGAL_USE_GLPK) || defined(CGAL_USE_SCIP)


#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include "polygonal_surface_reconstruction_test_framework.h"

// kernels:
typedef CGAL::Simple_cartesian<double>								Cartesian;
typedef CGAL::Exact_predicates_inexact_constructions_kernel			Epick;


int main(int argc, char* argv[])
{
    std::cerr << "Testing the Polygonal Surface Reconstruction method...\n";

#if 0
    argc = 2;
    argv[1] = "data/icosahedron.pwn";
    //argv[1] = "data/house.ply";
    //argv[1] = "data/chair.ply";
    //argv[1] = "data/square.ply";
#endif

	// usage
	if (argc - 1 == 0) {
		std::cerr << "For the input point cloud, reconstruct a water-tight polygonal surface.\n";
		std::cerr << "\n";
		std::cerr << "Usage: " << argv[0] << " point_cloud_file" << std::endl;
        std::cerr << "Input file formats are \'pwn\' and \'ply\'. No output.\n";
		return EXIT_FAILURE;
	}

	char* input_file = argv[1];

	//---------------------------------------------------------------------

    std::cerr << "--- Using Simple cartesian kernel";

	//---------------------------------------------------------------------

#ifdef CGAL_USE_GLPK
    std::cerr << "\n\t---- Using GLPK solver\n";

	std::cerr << "\t\t---- using provided planes\n";
	reconstruct<Cartesian, GLPK_Solver>(input_file, false);

    std::cerr << "\t\t---- re-extract planes\n";
	reconstruct<Cartesian, GLPK_Solver>(input_file, true);
#endif

#ifdef CGAL_USE_SCIP
    std::cerr << "\n\t---- Using SCIP solver\n";

	std::cerr << "\t\t---- using provided planes\n";
	reconstruct<Cartesian, SCIP_Solver>(input_file, false);

    std::cerr << "\t\t---- re-extract planes\n";
	reconstruct<Cartesian, SCIP_Solver>(input_file, true);
#endif

	//---------------------------------------------------------------------

    std::cerr << "\n--- Using Epick kernel";

	//---------------------------------------------------------------------

#ifdef CGAL_USE_GLPK
    std::cerr << "\n\t---- Using GLPK solver\n";

	std::cerr << "\t\t---- using provided planes\n";
	reconstruct<Epick, GLPK_Solver>(input_file, false);

    std::cerr << "\t\t---- re-extract planes\n";
	reconstruct<Epick, GLPK_Solver>(input_file, true);
#endif 

#ifdef CGAL_USE_SCIP
    std::cerr << "\n\t---- Using SCIP solver\n";

	std::cerr << "\t\t---- using provided planes\n";
	reconstruct<Epick, SCIP_Solver>(input_file, false);

    std::cerr << "\t\t---- re-extract planes\n";
	reconstruct<Epick, SCIP_Solver>(input_file, true);
#endif

    return EXIT_SUCCESS;
}

#else

int main(int argc, char* argv[])
{
    std::cerr << "This test requires either GLPK or SCIP.\n";
    return EXIT_FAILURE;
}

#endif  // defined(CGAL_USE_GLPK) || defined(CGAL_USE_SCIP)
