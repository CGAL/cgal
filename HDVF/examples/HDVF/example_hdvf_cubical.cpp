#include <iostream>
#include <fstream>
#include <CGAL/HDVF/Cub_object_io.h>
#include <CGAL/HDVF/Cubical_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

typedef int CoefficientType;
//typedef Z2 CoefficientType;

int main(int argc, char **argv)
{
    using ComplexType = CGAL::HDVF::Cubical_chain_complex<CoefficientType> ;
    using HDVFType = CGAL::HDVF::Hdvf<CoefficientType, ComplexType> ;

    if (argc != 2)
    {
        std::cout << "usage: example_hdvf_cubical pgm_file" << std::endl;
    }
    else
    {
        // Choose between PRIMAL and DUAL construction
        const ComplexType::typeComplexCube primal_dual = ComplexType::PRIMAL;
        // Adapt pgm loading into Cub_complex accordingly
        const bool khalimsky_coords = (primal_dual == ComplexType::PRIMAL) ? true : false ;

        CGAL::HDVF::Cub_object_io mesh ;

        // Load pgm into cub object
        mesh.read_pgm(argv[1], khalimsky_coords);
        mesh.print_infos();

        // Build simplicial chain complex
        ComplexType complex(mesh, primal_dual);

        std::cout << complex;

        // Build empty HDVF
        HDVFType hdvf(complex, CGAL::HDVF::OPT_FULL) ;

        // Compute a perfect HDVF
        hdvf.compute_perfect_hdvf();
        //        hdvf.compute_rand_perfect_hdvf();

        // Output HDVF to console
        hdvf.insert_matrices();
        hdvf.insert_reduction();

        // Output HDVF to vtk
        hdvf_geometric_chain_complex_output_vtk(hdvf, complex, "res", true) ;
    }

    return 0;
}
