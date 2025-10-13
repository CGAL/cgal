#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Cub_object_io.h>
#include <CGAL/HDVF/Cubical_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

//typedef int Coefficient_ring;
//typedef HDVF::Z2 Coefficient_ring;
typedef HDVF::Zp<5, char, true> Coefficient_ring;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;

int main(int argc, char **argv)
{
    using Complex = CGAL::Homological_discrete_vector_field::Cubical_chain_complex<Coefficient_ring, Traits> ;
    using HDVF_type = CGAL::Homological_discrete_vector_field::Hdvf<Complex> ;

    if (argc != 2)
    {
        std::cout << "usage: example_hdvf_cubical pgm_file" << std::endl;
    }
    else
    {
        // Choose between PRIMAL and DUAL construction
        const Complex::Cubical_complex_primal_dual primal_dual = Complex::PRIMAL;
        // Adapt pgm loading into Cub_complex accordingly
        const bool khalimsky_coords = (primal_dual == Complex::PRIMAL) ? true : false ;

        CGAL::Homological_discrete_vector_field::Cub_object_io mesh ;

        // Load pgm into cub object
        mesh.read_pgm(argv[1], khalimsky_coords);
        mesh.print_infos();

        // Build simplicial chain complex
        Complex complex(mesh, primal_dual);

        std::cout << complex;

        // Build empty HDVF
        HDVF_type hdvf(complex, CGAL::Homological_discrete_vector_field::OPT_FULL) ;

        // Compute a perfect HDVF
        hdvf.compute_perfect_hdvf();
        //        hdvf.compute_rand_perfect_hdvf();

        // Output HDVF to console
        hdvf.insert_matrices();
        hdvf.insert_reduction();

        // Output HDVF to vtk
        CGAL::IO::write_VTK(hdvf, complex, "res", true) ;
    }

    return 0;
}
