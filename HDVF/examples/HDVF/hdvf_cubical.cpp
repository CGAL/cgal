#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Hdvf_traits_2.h>
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
//typedef HDVF::Hdvf_traits_3<Kernel> Traits;
typedef HDVF::Hdvf_traits_2<Kernel> Traits;

int main(int argc, char **argv)
{
    using Complex = HDVF::Cubical_chain_complex<Coefficient_ring, Traits> ;
    using HDVF_type = HDVF::Hdvf<Complex> ;
    
    std::string filename ;
    if (argc > 2) std::cout << "usage: example_hdvf_cubical pgm_file" << std::endl;
    else if (argc == 1) filename  = "data/cub_data/Eight_3D.pgm";
    else filename = argv[1];
    
    // Choose between PRIMAL and DUAL construction
    const Complex::Cubical_complex_primal_dual primal_dual = Complex::PRIMAL;
    // Adapt pgm loading into Cub_complex accordingly
    const bool khalimsky_coords = (primal_dual == Complex::PRIMAL) ? true : false ;
    
    HDVF::Cub_object_io<Traits> mesh ;
    
    // Load pgm into cub object
    mesh.read_pgm(filename, khalimsky_coords);
    mesh.print_infos();
    
    // Build simplicial chain complex
    Complex complex(mesh, primal_dual);
    
    std::cout << complex;
    
    // Build empty HDVF
    HDVF_type hdvf(complex, HDVF::OPT_FULL) ;
    
    // Compute a perfect HDVF
    hdvf.compute_perfect_hdvf();
    //        hdvf.compute_rand_perfect_hdvf();
    
    // Output HDVF to console
    hdvf.write_matrices();
    hdvf.write_reduction();
    
    // Output HDVF to vtk
    CGAL::IO::write_VTK(hdvf, complex, "tmp/res", true) ;
    
    return 0;
}
