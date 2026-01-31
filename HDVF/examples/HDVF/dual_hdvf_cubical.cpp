#include <iostream>
#include <fstream>
#include <functional>
#include <CGAL/Zp.h>
#include <CGAL/Z2.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Cub_object_io.h>
#include <CGAL/HDVF/Cubical_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Filtration_lower_star.h>
#include <CGAL/HDVF/Hdvf_duality.h>
#include <CGAL/HDVF/Sub_chain_complex_mask.h>
#include <CGAL/OSM/OSM.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;

//typedef HDVF::Zp<5,int,true> Coefficient_ring;
typedef CGAL::Z2 Coefficient_ring;
typedef HDVF::Cubical_chain_complex<Coefficient_ring, Traits> Complex;
typedef HDVF::Hdvf_duality<Complex> HDVF_type;
typedef HDVF::Duality_cubical_complex_tools<Coefficient_ring,Traits> Tools_type;
typedef HDVF::Sub_chain_complex_mask<Complex> Sub_chain_complex;

// PRIMAL / DUAL construction

#define PRIMAL_DUAL true // or false

// Frame or not the complex (add 1 pixel around)

#define FRAME true // or false

int main(int argc, char **argv)
{
    std::string filename ;
    if (argc > 2) std::cout << "usage: dual_hdvf_cubical pgm_file" << std::endl;
    else if (argc == 1) filename  = "data/data_cubical/Eight_10.pgm";
    else filename = argv[1];

    // Set variables to chose loading mode (PRIMAL or DUAL associated complex, pgm read using Khalimsky coordinates or not)
    Complex::Cubical_complex_primal_dual primal_dual;
    bool khalimsky;
    if (PRIMAL_DUAL) {
        primal_dual = Complex::PRIMAL;
        khalimsky = true;
    }
    else {
        primal_dual = Complex::DUAL;
        khalimsky = false;
    }

    // Load pgm object
    HDVF::Cub_object_io<Traits> cub_object ;
    cub_object.read_pgm(filename, khalimsky);

    // Frame the object
    if (FRAME)
        cub_object.frame();

    cub_object.print_infos();

    // Build K and L
    typename Tools_type::Complex_duality_data t(Tools_type::dualize_complex(cub_object, primal_dual)) ;
    std::cout << "--- Triangulation built" << std::endl ;
    std::shared_ptr<Complex> L(t.L_complex) ;
    std::shared_ptr<Sub_chain_complex> K(t.K_complex) ;
    std::cout << "--- K,L built" << std::endl ;

    std::cout << "----> complex informations" << std::endl ;
    std::cout << "------> complex L" << std::endl ;
    std::cout << *L;
    std::cout << "------> subcomplex K" << std::endl ;
    std::cout << *K << std::endl ;

    // Create and compute a perfect HDVF
    HDVF_type hdvf(*L, *K, HDVF::OPT_FULL);
    hdvf.compute_perfect_hdvf();

    // Export K HDVF
    hdvf.set_mask_K();
    CGAL::IO::write_VTK(hdvf, *L, "tmp/res_complex_K", false) ;
    // Export L-K HDVF
    hdvf.set_mask_L_K();
    CGAL::IO::write_VTK(hdvf, *L, "tmp/res_cocomplex_L_K", false) ;
    // Compute pairing
    std::vector<HDVF::Cell_pair> pairs = hdvf.compute_pairing_hdvf();
    // Output pairing
    for (const auto& pair : pairs) {
        std::cout << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
    }

    return 0;
}
