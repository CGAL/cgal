#include <iostream>
#include <fstream>
#include <functional>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Surface_mesh_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Filtration_lower_star.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Hdvf_duality.h>
#include <CGAL/HDVF/Sub_chain_complex_mask.h>
#include <CGAL/OSM/OSM.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;

//typedef HDVF::Zp<5,int,true> Coefficient_ring;
typedef HDVF::Z2 Coefficient_ring;
typedef HDVF::Simplicial_chain_complex<Coefficient_ring, Traits> Complex;
typedef HDVF::Hdvf_duality<Complex> HDVF_type;
typedef HDVF::Duality_simplicial_complex_tools<Coefficient_ring,Traits> Tools_type;
typedef HDVF::Sub_chain_complex_mask<Complex> Sub_chain_complex;
typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

int main(int argc, char **argv)
{
    std::string filename ;
    if (argc > 2) std::cout << "usage: dual_hdvf_simplicial off_file" << std::endl;
    else if (argc == 1) filename  = "data/data_simplicial/two_rings.off";
    else filename = argv[1];
    
    // Load simplicial object into Surface_mesh
    Surface_mesh sm;
    CGAL::IO::read_polygon_mesh(filename, sm);
    
    if (sm.number_of_vertices() == 0)
    {
        std::cerr << "Emtpy mesh read" << std::endl;
        throw("Emtpy mesh read");
    }
    
//    std::cout << sm;
    
    // Build K and L
    typename Tools_type::Complex_duality_data t(Tools_type::simplicial_chain_complex_bb(sm, 1.5, 2));
    std::cout << "--- Triangulation built" << std::endl ;
    Complex& L(t.L);
    Sub_chain_complex& K(t.K);
    std::cout << "--- K,L built" << std::endl ;
    
    std::cout << "----> complex informations" << std::endl ;
    std::cout << "------> complex L" << std::endl ;
    std::cout << L;
    std::cout << "------> subcomplex K" << std::endl ;
    std::cout << K << std::endl ;
    
    // Create and compute a perfect HDVF
    HDVF_type hdvf(L, K, HDVF::OPT_FULL);
    hdvf.compute_perfect_hdvf();
    
    // Export K HDVF
    hdvf.set_mask_K();
    CGAL::IO::write_VTK(hdvf, L, "tmp/res_complex_K", false) ;
    // Export L-K HDVF
    hdvf.set_mask_L_K();
    CGAL::IO::write_VTK(hdvf, L, "tmp/res_cocomplex_L_K", false) ;
    // Compute pairing
    std::vector<HDVF::Cell_pair> pairs = hdvf.compute_pairing_hdvf();
    // Output pairing
    for (const auto& pair : pairs) {
        std::cout << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
    }

    return 0;
}
