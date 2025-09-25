#include <iostream>
#include <fstream>
#include <functional>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Filtration_lower_star.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Hdvf_duality.h>
#include <CGAL/HDVF/Sub_chain_complex_mask.h>
#include <CGAL/OSM/OSM.h>

//typedef CGAL::Homological_discrete_vector_field::Zp<5,int,true> CoefficientType;
typedef CGAL::Homological_discrete_vector_field::Z2 Coefficient_ring;
typedef CGAL::Homological_discrete_vector_field::Simplicial_chain_complex<Coefficient_ring> Complex;
typedef CGAL::Homological_discrete_vector_field::Hdvf_duality<Complex> HDVF_type;
typedef CGAL::Homological_discrete_vector_field::Duality_simplicial_complex_tools<Coefficient_ring> Tools_type;
typedef CGAL::Homological_discrete_vector_field::Sub_chain_complex_mask<Complex> Sub_chain_complex;

int main(int argc, char **argv)
{


    if (argc != 2)
    {
        std::cout << "usage: example_hdvf_simplicial off_file" << std::endl;
    }
    else
    {
        // Load cub object
        CGAL::Homological_discrete_vector_field::Mesh_object_io mesh ;
        mesh.read_off(argv[1]);

        mesh.print_infos();

        // Build simplicial chain complex
        Complex* complex = new Complex(mesh);

        std::cout << complex;

        // Build K and L
        typename Tools_type::Complex_duality_data t(Tools_type::simplicial_chain_complex_bb(*complex)) ;
        delete complex ;
        Complex& L(t.L) ;
        Sub_chain_complex& K(t.K) ;

        std::cout << "----> complex informations" << std::endl ;
        std::cout << "------> complex L" << std::endl ;
        std::cout << L;
        std::cout << "------> subcomplex K" << std::endl ;
        std::cout << K << std::endl ;

        // Create and compute a perfect HDVF
        HDVF_type hdvf(L, K, CGAL::Homological_discrete_vector_field::OPT_FULL);
        hdvf.compute_perfect_hdvf();

        // Export K HDVF
        hdvf.set_mask_K();
        CGAL::IO::write_VTK(hdvf, L, "res_complex_K", false) ;
        // Export L-K HDVF
        hdvf.set_mask_L_K();
        CGAL::IO::write_VTK(hdvf, L, "res_cocomplex_L_K", false) ;
        // Compute pairing
        std::vector<CGAL::Homological_discrete_vector_field::Cell_pair> pairs = hdvf.compute_pairing_hdvf();
        // Output pairing
        for (const auto& pair : pairs) {
            std::cout << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
        }

    }

    return 0;
}
