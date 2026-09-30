#include <iostream>
#include <fstream>
#include <functional>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Zp.h>
#include <CGAL/Z2.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Filtration_lower_star.h>
#include <CGAL/HDVF/Hdvf_duality.h>
#include <CGAL/HDVF/Sub_chain_complex_mask.h>
#include <CGAL/OSM/OSM.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;

typedef Kernel::Point_3 Point_3;

typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

//typedef CGAL::Zp<5,int,true> Coefficient_ring;
typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::OSM::Sub_sparse_matrix<CGAL::OSM::Sparse_chain> Sparse_matrix_struct;
typedef HDVF::Simplicial_chain_complex<Coefficient_ring, Traits, Sparse_matrix_struct> Complex;
typedef HDVF::Hdvf_duality<Complex> HDVF_type;
typedef HDVF::Duality_simplicial_complex_tools<Coefficient_ring,Traits> Tools_type;
typedef HDVF::Sub_chain_complex_mask<Complex> Sub_chain_complex;

int main(int argc, char **argv)
{
    std::string filename ;
    if (argc > 2) std::cout << "usage: dual_hdvf_simplicial mesh_file" << std::endl;
    else if (argc == 1) filename  = "data/data_simplicial/three_triangles.obj";
    else filename = argv[1];

    // Load mesh object
    Surface_mesh sm;
    CGAL::IO::read_polygon_mesh(filename, sm);
    HDVF::Surface_mesh_io<Surface_mesh,Traits> mesh(sm) ;

    mesh.print_infos();

    // Build K and L
    typename Tools_type::Complex_duality_data t = Tools_type::dualize_complex(mesh, 1.5, "tmp/file_K_closed.off", 1) ;
    std::cout << "--- Triangulation built" << std::endl ;
    std::shared_ptr<Complex> L(t.L_complex) ;
    std::shared_ptr<Sub_chain_complex> K(t.K_complex) ;
    std::cout << "--- K,L built" << std::endl ;

    std::cout << "----> complex informations" << std::endl ;
    std::cout << "------> complex L" << std::endl ;
    std::cout << *L;
    std::cout << "------> subcomplex K" << std::endl ;
    std::cout << *K << std::endl ;

    // Create and compute a perfect HDVF over K, L respectively
    HDVF_type hdvf(*L, *K, HDVF::OPT_FULL);
    hdvf.compute_perfect_hdvf();

    // Export K HDVF
    std::cout << "======> K" << std::endl;
    hdvf.set_mask_K();
    hdvf.print_reduction_sub();
    CGAL::IO::write_VTK(hdvf, *L, "tmp/res_complex_K", false) ;
    // Export L-K HDVF
    std::cout << "======> L-K" << std::endl;
    hdvf.set_mask_L_K();
    hdvf.print_reduction_sub();
    CGAL::IO::write_VTK(hdvf, *L, "tmp/res_cocomplex_L_K", false) ;
    // Compute pairing
    std::vector<HDVF::Cell_pair> pairs = hdvf.compute_pairing_hdvf();
    // Output pairing
    std::cout << "---> Alexander duality pairs" << std::endl;
    for (const auto& pair : pairs) {
        std::cout << "Sigma: " << pair.sigma << ", Tau: " << pair.tau << ", Dim: " << pair.dim << std::endl;
    }

    return 0;
}
