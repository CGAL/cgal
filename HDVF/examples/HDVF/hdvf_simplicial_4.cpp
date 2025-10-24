#include <iostream>
#include <fstream>
#include <CGAL/Epick_d.h>
#include <CGAL/HDVF/Hdvf_traits_d.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

//typedef int Coefficient_ring;
//typedef HDVF::Z2 Coefficient_ring;
typedef HDVF::Zp<5, char, true> Coefficient_ring;

typedef CGAL::Epick_d<CGAL::Dimension_tag<4>> Kernel;
typedef HDVF::Hdvf_traits_d<Kernel> Traits;

int main(int argc, char **argv)
{
#if 1
    using Complex = HDVF::Simplicial_chain_complex<Coefficient_ring,Traits> ;
    using HDVF_type = HDVF::Hdvf<Complex> ;

    std::string nodes_filename, simp_filename ;
    if (argc > 3) std::cerr << "usage: hdvf_simplicial_4 nodes_file simp_file" << std::endl;
    else if (argc == 1) {
        nodes_filename  = "data/dim_4/klein4.nodes";
        simp_filename = "data/dim_4/klein4.simp";
    }
    else {
        nodes_filename = argv[1];
        simp_filename = argv[2];
    }

    // Load nodes + simplicial data
    HDVF::Mesh_object_io<Traits> mesh ;
    mesh.read_nodes_file(nodes_filename);
    mesh.read_simp(simp_filename);

    mesh.print_infos();

    // Build simplicial chain complex
    Complex complex(mesh);

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
/   Traits::to_point3 = Traits::pca_frame_builder(complex.points());
    CGAL::IO::write_VTK(hdvf, complex, "tmp/res", true) ;

    // Save HDVF to .hdvf file
    hdvf.write_hdvf_reduction("tmp/test.hdvf") ;

#endif
    return 0;
}
