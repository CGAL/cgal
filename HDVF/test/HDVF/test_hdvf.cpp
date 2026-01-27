#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Zp.h>
#include <CGAL/Z2.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

//#define BUILD_TEST_DATA

namespace HDVF = CGAL::Homological_discrete_vector_field;

//typedef int Coefficient_ring;
//typedef CGAL::Z2 Coefficient_ring;
typedef CGAL::Zp<5, char, true> Coefficient_ring;

typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;

using Complex = HDVF::Simplicial_chain_complex<Coefficient_ring,Traits> ;
using HDVF_type = HDVF::Hdvf<Complex> ;

int main(int argc, char **argv) {
    std::string filename;
    if (argc > 2) {
        std::cerr << "usage: test_hdvf_core [off_file]" << std::endl;
    }
    else if (argc == 1) filename = "data/three_triangles.off" ;
    else filename = argv[1] ;

    // Load off into Mesh_object_io
    HDVF::Mesh_object_io<Traits> mesh ;
    mesh.read_off(filename);

//    mesh.print_infos();

    // Build simplicial chain complex
    Complex complex(mesh);

    std::cout << complex;

    // Build empty HDVF
    HDVF_type hdvf(complex, HDVF::OPT_FULL) ;

    // Compute a perfect HDVF
    hdvf.compute_perfect_hdvf();
    //        hdvf.compute_rand_perfect_hdvf();
    std::cerr << std::endl;

#ifdef BUILD_TEST_DATA
    // Save HDVF to .hdvf file
    hdvf.write_hdvf_reduction("data/test_hdvf/test_hdvf.hdvf") ;
#endif

    // Read HDVF from .hdvf
    HDVF_type hdvf2(complex, HDVF::OPT_FULL);
    hdvf2.read_hdvf_reduction("data/test_hdvf/test_hdvf.hdvf");

    // Compare
    bool test_true(hdvf.compare(hdvf2));
    std::cerr << "-- Test HDVF built: " << test_true << std::endl;
    assert(test_true);

//    // Test z
//    CGAL::IO::write_VTK(hdvf, complex, "tmp/test_hdvf");
//
//    // z over a critical cell
//    std::vector<size_t> criticals(hdvf.psc_flags(HDVF::CRITICAL,1));
//    HDVF_type::Column_chain test_z(hdvf.z(criticals.at(0),1));
//    CGAL::IO::write_VTK(complex, "tmp/z_critical.vtk", test_z, 1);
//
//    // z over a primary cell
//    std::vector<size_t> primary(hdvf.psc_flags(HDVF::PRIMARY,1));
//    HDVF_type::Column_chain test_z_p(hdvf.z(primary.at(0),1));
//    Complex::chain_to_vtk(complex, "tmp/z_primary.vtk", test_z_p, 1);
//
//    // Test co_z
//
//    // co_z over a critical cell
//    HDVF_type::Column_chain test_co_z(hdvf.co_z(criticals.at(0),1));
//    CGAL::IO::write_VTK(complex, "tmp/co_z_critical.vtk", test_co_z, 1);
//
//    // co_z over a secondary cell
//    std::vector<size_t> secondary(hdvf.psc_flags(HDVF::SECONDARY,1));
//    HDVF_type::Column_chain test_co_z_s(hdvf.co_z(secondary.at(0),1));
//    CGAL::IO::write_VTK(complex, "tmp/co_z_secondary.vtk", test_co_z_s, 1);

    return 0;
}


