#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include <ostream>
#include <cassert>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

//#define BUILD_TEST_DATA

namespace HDVF = CGAL::Homological_discrete_vector_field;

//typedef int Coefficient_ring;
//typedef HDVF::Z2 Coefficient_ring;
typedef HDVF::Zp<5, char, true> Coefficient_ring;

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
    
    return 0;
}


