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
    CGAL::IO::write_VTK(hdvf, complex, "tmp/res", true) ;
    
    // Save HDVF to .hdvf file
    hdvf.write_hdvf_reduction("tmp/test.hdvf") ;
    
    return 0;
}


