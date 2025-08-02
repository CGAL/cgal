#include <iostream>
#include <fstream>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Hdvf.h>
#include <CGAL/OSM/OSM.h>

//typedef int CoefficientType;
typedef CGAL::HDVF::Z2 CoefficientType;
//typedef CGAL::HDVF::Zp<5> CoefficientType;

int main(int argc, char **argv)
{
    using ComplexType = CGAL::HDVF::Simplicial_chain_complex<CoefficientType> ;
    using HDVFType = CGAL::HDVF::Hdvf<CoefficientType, ComplexType> ;
    
    if (argc != 2)
    {
        std::cout << "usage: example_hdvf_simplicial off_file" << std::endl;
    }
    else
    {
        // Load cub object
        CGAL::HDVF::Mesh_object_io mesh ;
        mesh.read_off(argv[1]);
        
        mesh.print_infos();
        
        // Build simplicial chain complex
        ComplexType complex(mesh, mesh.get_nodes());
        
        complex.print_complex();
        
        // Build empty HDVF
        HDVFType hdvf(complex, CGAL::HDVF::OPT_FULL) ;
        
        // Compute a perfect HDVF
        hdvf.compute_perfect_hdvf();
        //        hdvf.compute_rand_perfect_hdvf();
        
        // Output HDVF to console
        hdvf.print_matrices();
        hdvf.print_reduction();
        
        // Output HDVF to vtk
        hdvf_geometric_chain_complex_output_vtk(hdvf, complex, "res", true) ;
    }
    
    return 0;
}
