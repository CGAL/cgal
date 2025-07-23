#include <iostream>
#include <fstream>
#include "CGAL/HDVF/tools_io.hpp"
#include "CGAL/HDVF/Abstract_simplicial_chain_complex.hpp"
#include "CGAL/HDVF/Geometric_chain_complex_tools.h"
#include "CGAL/HDVF/Zp.hpp"
#include "CGAL/HDVF/Hdvf.h"

#include "CGAL/OSM/OSM.hpp"

using namespace CGAL;
using namespace HDVF;

//typedef int CoefficientType;
typedef Zp<2> CoefficientType;

int main(int argc, char **argv)
{
    using ComplexType = Simplicial_chain_complex<CoefficientType> ;
    using HDVFType = Hdvf<CoefficientType, ComplexType> ;
    
    if (argc != 2)
    {
        std::cout << "usage: example_hdvf_simplicial off_file" << std::endl;
    }
    else
    {
        // Load cub object
        Mesh_object mesh ;
        mesh.read_off(argv[1]);
        
        mesh.print_infos();
        
        // Build simplicial chain complex
        ComplexType complex(mesh, mesh.get_nodes());
        
        complex.print_complex();
        
        // Build empty HDVF
        HDVFType hdvf(complex, OPT_FULL) ;
        
        // Compute a perfect HDVF
        hdvf.compute_perfect_hdvf();
        //        hdvf.compute_rand_perfect_hdvf();
        
        // Output HDVF to console
        hdvf.print_matrices();
        hdvf.print_reduction();
        
        // Output HDVF to vtk
        hdvf_geometric_chain_complex_output_vtk<CoefficientType, ComplexType>(hdvf, complex, "res", true) ;
    }
    
    return 0;
}
