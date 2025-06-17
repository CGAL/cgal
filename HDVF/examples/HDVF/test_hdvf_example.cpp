#include <iostream>
#include <fstream>
#include <chrono>
#include "CGAL/Hdvf/tools_io.hpp"
#include "CGAL/Hdvf/Cubical_chain_complex.hpp"
#include "CGAL/HDVF/Geometric_chain_complex_tools.h"
#include "CGAL/Hdvf/Hdvf.h"

#include "CGAL/OSM/OSM.hpp"

typedef int CoefficientType;

using namespace CGAL;
using namespace HDVF;

int main(int argc, char **argv)
{
    using ComplexType = Cubical_chain_complex<CoefficientType> ;
    using HDVFType = Hdvf<CoefficientType, ComplexType> ;
    
    if (argc != 2)
    {
        std::cout << "usage: test_CubObject cub_file" << std::endl;
    }
    else
    {
        // Load cub object
        Cub_object mesh ;
        mesh.read_cub(argv[1], true);
            
        mesh.print_infos();
        
        // Build complex (PRIMAL construction)
        ComplexType complex(mesh, ComplexType::PRIMAL);
        
        complex.print_complex();
        
        // Build empty HDVF
        HDVFType hdvf(complex, OPT_FULL) ;
        
        // Build HDVF step by step
        hdvf.A(7,0,1);
        hdvf.A(2,5,0);
        hdvf.A(3,8,0);
        hdvf.compute_perfect_hdvf();
        
        // Output HDVF to vtk
        hdvf_geometric_chain_complex_output_vtk<CoefficientType, ComplexType>(hdvf, complex, "res") ;
        
    }
    return 0;
}
