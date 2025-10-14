#include <iostream>
#include <fstream>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/HDVF/Surface_mesh_io.h>
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

typedef CGAL::Surface_mesh<Kernel::Point_3> Surface_mesh;

int main(int argc, char **argv)
{
#if 1
    using Complex = HDVF::Simplicial_chain_complex<Coefficient_ring,Traits> ;
    using HDVF_type = HDVF::Hdvf<Complex> ;

    if (argc != 2)
    {
        std::cerr << "usage: example_hdvf_simplicial off_file" << std::endl;
    }
    else
    {
        // Load cub object
        Surface_mesh sm;
        HDVF::Surface_mesh_io<Surface_mesh,Traits> mesh(sm) ;

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
        CGAL::IO::write_VTK(hdvf, complex, "res", true) ;

        // Save HDVF to .hdvf file
        hdvf.write_hdvf_reduction("test.hdvf") ;
    }
#endif
    return 0;
}
