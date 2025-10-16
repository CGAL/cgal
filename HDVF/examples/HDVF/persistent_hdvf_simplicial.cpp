#include <iostream>
#include <fstream>
#include <functional>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/HDVF/Hdvf_traits_3.h>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Filtration_lower_star.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Hdvf_persistence.h>
#include <CGAL/OSM/OSM.h>

namespace HDVF = CGAL::Homological_discrete_vector_field;

//typedef int Coefficient_ring;
//typedef HDVF::Z2 Coefficient_ring;
typedef HDVF::Zp<5, char, true> Coefficient_ring;
typedef CGAL::Simple_cartesian<double> Kernel;
typedef HDVF::Hdvf_traits_3<Kernel> Traits;
typedef Kernel::Compute_z_3 Compute_z;
typedef Kernel::Point_3 Point_3;

typedef HDVF::Simplicial_chain_complex<Coefficient_ring,Traits> Complex;
typedef double Degree;
typedef HDVF::Filtration_lower_star<Complex, Degree> FiltrationType;
typedef HDVF::Hdvf_persistence<Complex, Degree, FiltrationType> HDVF_type;


int main(int argc, char **argv)
{
    
    std::string filename;
    if (argc > 2) std::cout << "usage: persistent_hdvf_simplicial off_file" << std::endl;
    else if (argc == 1) filename  = "data/mesh_data/two_rings.off";
    else filename = argv[1];
    
    // Load cub object
    HDVF::Mesh_object_io<Traits> mesh ;
    mesh.read_off(filename);
    
    mesh.print_infos();
    
    // Build simplicial chain complex
    Complex complex(mesh);
    
    std::cout << complex;
    
    /* Example 1 : build a lower star filtration "slicing" the object according to the "z" coordinate of vertices
     -> GeometricChainComplex required in this case
     */
    {
        // --- First: build the function mapping the index of a vertex to its degree (here the "z" coordinate of a vertex
        Compute_z compute_z = Kernel().compute_z_3_object() ;
        std::function<Degree(size_t)> f(HDVF::degree_function<Complex,Point_3>(complex, compute_z));
        
        // -- Second: build the associated lower star filtration
        FiltrationType filtration(complex, f);
        
        
        // Build empty persistent HDVF (with vtk export activated)
        HDVF_type hdvf(complex, filtration, HDVF::OPT_FULL, true);
        
        // Compute a perfect HDVF
        hdvf.compute_perfect_hdvf();
        
        // Output HDVF to console
        hdvf.write_matrices();
        hdvf.write_reduction();
        std::cout << hdvf ;
        
        // Output HDVF to vtk
        CGAL::IO::write_VTK(hdvf, complex, "tmp/res", true) ;
    }
    
    /* Example 2 : build a lower star filtration "slicing" the object according to the index of vertices (we can imagine any other filtration that does not depend on the geometry (e.g. color of vertices...)
     -> AbstractChainComplex in this case
     */
    {
        // --- First: build the function mapping the index of a vertex to its degree (here the index itself)
        std::function<Degree(size_t)> f = [&complex](size_t i)
        {
            return Degree(i) ;
        } ;
        
        // -- Second: build the associated lower star filtration
        FiltrationType filtration(complex, f);
        
        
        // Build empty persistent HDVF (with vtk export activated)
        HDVF_type hdvf(complex, filtration, HDVF::OPT_FULL, true);
        
        // Compute a perfect HDVF
        hdvf.compute_perfect_hdvf();
        
        // Output HDVF to console
        hdvf.write_matrices();
        hdvf.write_reduction();
        std::cout << hdvf ;
        
        // Output HDVF to vtk
        CGAL::IO::write_VTK(hdvf, complex, "tmp/res2", true) ;
    }
    
    return 0;
}
