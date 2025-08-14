#include <iostream>
#include <fstream>
#include <functional>
#include <CGAL/HDVF/Mesh_object_io.h>
#include <CGAL/HDVF/Simplicial_chain_complex.h>
#include <CGAL/HDVF/Geometric_chain_complex_tools.h>
#include <CGAL/HDVF/Filtration_lower_star.h>
#include <CGAL/HDVF/Zp.h>
#include <CGAL/HDVF/Z2.h>
#include <CGAL/HDVF/Hdvf_persistence.h>
#include <CGAL/OSM/OSM.h>

//typedef int CoefficientType;
//typedef CGAL::HDVF::Zp<5> CoefficientType;
typedef CGAL::HDVF::Z2 CoefficientType;
typedef CGAL::HDVF::Simplicial_chain_complex<CoefficientType> ComplexType;
typedef double DegreeType;
typedef CGAL::HDVF::Filtration_lower_star<CoefficientType,ComplexType, DegreeType> FiltrationType;
typedef CGAL::HDVF::Hdvf_persistence<CoefficientType, ComplexType, DegreeType, FiltrationType> HDVFType;


int main(int argc, char **argv)
{
   

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

        std::cout << complex;

        /* Example 1 : build a lower star filtration "slicing" the object according to the "z" coordinate of vertices
            -> GeometricChainComplex required in this case
         */
        {
            // --- First: build the function mapping the index of a vertex to its degree (here the "z" coordinate of a vertex
            std::function<DegreeType(size_t)> f(CGAL::HDVF::deg_fun(complex, CGAL::HDVF::f_z));
            
            // -- Second: build the associated lower star filtration
            FiltrationType filtration(complex, f);
            
            
            // Build empty persistant HDVF (with vtk export activated)
            HDVFType hdvf(complex, filtration, CGAL::HDVF::OPT_FULL, true);
            
            // Compute a perfect HDVF
            hdvf.compute_perfect_hdvf();
            
            // Output HDVF to console
            hdvf.insert_matrices();
            hdvf.insert_reduction();
            std::cout << hdvf ;
            
            // Output HDVF to vtk
            hdvf_persistence_geometric_chain_complex_output_vtk(hdvf, complex, "res", true) ;
        }
        
        /* Example 2 : build a lower star filtration "slicing" the object according to the index of vertices (we can imagine any other filtration that does not depend on the geometry (e.g. color of vertices...)
            -> AbstractChainComplex in this case
         */
        {
            // --- First: build the function mapping the index of a vertex to its degree (here the index itself)
            std::function<DegreeType(size_t)> f = [&complex](size_t i)
            {
                return DegreeType(i) ;
            } ;
            
            // -- Second: build the associated lower star filtration
            FiltrationType filtration(complex, f);
            
            
            // Build empty persistant HDVF (with vtk export activated)
            HDVFType hdvf(complex, filtration, CGAL::HDVF::OPT_FULL, true);
            
            // Compute a perfect HDVF
            hdvf.compute_perfect_hdvf();
            
            // Output HDVF to console
            hdvf.insert_matrices();
            hdvf.insert_reduction();
            std::cout << hdvf ;
            
            // Output HDVF to vtk
            hdvf_persistence_geometric_chain_complex_output_vtk(hdvf, complex, "res2", true) ;
        }
    }

    return 0;
}
