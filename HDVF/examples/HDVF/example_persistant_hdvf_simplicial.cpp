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

namespace HDVF = CGAL::Homological_discrete_vector_field;

//typedef int Coefficient_ring;
//typedef HDVF::Z2 Coefficient_ring;
typedef HDVF::Zp<5, char, true> Coefficient_ring;

typedef CGAL::Homological_discrete_vector_field::Simplicial_chain_complex<Coefficient_ring> Complex;
typedef double Degree;
typedef CGAL::Homological_discrete_vector_field::Filtration_lower_star<Complex, Degree> FiltrationType;
typedef CGAL::Homological_discrete_vector_field::Hdvf_persistence<Complex, Degree, FiltrationType> HDVF_type;


int main(int argc, char **argv)
{


    if (argc != 2)
    {
        std::cout << "usage: example_hdvf_simplicial off_file" << std::endl;
    }
    else
    {
        // Load cub object
        CGAL::Homological_discrete_vector_field::Mesh_object_io mesh ;
        mesh.read_off(argv[1]);

        mesh.print_infos();

        // Build simplicial chain complex
        Complex complex(mesh);

        std::cout << complex;

        /* Example 1 : build a lower star filtration "slicing" the object according to the "z" coordinate of vertices
            -> GeometricChainComplex required in this case
         */
        {
            // --- First: build the function mapping the index of a vertex to its degree (here the "z" coordinate of a vertex
            std::function<Degree(size_t)> f(CGAL::Homological_discrete_vector_field::degree_function(complex, CGAL::Homological_discrete_vector_field::f_z));

            // -- Second: build the associated lower star filtration
            FiltrationType filtration(complex, f);


            // Build empty persistant HDVF (with vtk export activated)
            HDVF_type hdvf(complex, filtration, CGAL::Homological_discrete_vector_field::OPT_FULL, true);

            // Compute a perfect HDVF
            hdvf.compute_perfect_hdvf();

            // Output HDVF to console
            hdvf.insert_matrices();
            hdvf.insert_reduction();
            std::cout << hdvf ;

            // Output HDVF to vtk
            CGAL::IO::write_VTK(hdvf, complex, "res", true) ;
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


            // Build empty persistant HDVF (with vtk export activated)
            HDVF_type hdvf(complex, filtration, CGAL::Homological_discrete_vector_field::OPT_FULL, true);

            // Compute a perfect HDVF
            hdvf.compute_perfect_hdvf();

            // Output HDVF to console
            hdvf.insert_matrices();
            hdvf.insert_reduction();
            std::cout << hdvf ;

            // Output HDVF to vtk
            CGAL::IO::write_VTK(hdvf, complex, "res2", true) ;
        }
    }

    return 0;
}
