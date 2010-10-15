#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Square_border_parameterizer_3.h>

#include <iostream>
#include <fstream>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

typedef CGAL::Cartesian<double>             Kernel;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "PARAMETERIZATION" << std::endl;
    std::cerr << "  Floater parameterization" << std::endl;
    std::cerr << "  square border" << std::endl;
    std::cerr << "  OpenNL solver" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 != 1)
    {
        std::cerr << "Usage: " << argv[0] << " input_file.off" << std::endl;
        return(EXIT_FAILURE);
    }

    // File name is:
    const char* input_filename  = argv[1];

    //***************************************
    // Read the mesh
    //***************************************

    // Read the mesh
    std::ifstream stream(input_filename);
    Polyhedron mesh;
    stream >> mesh;
    if(!stream || !mesh.is_valid() || mesh.empty())
    {
        std::cerr << "Error: cannot read OFF file " << input_filename << std::endl;
        return EXIT_FAILURE;
    }

    //***************************************
    // Create Polyhedron adaptor
    // Note: no cutting => we support only
    // meshes that are topological disks
    //***************************************

    typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron>
                                            Parameterization_polyhedron_adaptor;
    Parameterization_polyhedron_adaptor mesh_adaptor(mesh);

    //***************************************
    // Floater Mean Value Coordinates parameterization
    // with square border
    //***************************************

    // Square border parameterizer
    typedef CGAL::Square_border_arc_length_parameterizer_3<Parameterization_polyhedron_adaptor>
                                                       Border_parameterizer;

    // Floater Mean Value Coordinates parameterizer with square border
    typedef CGAL::Mean_value_coordinates_parameterizer_3<Parameterization_polyhedron_adaptor,
                                                         Border_parameterizer>
                                                        Parameterizer;

    Parameterizer::Error_code err = CGAL::parameterize(mesh_adaptor, Parameterizer());
    switch(err) {
    case Parameterizer::OK: // Success
        break;
    case Parameterizer::ERROR_EMPTY_MESH: // Input mesh not supported
    case Parameterizer::ERROR_NON_TRIANGULAR_MESH:   
    case Parameterizer::ERROR_NO_TOPOLOGICAL_DISC:     
    case Parameterizer::ERROR_BORDER_TOO_SHORT:    
        std::cerr << "Input mesh not supported: " << Parameterizer::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
        break;
    default: // Error
        std::cerr << "Error: " << Parameterizer::get_error_message(err) << std::endl;
        return EXIT_FAILURE;
        break;
    };

    //***************************************
    // Output
    //***************************************

    // Raw output: dump (u,v) pairs
    Polyhedron::Vertex_const_iterator pVertex;
    for (pVertex = mesh.vertices_begin();
        pVertex != mesh.vertices_end();
        pVertex++)
    {
        // (u,v) pair is stored in any halfedge
        double u = mesh_adaptor.info(pVertex->halfedge())->uv().x();
        double v = mesh_adaptor.info(pVertex->halfedge())->uv().y();
        std::cout << "(u,v) = (" << u << "," << v << ")" << std::endl;
    }

    return EXIT_SUCCESS;
}
