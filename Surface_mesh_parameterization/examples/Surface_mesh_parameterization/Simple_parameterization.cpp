#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>

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
    std::cerr << "  Circle border" << std::endl;
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
        std::cerr << "FATAL ERROR: cannot read OFF file " << input_filename << std::endl;
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
    // (defaults are circular border and OpenNL solver)
    //***************************************

    // Type that defines the error codes
    typedef CGAL::Parameterizer_traits_3<Parameterization_polyhedron_adaptor>
                                            Parameterizer;

    Parameterizer::Error_code err = CGAL::parameterize(mesh_adaptor);
    if (err != Parameterizer::OK)
        std::cerr << "FATAL ERROR: " << Parameterizer::get_error_message(err) << std::endl;

    //***************************************
    // Output
    //***************************************

    if (err == Parameterizer::OK)
    {
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
    }

    return (err == Parameterizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}
