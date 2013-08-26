#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Parameterization_mesh_patch_3.h>

#include <iostream>
#include <fstream>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

typedef CGAL::Simple_cartesian<double>      Kernel;
typedef CGAL::Polyhedron_3<Kernel>          Polyhedron;

// Polyhedron adaptor
typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron>
                                            Parameterization_polyhedron_adaptor;

// Type describing a border or seam as a vertex list
typedef std::list<Parameterization_polyhedron_adaptor::Vertex_handle>
                                            Seam;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// If the mesh is a topological disk, extract its longest border,
// else compute a very simple cut to make it homeomorphic to a disk.
// Return the border of this region (empty on error)
//
// CAUTION: this cutting algorithm is very naive. Write your own!
static Seam cut_mesh(Parameterization_polyhedron_adaptor& mesh_adaptor)
{
    // Helper class to compute genus or extract borders
    typedef CGAL::Parameterization_mesh_feature_extractor<Parameterization_polyhedron_adaptor>
                                            Mesh_feature_extractor;

    Seam seam;              // returned list

    // Get reference to Polyhedron_3 mesh
    Polyhedron& mesh = mesh_adaptor.get_adapted_mesh();

    // Extract mesh borders and compute genus
    Mesh_feature_extractor feature_extractor(mesh_adaptor);
    int nb_borders = feature_extractor.get_nb_borders();
    int genus = feature_extractor.get_genus();

    // If mesh is a topological disk
    if (genus == 0 && nb_borders > 0)
    {
        // Pick the longest border
        seam = feature_extractor.get_longest_border();
    }
    else // if mesh is *not* a topological disk, create a virtual cut
    {
        const int CUT_LENGTH = 6;

        // Build consecutive halfedges array
        Polyhedron::Halfedge_handle seam_halfedges[CUT_LENGTH];
        seam_halfedges[0] = mesh.halfedges_begin();
        if (seam_halfedges[0] == NULL)
            return seam;                    // return empty list
        int i;
        for (i=1; i<CUT_LENGTH; i++)
        {
            seam_halfedges[i] = seam_halfedges[i-1]->next()->opposite()->next();
            if (seam_halfedges[i] == NULL)
                return seam;                // return empty list
        }

        // Convert halfedges array to two-ways vertices list
        for (i=0; i<CUT_LENGTH; i++)
            seam.push_back(seam_halfedges[i]->vertex());
        for (i=CUT_LENGTH-1; i>=0; i--)
            seam.push_back(seam_halfedges[i]->opposite()->vertex());
    }

    return seam;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "PARAMETERIZATION" << std::endl;
    std::cerr << "  Floater parameterization" << std::endl;
    std::cerr << "  Circle border" << std::endl;
    std::cerr << "  OpenNL solver" << std::endl;
    std::cerr << "  Very simple cut if model is not a topological disk" << std::endl;

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
    //***************************************

    Parameterization_polyhedron_adaptor mesh_adaptor(mesh);

    //***************************************
    // Virtually cut mesh
    //***************************************

    // The parameterization methods support only meshes that
    // are topological disks => we need to compute a "cutting" of the mesh
    // that makes it homeomorphic to a disk
    Seam seam = cut_mesh(mesh_adaptor);
    if (seam.empty())
    {
        std::cerr << "Input mesh not supported: the example cutting algorithm is too simple to cut this shape" << std::endl;
        return EXIT_FAILURE;
    }

    // Create a second adaptor that virtually "cuts" the mesh following the 'seam' path
    typedef CGAL::Parameterization_mesh_patch_3<Parameterization_polyhedron_adaptor>
                                            Mesh_patch_polyhedron;
    Mesh_patch_polyhedron   mesh_patch(mesh_adaptor, seam.begin(), seam.end());
    if (!mesh_patch.is_valid())
    {
        std::cerr << "Input mesh not supported: non manifold shape or invalid cutting" << std::endl;
        return EXIT_FAILURE;
    }

    //***************************************
    // Floater Mean Value Coordinates parameterization
    //***************************************

    typedef CGAL::Parameterizer_traits_3<Mesh_patch_polyhedron>
                                            Parameterizer; // Type that defines the error codes

    Parameterizer::Error_code err = CGAL::parameterize(mesh_patch);
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
