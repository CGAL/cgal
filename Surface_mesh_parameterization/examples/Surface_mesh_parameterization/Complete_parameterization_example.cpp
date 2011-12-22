#include <CGAL/basic.h> // include basic.h before testing #defines

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Square_border_parameterizer_3.h>
#include <CGAL/Parameterization_mesh_patch_3.h>

#include <CGAL/Eigen_solver_traits.h>

#include <iostream>
#include <fstream>
#include <cstdlib>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

typedef CGAL::Cartesian<double>             Kernel;
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

// Dump parameterized mesh to an eps file
static bool write_file_eps(const Parameterization_polyhedron_adaptor& mesh_adaptor,
                           const char *pFilename,
                           double scale = 500.0)
{
    const Polyhedron& mesh = mesh_adaptor.get_adapted_mesh();

    std::ofstream out(pFilename);
    if(!out)
        return false;
    CGAL::set_ascii_mode(out);

    // compute bounding box
    double xmin,xmax,ymin,ymax;
    xmin = ymin = xmax = ymax = 0;
    Polyhedron::Halfedge_const_iterator pHalfedge;
    for (pHalfedge = mesh.halfedges_begin();
         pHalfedge != mesh.halfedges_end();
         pHalfedge++)
    {
        double x1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().x();
        double y1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().y();
        double x2 = scale * mesh_adaptor.info(pHalfedge)->uv().x();
        double y2 = scale * mesh_adaptor.info(pHalfedge)->uv().y();
        xmin = (std::min)(xmin,x1);
        xmin = (std::min)(xmin,x2);
        xmax = (std::max)(xmax,x1);
        xmax = (std::max)(xmax,x2);
        ymax = (std::max)(ymax,y1);
        ymax = (std::max)(ymax,y2);
        ymin = (std::min)(ymin,y1);
        ymin = (std::min)(ymin,y2);
    }

    out << "%!PS-Adobe-2.0 EPSF-2.0" << std::endl;
    out << "%%BoundingBox: " << int(xmin+0.5) << " "
                                << int(ymin+0.5) << " "
                                << int(xmax+0.5) << " "
                                << int(ymax+0.5) << std::endl;
    out << "%%HiResBoundingBox: " << xmin << " "
                                    << ymin << " "
                                    << xmax << " "
                                    << ymax << std::endl;
    out << "%%EndComments" << std::endl;
    out << "gsave" << std::endl;
    out << "0.1 setlinewidth" << std::endl;

    // color macros
    out << std::endl;
    out << "% RGB color command - r g b C" << std::endl;
    out << "/C { setrgbcolor } bind def" << std::endl;
    out << "/white { 1 1 1 C } bind def" << std::endl;
    out << "/black { 0 0 0 C } bind def" << std::endl;

    // edge macro -> E
    out << std::endl;
    out << "% Black stroke - x1 y1 x2 y2 E" << std::endl;
    out << "/E {moveto lineto stroke} bind def" << std::endl;
    out << "black" << std::endl << std::endl;

    // for each halfedge
    for (pHalfedge = mesh.halfedges_begin();
         pHalfedge != mesh.halfedges_end();
         pHalfedge++)
    {
        double x1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().x();
        double y1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().y();
        double x2 = scale * mesh_adaptor.info(pHalfedge)->uv().x();
        double y2 = scale * mesh_adaptor.info(pHalfedge)->uv().y();
        out << x1 << " " << y1 << " " << x2 << " " << y2 << " E" << std::endl;
    }

    /* Emit EPS trailer. */
    out << "grestore" << std::endl;
    out << std::endl;
    out << "showpage" << std::endl;

    return true;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "PARAMETERIZATION" << std::endl;
    std::cerr << "  Discrete Authalic Parameterization" << std::endl;
    std::cerr << "  Square border" << std::endl;
    std::cerr << "  Eigen solver" << std::endl;
    std::cerr << "  Very simple cut if model is not a topological disk" << std::endl;
    std::cerr << "  Output: EPS" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 != 2)
    {
        std::cerr << "Usage: " << argv[0] << " input_file.off output_file.eps" << std::endl;
        return(EXIT_FAILURE);
    }

    // File names are:
    const char* input_filename  = argv[1];
    const char* output_filename = argv[2];

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
    // Discrete Authalic Parameterization (square border)
    // with Eigen solver
    //***************************************

    // Border parameterizer
    typedef CGAL::Square_border_arc_length_parameterizer_3<Mesh_patch_polyhedron>
                                                            Border_parameterizer;
    // Eigen solver
    typedef CGAL::Eigen_solver_traits<>                Solver;

    // Discrete Authalic Parameterization (square border)
    // with Eigen solver
    typedef CGAL::Discrete_authalic_parameterizer_3<Mesh_patch_polyhedron,
                                                    Border_parameterizer,
                                                    Solver> Parameterizer;

    Parameterizer::Error_code err = CGAL::parameterize(mesh_patch, Parameterizer());
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

    // Write Postscript file
    if ( ! write_file_eps(mesh_adaptor, output_filename) )
    {
        std::cerr << "Error: cannot write file " << output_filename << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
