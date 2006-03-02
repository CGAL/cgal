// Complete_parameterization_example.C

#ifdef CGAL_USE_TAUCS


#include "short_names.h"                    // must be included first

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>
#include <CGAL/Square_border_parameterizer_3.h>
#include <CGAL/Parameterization_mesh_patch_3.h>

#include <CGAL/Taucs_solver_traits.h>

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

typedef CGAL::Cartesian<double>                         Kernel;
typedef CGAL::Polyhedron_3<Kernel>                      Polyhedron;

// Mesh adaptors
typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron>
                                                        Parameterization_polyhedron_adaptor;
typedef CGAL::Parameterization_mesh_patch_3<Parameterization_polyhedron_adaptor>
                                                        Mesh_patch_polyhedron;

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
// CAUTION:
// This method is provided "as is". It is very buggy and simply part of this example.
// Developers using this package should implement a more robust cut algorithm!
static Seam cut_mesh(Parameterization_polyhedron_adaptor* mesh_adaptor)
{
    // Helper class to compute genus or extract borders
    typedef CGAL::Parameterization_mesh_feature_extractor<Parameterization_polyhedron_adaptor>
                                            Mesh_feature_extractor;

    Seam seam;              // returned list

    // Get pointer to Polyhedron_3 mesh
    assert(mesh_adaptor != NULL);
    Polyhedron* mesh = mesh_adaptor->get_adapted_mesh();
    assert(mesh != NULL);

    // Extract mesh borders and compute genus
    Mesh_feature_extractor feature_extractor(mesh_adaptor);
    int nb_borders = feature_extractor.get_nb_borders();
    int genus = feature_extractor.get_genus();

    // If mesh is a topological disk
    if (genus == 0 && nb_borders > 0)
    {
        // Pick the longest border
        const Mesh_feature_extractor::Border* pBorder = feature_extractor.get_longest_border();
        seam = *pBorder;
    }
    else // if mesh is NOT a topological disk, create a virtual cut
    {
        const int CUT_LENGTH = 6;

        // Build consecutive halfedges array
        Polyhedron::Halfedge_handle seam_halfedges[CUT_LENGTH];
        seam_halfedges[0] = mesh->halfedges_begin();
        if (seam_halfedges[0] == NULL)
            return seam;                    // return empty list
        int i;
        for (i=1; i<CUT_LENGTH; i++)
        {
            seam_halfedges[i] = seam_halfedges[i-1]->next()->opposite()->next();
            if (seam_halfedges[i] == NULL)
                return seam;                // return empty list
        }

        // Convert halfedges array to 2-ways vertices list
        for (i=0; i<CUT_LENGTH; i++)
            seam.push_back(seam_halfedges[i]->vertex());
        for (i=CUT_LENGTH-1; i>=0; i--)
            seam.push_back(seam_halfedges[i]->opposite()->vertex());
    }

    return seam;
}

// Dump parameterized mesh to a Wavefront OBJ file
// v x y z
// f 1 2 3 4 (1-based)
//
// Implementation note: the UV is meaningless for a NON parameterized halfedge
static bool write_file_obj(Parameterization_polyhedron_adaptor* mesh_adaptor,
                           const char *pFilename)
{
    assert(mesh_adaptor != NULL);
    Polyhedron* mesh = mesh_adaptor->get_adapted_mesh();
    assert(mesh != NULL);
    assert(pFilename != NULL);

    std::ofstream out(pFilename);
    if(!out)
        return false;
    CGAL::set_ascii_mode(out);

    // Index all mesh vertices following the order of vertices_begin() iterator
    Polyhedron::Vertex_const_iterator pVertex;
    unsigned int i = 0;
    for(pVertex = mesh->vertices_begin(); pVertex != mesh->vertices_end(); pVertex++)
        mesh_adaptor->info(pVertex)->index(i++);

    // Index all mesh half edges following the order of halfedges_begin() iterator
    Polyhedron::Halfedge_const_iterator pHalfedge;
    i = 0;
    for(pHalfedge = mesh->halfedges_begin(); pHalfedge != mesh->halfedges_end(); pHalfedge++)
        mesh_adaptor->info(pHalfedge)->index(i++);

    // write the name of material file
    out <<  "mtllib parameterization.mtl" << std::endl ;

    // output coordinates
    out <<  "# vertices" << std::endl ;
    for(pVertex = mesh->vertices_begin(); pVertex != mesh->vertices_end(); pVertex++)
        out << "v " << pVertex->point().x() << " "
                    << pVertex->point().y() << " "
                    << pVertex->point().z() << std::endl;

    // Write UVs (1 UV / halfedge)
    out <<  "# uv coordinates" << std::endl ;
    for(pHalfedge = mesh->halfedges_begin(); pHalfedge != mesh->halfedges_end(); pHalfedge++)
    {
        Parameterization_polyhedron_adaptor::Halfedge_info* he_info = mesh_adaptor->info(pHalfedge);
        if (he_info->is_parameterized())
            out << "vt " << he_info->uv().x() << " " << he_info->uv().y() << std::endl;
        else
            out << "vt " << 0.0 << " " << 0.0 << std::endl;
    }

    // Write facets using the unique material # 1
    out << "# facets" << std::endl;
    out << "usemtl Mat_1" << std::endl;
    Polyhedron::Facet_const_iterator pFacet;
    for(pFacet = mesh->facets_begin(); pFacet != mesh->facets_end(); pFacet++)
    {
        Polyhedron::Halfedge_around_facet_const_circulator h = pFacet->facet_begin();
        out << "f";
        do {
            Parameterization_polyhedron_adaptor::Halfedge_info* he_info  = mesh_adaptor->info(h);
            Parameterization_polyhedron_adaptor::Vertex_info*   vtx_info = mesh_adaptor->info(h->vertex());
            out << " " << vtx_info->index()+1;
            if (he_info->is_parameterized())
                out <<  "/" << he_info->index()+1;
        }
        while(++h != pFacet->facet_begin());
        out << std::endl;
    }

    return true;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "PARAMETERIZATION" << std::endl;
    std::cerr << "  Discrete Authalic Parameterization" << std::endl;
    std::cerr << "  Square border" << std::endl;
    std::cerr << "  TAUCS solver" << std::endl;
    std::cerr << "  Very simple cut if model is not a topological disk" << std::endl;
    std::cerr << "  Output: OBJ" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 != 2)
    {
        std::cerr << "Usage: " << argv[0] << " input_file.off output_file.obj" << std::endl;
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
    if(!stream)
    {
        std::cerr << "FATAL ERROR: cannot open file " << input_filename << std::endl;
        return EXIT_FAILURE;
    }
    Polyhedron mesh;
    stream >> mesh;

    //***************************************
    // Create mesh adaptors
    //***************************************

    // The Surface_mesh_parameterization package needs an adaptor to handle Polyhedron_3 meshes
    Parameterization_polyhedron_adaptor mesh_adaptor(&mesh);

    // The parameterization methods support only meshes that
    // are topological disks => we need to compute a "cutting" of the mesh
    // that makes it it homeomorphic to a disk
    Seam seam = cut_mesh(&mesh_adaptor);
    if (seam.empty())
    {
        fprintf(stderr, "\nFATAL ERROR: an unexpected error occurred while cutting the shape!\n\n");
        return EXIT_FAILURE;
    }

    // Create adaptor that virtually "cuts" the mesh following the 'seam' path
    Mesh_patch_polyhedron   mesh_patch(&mesh_adaptor, seam.begin(), seam.end());

    //***************************************
    // Discrete Authalic Parameterization (square border)
    // with TAUCS solver
    //***************************************

    // Border parameterizer
    typedef CGAL::Square_border_arc_length_parameterizer_3<Mesh_patch_polyhedron>
                                                            Border_parameterizer;
    // TAUCS solver
    typedef CGAL::Taucs_solver_traits<double>               Solver;

    // Discrete Authalic Parameterization (square border)
    // with TAUCS solver
    typedef CGAL::Discrete_authalic_parameterizer_3<Mesh_patch_polyhedron,
                                                    Border_parameterizer,
                                                    Solver> Parameterizer;

    Parameterizer::Error_code err = CGAL::parameterize(&mesh_patch, Parameterizer());
    if (err != Parameterizer::OK)
        std::cerr << "FATAL ERROR: " << Parameterizer::get_error_message(err) << std::endl;

    //***************************************
    // Output
    //***************************************

    // Write Wavefront OBJ file
    if (err == Parameterizer::OK)
    {
        if ( ! write_file_obj(&mesh_adaptor, output_filename) )
        {
            std::cerr << "FATAL ERROR: cannot write file " << output_filename << std::endl;
            return EXIT_FAILURE;
        }
    }

    return (err == Parameterizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}


#else // CGAL_USE_TAUCS


#include <iostream>
#include <stdlib.h>
#include <stdio.h>

// ----------------------------------------------------------------------------
// Empty main() if TAUCS is not installed
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "Skip test as TAUCS is not installed" << std::endl;
    return EXIT_SUCCESS;
}


#endif // CGAL_USE_TAUCS

