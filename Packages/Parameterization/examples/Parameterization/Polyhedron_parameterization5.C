#include "short_names.h"                    // must be included first
#include <CGAL/basic.h>                     // must be included second
#include <CGAL/Cartesian.h>                 // must be included third

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Mesh_adaptor_polyhedron_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Mesh_adaptor_patch_3.h>

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// CGAL kernel
typedef CGAL::Cartesian<double>                             Kernel;

// Mesh true type and parameterization adaptors
typedef CGAL::Polyhedron_3<Kernel>                          Polyhedron;
typedef CGAL::Mesh_adaptor_polyhedron_3<Polyhedron>         Mesh_adaptor_polyhedron;
typedef CGAL::Mesh_adaptor_patch_3<Mesh_adaptor_polyhedron> Mesh_patch_polyhedron;

// Parametizers base class for this kind of mesh
typedef CGAL::Parametizer_traits_3<Mesh_patch_polyhedron>   Parametizer;

// Type describing a border or seam as a vertex list
typedef std::list<Mesh_adaptor_polyhedron::Vertex_handle>   Seam;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// If the mesh is a topological disk, extract its longest boundary,
// else compute a very simple cut to make it homeomorphic to a disk.
// Return the border/seam (empty on error)
static Seam cut_mesh(Mesh_adaptor_polyhedron* mesh_adaptor)
{
    // Helper class to compute genus or extract boundaries
    typedef CGAL::Mesh_adaptor_feature_extractor<Mesh_adaptor_polyhedron>
                                            Mesh_feature_extractor;

    Seam seam;              // returned list

    // Get pointer to Polyhedron_3 mesh
    assert(mesh_adaptor != NULL);
    Polyhedron* mesh = mesh_adaptor->get_adapted_mesh();
    assert(mesh != NULL);

    // Extract mesh boundaries and compute genus
    Mesh_feature_extractor feature_extractor(mesh_adaptor);
    int nb_boundaries = feature_extractor.get_nb_boundaries();
    int genus = feature_extractor.get_genus();

    // If mesh is a topological disk
    if (genus == 0 && nb_boundaries > 0)
    {
        // Pick the longest boundary
        const Mesh_feature_extractor::Boundary* pBoundary = feature_extractor.get_longest_boundary();
        seam = *pBoundary;
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


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "\nPARAMETERIZATION" << std::endl;
    std::cerr << "  Floater parameterization" << std::endl;
    std::cerr << "  circle boundary" << std::endl;
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

    // File names are:
    const char* input_filename  = argv[1];


    //***************************************
    // read the mesh
    //***************************************

    fprintf(stderr, "\n  read file...%s...", input_filename);
    std::ifstream stream(input_filename);
    if(!stream) {
        fprintf(stderr, "\nFATAL ERROR: cannot open file!\n\n");
        return EXIT_FAILURE;
    }

    // read the mesh
    Polyhedron mesh;
    fprintf(stderr, "ok\n  fill mesh\n");
    stream >> mesh;


    //***************************************
    // Create mesh adaptor
    //***************************************

    // The parameterization package needs an adaptor to handle Polyhedron_3 meshes
    Mesh_adaptor_polyhedron mesh_adaptor(&mesh);

    // The parameterization package supports only meshes that
    // are toplogical disks => we need to virtually "cut" the mesh
    // to make it homeomorphic to a disk
    //
    // 1) Cut the mesh
    Seam seam = cut_mesh(&mesh_adaptor);
    if (seam.empty())
    {
        fprintf(stderr, "\nFATAL ERROR: an unexpected error occurred while cutting the shape!\n\n");
        return EXIT_FAILURE;
    }
    //
    // 2) Create adaptor that virtually "cuts" a patch in a Polyhedron_ex mesh
    Mesh_patch_polyhedron   mesh_patch(&mesh_adaptor,
                                       seam.begin(),
                                       seam.end());

    //***************************************
    // Floater's mean value coordinates parameterization
    //***************************************

    Parametizer::Error_code err = CGAL::parameterize(&mesh_patch);
    if (err != Parametizer::OK)
        fprintf(stderr, "\nFATAL ERROR: parameterization error # %d\n", (int)err);


    fprintf(stderr, "\n");

    return (err == Parametizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}


