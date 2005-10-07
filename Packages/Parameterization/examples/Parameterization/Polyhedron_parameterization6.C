// Polyhedron_parameterization6.C

#ifdef CGAL_USE_TAUCS


#include "short_names.h"                    // must be included first
#include <CGAL/basic.h>                     // must be included second
#include <CGAL/Cartesian.h>                 // must be included third

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Mesh_adaptor_polyhedron_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_authalic_parametizer_3.h>
#include <CGAL/Square_border_parametizer_3.h>
#include <CGAL/Mesh_adaptor_patch_3.h>

#include <CGAL/Taucs_solver_traits.h>

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cassert>
#ifdef WIN32
    #include <Windows.h>
#endif


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// CGAL kernel
typedef CGAL::Cartesian<double>                         Kernel;

// Mesh true type and parameterization adaptors
typedef CGAL::Polyhedron_3<Kernel>                      Polyhedron;
typedef CGAL::Mesh_adaptor_polyhedron_3<Polyhedron>     Mesh_adaptor_polyhedron;
typedef CGAL::Mesh_adaptor_patch_3<Mesh_adaptor_polyhedron> Mesh_patch_polyhedron;

// Border parametizer
typedef CGAL::Square_border_arc_length_parametizer_3<Mesh_patch_polyhedron>
                                                        Border_parametizer;
// TAUCS solver
typedef CGAL::Taucs_solver_traits<double>               Solver;

// Discrete Authalic parametizer with square border
typedef CGAL::Discrete_authalic_parametizer_3<Mesh_patch_polyhedron,
                                              Border_parametizer,
                                              Solver>
                                                        Parametizer;

// Type describing a border or seam as a vertex list
typedef std::list<Mesh_adaptor_polyhedron::Vertex_handle> Seam;


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

// Dump parameterized mesh to a Wavefront OBJ file
// v x y z
// f 1 2 3 4 (1-based)
static bool write_file_obj(Mesh_adaptor_polyhedron* mesh_adaptor,
                           const char *pFilename)
{
    Polyhedron* mesh = mesh_adaptor->get_adapted_mesh();
    assert(mesh != NULL);
    assert(pFilename != NULL);

    std::cerr << "  dump mesh to " << pFilename << "..." << std::endl;
    FILE *pFile = fopen(pFilename,"wt");
    if(pFile == NULL)
    {
        std::cerr << "  unable to open file " << pFilename <<  " for writing\n";
        return false;
    }

    // Index all mesh vertices
    Polyhedron::Vertex_const_iterator pVertex;
    unsigned int i = 0;
    for(pVertex = mesh->vertices_begin(); pVertex != mesh->vertices_end(); pVertex++)
        mesh_adaptor->info(pVertex)->index(i++);

    // Index all mesh halfedges
    Polyhedron::Halfedge_const_iterator pHalfedge;
    i = 0;
    for(pHalfedge = mesh->halfedges_begin(); pHalfedge != mesh->halfedges_end(); pHalfedge++)
        mesh_adaptor->info(pHalfedge)->index(i++);

    // write the name of material file
    fprintf(pFile, "mtllib parameterization.mtl\n") ;

    // output coordinates
    fprintf(pFile, "# vertices\n") ;
    for(pVertex = mesh->vertices_begin(); pVertex != mesh->vertices_end(); pVertex++)
        fprintf(pFile,"v %g %g %g\n",
                (double)pVertex->point().x(),
                (double)pVertex->point().y(),
                (double)pVertex->point().z());

    // Write UVs (1 UV / halfedge)
    fprintf(pFile, "# uv coordinates\n") ;
    for(pHalfedge = mesh->halfedges_begin(); pHalfedge != mesh->halfedges_end(); pHalfedge++)
    {
        Mesh_adaptor_polyhedron::Halfedge_info* he_info = mesh_adaptor->info(pHalfedge);
        if (he_info->is_parameterized())
            fprintf(pFile, "vt %f %f\n", he_info->uv().x(), he_info->uv().y());
        else
            fprintf(pFile, "vt %f %f\n", 0.0, 0.0);
    }

    // Write facets using the unique material # 1
    fprintf(pFile, "# facets\nusemtl Mat_1\n");
    Polyhedron::Facet_const_iterator pFacet;
    for(pFacet = mesh->facets_begin(); pFacet != mesh->facets_end(); pFacet++)
    {
        Polyhedron::Halfedge_around_facet_const_circulator h = pFacet->facet_begin();
        fprintf(pFile,"f");
        do {
            Mesh_adaptor_polyhedron::Halfedge_info* he_info  = mesh_adaptor->info(h);
            Mesh_adaptor_polyhedron::Vertex_info*   vtx_info = mesh_adaptor->info(h->vertex());
            fprintf(pFile, " %d", (int)vtx_info->index()+1);
            if (he_info->is_parameterized())
                fprintf(pFile, "/%d", (int)he_info->index()+1);
        }
        while(++h != pFacet->facet_begin());
        fprintf(pFile,"\n");
    }

    fclose(pFile);
    return true;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
#if _WIN32_WINNT >= 0x0400
    // Trick to be prompted by VisualC++ debugger when an assertion
    // fails even though we use NON debug runtime libraries
    // (the only ones compatible with TAUCS)
    if (IsDebuggerPresent())
        _set_error_mode(_OUT_TO_MSGBOX);
#endif

    std::cerr << "\nPARAMETERIZATION" << std::endl;
    std::cerr << "  Authalic parameterization" << std::endl;
    std::cerr << "  square boundary" << std::endl;
    std::cerr << "  TAUCS solver" << std::endl;
    std::cerr << "  Very simple cut if model is not a topological disk" << std::endl;
    std::cerr << "  output: OBJ" << std::endl;


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
    // Discrete Authalic Parameterization
    // with square border.
    // TAUCS solver.
    //***************************************

    Parametizer::Error_code err = CGAL::parameterize(&mesh_patch, Parametizer());
    if (err != Parametizer::OK)
        fprintf(stderr, "\nFATAL ERROR: parameterization error # %d\n", (int)err);


    //***************************************
    // output
    //***************************************

    // Write Wavefront OBJ file
    if (err == Parametizer::OK)
        write_file_obj(&mesh_adaptor, output_filename);

    fprintf(stderr, "\n");

    return (err == Parametizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}


#else // CGAL_USE_TAUCS


#include <stdio.h>
#include <stdlib.h>

// ----------------------------------------------------------------------------
// Empty main() if TAUCS is not installed
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    fprintf(stderr, "\nSkip test as TAUCS is not installed\n\n");
    return EXIT_SUCCESS;
}


#endif // CGAL_USE_TAUCS

