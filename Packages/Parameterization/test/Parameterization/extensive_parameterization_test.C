// Copyright (c) 2005  INRIA (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$
// $Name$
//
// Author(s)     : Laurent Saboret, Pierre Alliez


// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// Test all parameterization methods and all solvers
// No output
// Input files are .off
//----------------------------------------------------------
// extensive_parameterization_test mesh1.off mesh2.off


#include "short_names.h"                    // must be included first

#include <CGAL/basic.h>                     // must be included second

#include <CGAL/parameterization.h>
#include <CGAL/Mesh_adaptor_patch_3.h>
#include <CGAL/Circular_border_parametizer_3.h>
#include <CGAL/Square_border_parametizer_3.h>
#include <CGAL/Two_vertices_parametizer_3.h>
#include <CGAL/Barycentric_mapping_parametizer_3.h>
#include <CGAL/Discrete_conformal_map_parametizer_3.h>
#include <CGAL/Discrete_authalic_parametizer_3.h>
#include <CGAL/Mean_value_coordinates_parametizer_3.h>
#include <CGAL/LSCM_parametizer_3.h>
#include <CGAL/Mesh_adaptor_feature_extractor.h>

#include <OpenNL/linear_solver.h>
#ifdef CGAL_USE_TAUCS
    #include <CGAL/Taucs_solver_traits.h>
#endif

#include "Polyhedron_ex.h"
#include "Mesh_cutter.h"
#include "Mesh_adaptor_polyhedron_ex.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <fstream>
#include <cassert>
#ifdef WIN32
    #include <Windows.h>
#endif


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// Mesh true type and parameterization adaptors
typedef Polyhedron_ex                                       Polyhedron;
typedef Mesh_adaptor_polyhedron_ex                          Mesh_adaptor_polyhedron;
typedef CGAL::Mesh_adaptor_patch_3<Mesh_adaptor_polyhedron> Mesh_patch_polyhedron;

// Parametizer for this kind of mesh
typedef CGAL::Parametizer_3<Mesh_patch_polyhedron>          Parametizer;

// Type describing a border or seam as a vertex list
typedef std::list<Mesh_adaptor_polyhedron::Vertex_handle>   Seam;

// Sparse matrix solver 1 to test: OpenNL
typedef OpenNL::DefaultLinearSolverTraits<double>           Solver1;
const char*                                                 Solver1_name = "OpenNL";

// Sparse matrix solver 2 to test: TAUCS (falls back to OpenNL if not installed)
#ifdef CGAL_USE_TAUCS
    typedef CGAL::Taucs_solver_traits<double>               Solver2;
    const char*                                             Solver2_name = "TAUCS";
#else
    typedef OpenNL::DefaultLinearSolverTraits<double>       Solver2;
    const char*                                             Solver2_name = "OpenNL";
#endif


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// Cut the mesh to make it homeomorphic to a disk
// or extract a region homeomorphic to a disc.
// Return the border of this region (empty on error)
//
// CAUTION:
// This method is provided "as is". It is very buggy and simply part of this example.
// Developers using this package should implement a more robust cut algorithm!
static Seam cut_mesh(Mesh_adaptor_polyhedron* mesh_adaptor)
{
    // Type describing a border or seam as an halfedge list
    typedef CGAL::Mesh_adaptor_feature_extractor<Mesh_adaptor_polyhedron_ex>
                                            Mesh_feature_extractor;
    typedef Mesh_feature_extractor::Boundary Boundary;
    typedef Mesh_cutter::Backbone           Backbone;

    Seam seam;              // returned list

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
        const Boundary* pBoundary = feature_extractor.get_longest_boundary();
        seam = *pBoundary;
    }
    else // if mesh is NOT a topological disk
    {
        Backbone seamingBackbone;           // result of cutting
        Backbone::iterator he;

        // Virtually "cut" mesh to make it a topological disk
        mesh->compute_facet_centers();
        Mesh_cutter cutter(mesh);
        if (genus == 0)
        {
            // no boundary, we need to cut the mesh
            assert (nb_boundaries == 0);
            cutter.cut(&seamingBackbone);   // simple cut
        }
        else // genus > 0 -> cut the mesh
        {
            cutter.cut_genus(&seamingBackbone);
        }

        // The Mesh_cutter class is quite buggy
        // => we check that seamingBackbone is valid
        //
#ifdef DEBUG_TRACE
        // Dump seam (for debug purpose)
        mesh->precompute_halfedge_indices();
        mesh->precompute_vertex_indices();
        fprintf(stderr,"  HE seam is is: ");
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
        {
          fprintf(stderr, "H%d=%d->%d ",
                          (int)(*he)->index(),
                          (int)(*he)->opposite()->vertex()->index(),
                          (int)(*he)->vertex()->index());
        }
        fprintf(stderr,"ok\n");
#endif
        //
        // 1) Check that seamingBackbone is not empty
        if (seamingBackbone.begin() == seamingBackbone.end())
            return seam;                    // return empty list
        //
        // 2) Check that seamingBackbone is a loop and
        //    count occurences of seam halfedges
        mesh->tag_halfedges(0);             // Reset counters
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
        {
            // Get next halfedge iterator (looping)
            Backbone::iterator next_he = he;
            next_he++;
            if (next_he == seamingBackbone.end())
                next_he = seamingBackbone.begin();

            // Check that seamingBackbone is a loop: check that
            // end of current HE == start of next one
            if ((*he)->vertex() != (*next_he)->opposite()->vertex())
                return seam;                // return empty list

            // Increment counter (in "tag" field) of seam halfedges
            (*he)->tag( (*he)->tag()+1 );
        }
        //
        // 3) check that the seamingBackbone is a 2-way list
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
        {
            // Counter of halfedge and opposite halfedge must be 1
            if ((*he)->tag() != 1 || (*he)->opposite()->tag() != 1)
                return seam;                // return empty list
        }

        // Convert list of halfedges to a list of vertices
        for (he = seamingBackbone.begin(); he != seamingBackbone.end(); he++)
            seam.push_back((*he)->vertex());
    }

    return seam;
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
    std::cerr << "Test all parameterization methods and all solvers" << std::endl;
    std::cerr << "No output" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 == 0)
    {
        std::cerr << "Usage: " << argv[0] << " input_file1.off input_file2.obj ..." << std::endl;
        return(EXIT_FAILURE);
    }

    // Accumulated errors
    int accumulated_err = EXIT_SUCCESS;

    // Parameterize each input file and accumulate errors
    for (int arg_index = 1; arg_index <= argc-1; arg_index++)
    {
        // File name is:
        const char* input_filename  = argv[arg_index];

        //***************************************
        // read the mesh
        //***************************************

        fprintf(stderr, "\nRead file...%s...", input_filename);
        std::ifstream stream(input_filename);
        if(!stream) {
            fprintf(stderr, "\nFATAL ERROR: cannot open file!\n\n");
            accumulated_err = EXIT_FAILURE;
            continue;
        }

        // read the mesh
        Polyhedron mesh;
        fprintf(stderr, "ok\n  fill mesh...");
        stream >> mesh;

        // print mesh info
        fprintf(stderr, "(%d facets, ",mesh.size_of_facets());
        fprintf(stderr, "%d vertices)\n",mesh.size_of_vertices());

        //***************************************
        // Create mesh adaptor
        //***************************************

        // The parameterization package needs an adaptor to handle Polyhedron_ex meshes
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
            accumulated_err = EXIT_FAILURE;
            continue;
        }
        //
        // 2) Create adaptor that virtually "cuts" a patch in a Polyhedron_ex mesh
        Mesh_patch_polyhedron   mesh_patch(&mesh_adaptor,
                                        seam.begin(),
                                        seam.end());

        std::cerr << std::endl;

        //***************************************
        // Tutte uniform parameterization
        // with square uniform boundary parameterization
        // OpenNL solver
        //***************************************

        Parametizer::Error_code err;

        std::cerr << "Tutte uniform parameterization" << std::endl;
        std::cerr << "with square uniform boundary parameterization" << std::endl;
        std::cerr << Solver1_name << " solver" << std::endl;

        err = CGAL::parameterize(
            &mesh_patch,
            CGAL::Barycentric_mapping_parametizer_3<
                Mesh_patch_polyhedron,
                CGAL::Square_border_uniform_parametizer_3<Mesh_patch_polyhedron>,
                Solver1
            >());

        if (err == Parametizer::OK)
            fprintf(stderr, "  Parameterization success\n\n");
        else
            fprintf(stderr, "\nMINOR ERROR: parameterization error # %d\n\n", (int)err);

        //***************************************
        // Floater's mean value coordinates parameterization
        // with circular arc length boundary parameterization
        // TAUCS solver (if installed)
        //***************************************

        std::cerr << "Floater's mean value coordinates parameterization" << std::endl;
        std::cerr << "with circular arc length boundary parameterization" << std::endl;
        std::cerr << Solver2_name << " solver" << std::endl;

        err = CGAL::parameterize(
            &mesh_patch,
            CGAL::Mean_value_coordinates_parametizer_3<
                Mesh_patch_polyhedron,
                CGAL::Circular_border_arc_length_parametizer_3<Mesh_patch_polyhedron>,
                Solver2
            >());

        if (err == Parametizer::OK)
            fprintf(stderr, "  Parameterization success\n\n");
        else
            fprintf(stderr, "\nMINOR ERROR: parameterization error # %d\n\n", (int)err);

        //***************************************
        // Discrete Conformal Map parameterization
        // with circular arc length boundary parameterization
        // TAUCS solver (if installed)
        //***************************************

        std::cerr << "Discrete Conformal Map parameterization" << std::endl;
        std::cerr << "with circular arc length boundary parameterization" << std::endl;
        std::cerr << Solver2_name << " solver" << std::endl;

        err = CGAL::parameterize(
            &mesh_patch,
            CGAL::Discrete_conformal_map_parametizer_3<
                Mesh_patch_polyhedron,
                CGAL::Circular_border_arc_length_parametizer_3<Mesh_patch_polyhedron>,
                Solver2
            >());

        if (err == Parametizer::OK)
            fprintf(stderr, "  Parameterization success\n\n");
        else
            fprintf(stderr, "\nMINOR ERROR: parameterization error # %d\n\n", (int)err);

        //***************************************
        // Authalic parameterization
        // with square arc length boundary parameterization
        // TAUCS solver (if installed)
        //***************************************

        std::cerr << "Authalic parameterization" << std::endl;
        std::cerr << "with square arc length boundary parameterization" << std::endl;
        std::cerr << Solver2_name << " solver" << std::endl;

        err = CGAL::parameterize(
            &mesh_patch,
            CGAL::Discrete_authalic_parametizer_3<
                Mesh_patch_polyhedron,
                CGAL::Square_border_arc_length_parametizer_3<Mesh_patch_polyhedron>,
                Solver2
            >());

        if (err == Parametizer::OK)
            fprintf(stderr, "  Parameterization success\n\n");
        else
            fprintf(stderr, "\nMINOR ERROR: parameterization error # %d\n\n", (int)err);

        //***************************************
        // LSCM parameterization
        // TAUCS solver (if installed)
        //***************************************

        std::cerr << "LSCM parameterization" << std::endl;
        std::cerr << Solver2_name << " solver" << std::endl;

        err = CGAL::parameterize(
            &mesh_patch,
            CGAL::LSCM_parametizer_3<
                Mesh_patch_polyhedron,
                CGAL::Two_vertices_parametizer_3<Mesh_patch_polyhedron>,
                Solver2
            >());

        if (err == Parametizer::OK)
            fprintf(stderr, "  Parameterization success\n\n");
        else
            fprintf(stderr, "\nMINOR ERROR: parameterization error # %d\n\n", (int)err);

    } // for each input file

    // Return accumulated fatal error
    fprintf(stderr, "\nTool returns %d\n\n", (int)accumulated_err);
    return accumulated_err;
}


