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
// conformal parameterization
// circle boundary
// OpenNL solver
// output is a ps map
// input file is mesh.off (the latter must be at the very end)
//----------------------------------------------------------
// polyhedron_ex_parameterization -t conformal -b circle -o eps mesh.off mesh.eps

//----------------------------------------------------------
// floater parameterization
// square boundary
// TAUCS solver
// output is a ps map
// input file is mesh.off
//----------------------------------------------------------
// polyhedron_ex_parameterization -t floater -b square -s taucs -o eps mesh.off mesh.eps

//----------------------------------------------------------
// natural parameterization
// no explicitly pinned vertices
// OpenNL solver
// output is a ps map
// input file is mesh.off
//----------------------------------------------------------
// polyhedron_ex_parameterization -t natural -s opennl -o eps mesh.off mesh.eps

//----------------------------------------------------------
// LSCM parameterization
// 2 pinned vertices (automatically picked)
// OpenNL solver
// output is a .obj
// input file is mesh.off
//----------------------------------------------------------
// polyhedron_ex_parameterization -t lscm -o obj mesh.off mesh.obj


#include <CGAL/basic.h>

#include <CGAL/parameterization.h>
#include <CGAL/Circular_border_parametizer_3.h>
#include <CGAL/Square_border_parametizer_3.h>
#include <CGAL/Two_vertices_parametizer_3.h>
#include <CGAL/Barycentric_mapping_parametizer_3.h>
#include <CGAL/Discrete_conformal_map_parametizer_3.h>
#include <CGAL/Discrete_authalic_parametizer_3.h>
#include <CGAL/Mean_value_coordinates_parametizer_3.h>
#include <CGAL/LSCM_parametizer_3.h>

#include "options.h"
#include "cgal_types.h"
#include "Feature_skeleton.h"
#include "Mesh_cutter.h"
#include "Mesh_feature_extractor.h"
#include "Mesh_adaptor_polyhedron_ex.h"

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <fstream>
#include <cassert>

#include <OpenNL/linear_solver.h>
#include <CGAL/Taucs_solver_traits.h>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

// Type describing a seam in a Polyhedron_ex mesh
typedef Feature_backbone<Polyhedron_ex::Vertex_handle,
                         Polyhedron_ex::Halfedge_handle> Backbone;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// Cut the mesh to make it homeomorphic to a disk
// or extract a region homeomorphic to a disc.
// Return the border of this region.
static const Backbone* cut_mesh(Polyhedron_ex* mesh)
{
    // init
    mesh->compute_facet_centers();
    Backbone *pSeamingBackbone = mesh->get_seaming_backbone();
    mesh->free_skeleton();
    pSeamingBackbone->clear();

    Mesh_cutter cutter(mesh);
    Mesh_feature_extractor feature_extractor(mesh);

    // compute genus
    int genus = mesh->genus();
    if(genus == 0)
    {
        int nb_boundaries = feature_extractor.extract_boundaries(true);

        // no boundary, we need to cut the mesh
        if(nb_boundaries == 0)
            cutter.cut(pSeamingBackbone); // simple cut
        else
        {
            // must have one boundary, pick the first_file_arg
            Backbone *pBackbone = (*mesh->get_skeleton()->backbones())[0];
            pSeamingBackbone->copy_from(pBackbone);

            // cleanup this one (has to be done later if has corners)
            mesh->free_skeleton();
        }
    }
    else // genus > 0 -> cut the mesh
    cutter.cut_genus(pSeamingBackbone);

    return pSeamingBackbone;
}

// Switch parameterization
template<class SparseLinearAlgebraTraits_d> 
                                    // Traits class to solve a sparse linear system
CGAL::Parametizer_3<Mesh_adaptor_polyhedron_ex>::ErrorCode
parameterize(Mesh_adaptor_polyhedron_ex* mesh_adaptor,
                                    // Mesh parameterization adaptor
             const char *type,      // type of parameterization (see usage)
             const char *boundary)  // type of boundary parameterization (see usage)
{
    CGAL::Parametizer_3<Mesh_adaptor_polyhedron_ex>::ErrorCode err;

    if ( (strcmp(type,"floater") == 0) && (strcmp(boundary,"circle") == 0) )
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::Mean_value_coordinates_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Circular_border_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (strcmp(type,"floater") == 0) && (strcmp(boundary,"square") == 0) )
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::Mean_value_coordinates_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (strcmp(type,"uniform") == 0) && (strcmp(boundary,"circle") == 0) )
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::Barycentric_mapping_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Circular_border_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (strcmp(type,"uniform") == 0) && (strcmp(boundary,"square") == 0) )
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::Barycentric_mapping_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (strcmp(type,"conformal") == 0) && (strcmp(boundary,"circle") == 0) )
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::Discrete_conformal_map_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Circular_border_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (strcmp(type,"conformal") == 0) && (strcmp(boundary,"square") == 0) )
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::Discrete_conformal_map_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (strcmp(type,"authalic") == 0) && (strcmp(boundary,"circle") == 0) )
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::Discrete_authalic_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Circular_border_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if ( (strcmp(type,"authalic") == 0) && (strcmp(boundary,"square") == 0) )
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::Discrete_authalic_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else if (strcmp(type,"lscm") == 0)
    {
        err = CGAL::parameterize(
            mesh_adaptor,
            CGAL::LSCM_parametizer_3<
                Mesh_adaptor_polyhedron_ex,
                CGAL::Two_vertices_parametizer_3<Mesh_adaptor_polyhedron_ex>,
                SparseLinearAlgebraTraits_d
            >());
    }
    else
    {
        fprintf(stderr,"invalid choice\n");
        exit(1);
    }

    return err;
}


// ----------------------------------------------------------------------------
// Usage
// ----------------------------------------------------------------------------

// Parameters description for Options library
static const char *  optv[] =
{
    //     '|' -- indicates that the option takes NO argument;
    //     '?' -- indicates that the option takes an OPTIONAL argument;
    //     ':' -- indicates that the option takes a REQUIRED argument;
    //     '*' -- indicates that the option takes 0 or more arguments;
    //     '+' -- indicates that the option takes 1 or more arguments;

    "t:type <string>", // -t or --type
    // -t    floater     -> mean coordinate values (default)
    //       conformal   -> conformal
    //       natural     -> free boundaries
    //       uniform     -> weight = 1
    //       authalic    -> weak area-preserving
    //       lscm        -> Least Squares Conformal Maps


    "b:boundary <string>", // -b or --boundary (for fixed border param.)
    // -b circle        -> map mesh boundary onto a circle
    //    square        -> map mesh boundary onto a square

    "s:solver <string>", // -s or --solver
    // -s opennl        -> OpenNL solver (GPL)
    //    taucs         -> TAUCS solver

    "o:output <string>", // -o or --output
    // -o eps   -> eps map       (-o eps file.off > file.eps)
    //    obj   -> Wavefront obj (-o obj file.off > file.obj)

    NULL
} ;

// Parameters description for usage
static const char * usage = "off-input-file [output-file]\n\
where type is:     floater (default), conformal, natural, uniform, authalic or lscm\n\
      boundary is: circle (default) or square\n\
      solver is:   opennl (default) or taucs\n\
      output is:   eps or obj (default is no output)";


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "\nPARAMETERIZATION" << std::endl;

    // options
    const char *type = "floater";       // default: Floater param
    const char *boundary = "circle";    // default: circular boundary param.
    const char *solver = "opennl";      // default: OpenNL solver
    const char *output = "";            // default: no output 

    // misc
    char  optchar;
    const char *optarg;
    int  errors = 0;
    int  npos = 0;

    //***************************************
    // decode parameters
    //***************************************

    Options opts(*argv, optv);
    OptArgvIter iter(--argc,++argv);
    while ((optchar = opts(iter,optarg)) != 0)
    {
        switch(optchar)
        {
        // type of param.
        case 't' :
            type = optarg;
            std::cerr << "  " << type << " parameterization" << std::endl;
            break;

        // output: "eps" or "obj"
        case 'o' :
            output = optarg;
            std::cerr << "  " << "output: " << output << std::endl;
            break;

        // boundary parameterization (for fixed border algorithms)
        case 'b' :
            boundary = optarg;
            std::cerr << "  " << boundary << " boundary" << std::endl;
            break;

        // solver
        case 's' :
            solver = optarg;
            std::cerr << "  " << solver << " solver" << std::endl;
            break;

        // help
        case '?' :
        case 'H' :
            opts.usage(cerr, usage);
            ::exit(0);
            break;

        default :
            ++errors;
            break;
        } // end switch
    }

    // Get file name arguments
    int index = iter.index();
    int first_file_arg = 0;             // index of 1st file name argument
    int file_args_end   = 0;            // index of last file name argument + 1
    if ((npos > 0) || (index < argc))
    {
        first_file_arg = (npos > 0) ? 0    : index;
        file_args_end   = (npos > 0) ? npos : argc;
    }
    int nb_filename_arguments = file_args_end - first_file_arg;
    assert(nb_filename_arguments >= 0);

    // check options
    int nb_filenames_needed = 1 /* input*/ + (strlen(output) > 0) /* output? */;
    if (errors || nb_filename_arguments != nb_filenames_needed)
    {
        opts.usage(cerr, usage);
        ::exit(1);
    }

    // File names are:
    const char* input_filename  = argv[first_file_arg];
    const char* output_filename = (strlen(output) > 0) ? argv[first_file_arg+1]
                                                       : NULL;

    //***************************************
    // read the mesh
    //***************************************

    fprintf(stderr,"\n  read file...%s...", input_filename);
    std::ifstream stream(input_filename);
    if(!stream) {
        fprintf(stderr,"failed: cannot open file!\n");
        return 1;
    }

    // read the mesh
    Polyhedron_ex mesh;
    fprintf(stderr,"ok\n  fill mesh...");
    stream >> mesh;

    // print mesh info
    fprintf(stderr,"(%d facets, ",mesh.size_of_facets());
    fprintf(stderr,"%d vertices)\n",mesh.size_of_vertices());

    //***************************************
    // switch parameterization
    //***************************************

    // Cut the mesh to make it homeomorphic to a disk
    // or extract a region homeomorphic to a disc
    const Backbone* seam = cut_mesh(&mesh);
    assert(seam != NULL);

    // The parameterization package needs an adaptor to handle Polyhedrons
    Mesh_adaptor_polyhedron_ex mesh_adaptor(&mesh,
                                            seam->halfedges()->begin(),
                                            seam->halfedges()->end());

    // Switch parameterization
    CGAL::Parametizer_3<Mesh_adaptor_polyhedron_ex>::ErrorCode err;
    if (strcmp(solver,"opennl") == 0) 
    {
        err = parameterize<OpenNL::DefaultLinearSolverTraits<double> >(&mesh_adaptor, type, boundary);
    }
    else if (strcmp(solver,"taucs") == 0)
    {
        err = parameterize<CGAL::Taucs_solver_traits<double> >(&mesh_adaptor, type, boundary);
    }
    else
    {
        fprintf(stderr,"invalid choice\n");
        return 1;
    }

    // On parameterization error
    if (err != CGAL::Parametizer_3<Mesh_adaptor_polyhedron_ex>::OK)
    {
        fprintf(stderr,"parameterization failure: error = %d\n", (int)err);
        return err;
    }

    //***************************************
    // output
    //***************************************

    if(strcmp(output,"eps") == 0)
    {
        mesh.dump_param(output_filename);               // write EPS file
    }
    else if(strcmp(output,"obj") == 0)
    {
        mesh.write_file_obj(output_filename);           // write Wavefront obj file
    }
    else
    {
        fprintf(stderr,"invalid choice\n");
        return 1;
    }

    //***************************************
    // cleanup
    //***************************************

    // Flush cout and cerr
    cout << std::endl;
    cerr << std::endl;

    return EXIT_SUCCESS;
}


