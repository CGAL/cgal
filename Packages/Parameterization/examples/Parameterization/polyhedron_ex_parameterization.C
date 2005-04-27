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
// Author(s)     : Laurent Saboret, Bruno Levy, Pierre Alliez


// ----------------------------------------------------------------------------
// USAGE EXAMPLES
// ----------------------------------------------------------------------------

//----------------------------------------------------------
// conformal parameterization
// circle boundary
// output is a ps map
// input file is mesh.off (the latter must be at the very end)
//----------------------------------------------------------
// polyhedron_ex_parameterization -t conformal -b circle -o eps mesh.off mesh.eps

//----------------------------------------------------------
// floater parameterization
// square boundary
// output is a ps map
// input file is mesh.off
//----------------------------------------------------------
// polyhedron_ex_parameterization -t floater -b square -o eps mesh.off mesh.eps

//----------------------------------------------------------
// natural parameterization
// no explicitly pinned vertices
// output is a ps map
// input file is mesh.off
//----------------------------------------------------------
// polyhedron_ex_parameterization -t natural -o eps mesh.off mesh.eps

//----------------------------------------------------------
// LSCM parameterization
// no explicitly pinned vertices
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
    Backbone *pSeamingBackbone = mesh->get_seaming_backbone();
    mesh->free_skeleton();
    pSeamingBackbone->clear();
    //mesh->flag_halfedges_seaming(false);

    Mesh_cutter cutter(mesh);
    // cutter.keep_one_connected_component();
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
            //cutter.convert_to_seaming(pSeamingBackbone);

            // cleanup this one (has to be done later if has corners)
            mesh->free_skeleton();
        }
    }
    else // genus > 0 -> cut the mesh
    cutter.cut_genus(pSeamingBackbone);

    return pSeamingBackbone;
}


// ----------------------------------------------------------------------------
// main
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
    // -t    conformal   -> default
    //       natural     -> free boundaries
    //       floater     -> mean coordinate values
    //       uniform     -> weight = 1
    //       authalic    -> weak area-preserving
    //       lscm        -> Least Squares Conformal Maps


    "b:boundary <string>", // -b or --boundary (for fixed border param.)
    // -b circle        -> map mesh boundary onto a circle
    //    square        -> map mesh boundary onto a square

    "o:output <string>", // -o or --output
    // -o eps   -> eps map       (-o eps file.off > file.eps)
    //    obj   -> Wavefront obj (-o obj file.off > file.obj)

    NULL
} ;

// Parameters description for usage
static const char * usage = "off-input-file [output-file]\n\
where type is:     conformal, natural, floater, uniform, authalic or lscm\n\
      boundary is: circle or square\n\
      output is:   eps or obj";

// Border parameterization methods
enum boundary_type {BOUNDARY_CIRCLE, BOUNDARY_SQUARE};

int main(int argc,char * argv[])
{
    std::cerr << "\nPARAMETERIZATION" << std::endl;

    // options
    const char *type = "conformal";  // default: conformal param
    const char *boundary = "circle"; // default: circle boundary
    boundary_type boundary_type = BOUNDARY_CIRCLE;
    const char *output = "";         // default:  output nothing

    // misc
    char  optchar;
    const char *optarg;
    int  errors = 0;
    int  npos = 0;

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
            if(strcmp(boundary,"circle") == 0)
            boundary_type = BOUNDARY_CIRCLE;
            else
            boundary_type = BOUNDARY_SQUARE;
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
    // check input file
    //***************************************
    fprintf(stderr,"\n  read file...%s...", input_filename);
    std::ifstream stream(input_filename);
    if(!stream) {
        fprintf(stderr,"failed: cannot open file!\n");
        return 1;
    }

    //***************************************
    // read the mesh
    //***************************************
    Polyhedron_ex mesh;
    fprintf(stderr,"ok\n  fill mesh...");
    stream >> mesh;

    // print mesh info
    fprintf(stderr,"(%d faces, ",mesh.size_of_facets());
    fprintf(stderr,"%d vertices)\n",mesh.size_of_vertices());

    // compute misc.
    //mesh.compute_normals();
    mesh.compute_facet_centers();
    //mesh.compute_mean_curvature_normal();

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

    // Switch parameterization. Defaults are:
    // - floater mean value coordinates
    // - fixing border on a circle (for fixed parameterization)
    // - OpenNL sparse linear solver
    CGAL::Parametizer_3<Mesh_adaptor_polyhedron_ex>::ErrorCode err;
    if ( (strcmp(type,"floater") == 0) && (strcmp(boundary,"circle") == 0) )
    {
        // Floater/circle is the default parameterization algorithm
        err = CGAL::parameterize(&mesh_adaptor);
    }
    else if ( (strcmp(type,"floater") == 0) && (strcmp(boundary,"square") == 0) )
    {
        err = CGAL::parameterize(
            &mesh_adaptor,
            CGAL::Mean_value_coordinates_parametizer_3<Mesh_adaptor_polyhedron_ex,
            CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
    }
    else if ( (strcmp(type,"uniform") == 0) && (strcmp(boundary,"circle") == 0) )
    {
        err = CGAL::parameterize(
            &mesh_adaptor,
            CGAL::Barycentric_mapping_parametizer_3<Mesh_adaptor_polyhedron_ex>());
    }
    else if ( (strcmp(type,"uniform") == 0) && (strcmp(boundary,"square") == 0) )
    {
        err = CGAL::parameterize(
            &mesh_adaptor,
            CGAL::Barycentric_mapping_parametizer_3<Mesh_adaptor_polyhedron_ex,
            CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
    }
    else if ( (strcmp(type,"conformal") == 0) && (strcmp(boundary,"circle") == 0) )
    {
        err = CGAL::parameterize(
            &mesh_adaptor,
            CGAL::Discrete_conformal_map_parametizer_3<Mesh_adaptor_polyhedron_ex>());
    }
    else if ( (strcmp(type,"conformal") == 0) && (strcmp(boundary,"square") == 0) )
    {
        err = CGAL::parameterize(
            &mesh_adaptor,
            CGAL::Discrete_conformal_map_parametizer_3<Mesh_adaptor_polyhedron_ex,
            CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
    }
    else if ( (strcmp(type,"authalic") == 0) && (strcmp(boundary,"circle") == 0) )
    {
        err = CGAL::parameterize(
            &mesh_adaptor,
            CGAL::Discrete_authalic_parametizer_3<Mesh_adaptor_polyhedron_ex>());
    }
    else if ( (strcmp(type,"authalic") == 0) && (strcmp(boundary,"square") == 0) )
    {
        err = CGAL::parameterize(
            &mesh_adaptor,
            CGAL::Discrete_authalic_parametizer_3<Mesh_adaptor_polyhedron_ex,
            CGAL::Square_border_parametizer_3<Mesh_adaptor_polyhedron_ex> >());
    }
    else if (strcmp(type,"lscm") == 0)
    {
        err = CGAL::parameterize(
            &mesh_adaptor,
            CGAL::LSCM_parametizer_3<Mesh_adaptor_polyhedron_ex>());
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
        mesh.dump_param(output_filename);               // write EPS file
    else if(strcmp(output,"obj") == 0)
        mesh.write_file_obj(output_filename);           // write Wavefront obj file

    // Flush cout and cerr
    cout << std::endl;
    cerr << std::endl;

    return EXIT_SUCCESS;
}

