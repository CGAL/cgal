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
// Floater parameterization
// circle boundary
// OpenNL solver
// Very simple cut if model is not a topological disk
// output is a eps map
// input file is mesh.off
//----------------------------------------------------------
// Polyhedron_parameterization5 mesh.off mesh.eps


#include "short_names.h"                    // must be included first

#include <CGAL/basic.h>                     // must be included second

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Mesh_adaptor_polyhedron_3.h>
#include <CGAL/parameterization.h>
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
typedef CGAL::Cartesian<double>                         Kernel;

// Mesh true type and parameterization adaptors
typedef CGAL::Polyhedron_3<Kernel>                      Polyhedron;
typedef CGAL::Mesh_adaptor_polyhedron_3<Polyhedron>     Mesh_adaptor_polyhedron;
typedef CGAL::Mesh_adaptor_patch_3<Mesh_adaptor_polyhedron> Mesh_patch_polyhedron;

// Parametizers base class for this kind of mesh
typedef CGAL::Parametizer_3<Mesh_patch_polyhedron>      Parametizer;

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
    // Type describing a border or seam as an halfedge list
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

// Dump parameterized mesh to an eps file
static bool write_file_eps(const Mesh_adaptor_polyhedron& mesh_adaptor, 
                           const char *pFilename, 
                           double scale = 500.0)
{
    const Polyhedron* mesh = mesh_adaptor.get_adapted_mesh();
    assert(mesh != NULL);
    assert(pFilename != NULL);

    std::cerr << "  dump mesh to " << pFilename << "..." << std::endl;
    FILE *pFile = fopen(pFilename,"wt");
    if(pFile == NULL)
    {
        std::cerr << "  unable to open file " << pFilename <<  " for writing" << std::endl;
        return false;
    }

    // compute bounding box
    double xmin,xmax,ymin,ymax;
    xmin = ymin = xmax = ymax = 0;
    Polyhedron::Halfedge_const_iterator pHalfedge;
    for (pHalfedge = mesh->halfedges_begin();
         pHalfedge != mesh->halfedges_end();
         pHalfedge++)
    {
        double x1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().x();
        double y1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().y();
        double x2 = scale * mesh_adaptor.info(pHalfedge)->uv().x();
        double y2 = scale * mesh_adaptor.info(pHalfedge)->uv().y();
        xmin = std::min(xmin,x1);
        xmin = std::min(xmin,x2);
        xmax = std::max(xmax,x1);
        xmax = std::max(xmax,x2);
        ymax = std::max(ymax,y1);
        ymax = std::max(ymax,y2);
        ymin = std::min(ymin,y1);
        ymin = std::min(ymin,y2);
    }

    fprintf(pFile,"%%!PS-Adobe-2.0 EPSF-2.0\n");
    fprintf(pFile,"%%%%BoundingBox: %d %d %d %d\n", int(xmin+0.5), int(ymin+0.5), int(xmax+0.5), int(ymax+0.5));
    fprintf(pFile,"%%%%HiResBoundingBox: %g %g %g %g\n",xmin,ymin,xmax,ymax);
    fprintf(pFile,"%%%%EndComments\n");
    fprintf(pFile,"gsave\n");
    fprintf(pFile,"0.1 setlinewidth\n");

    // color macros
    fprintf(pFile,"\n%% RGB color command - r g b C\n");
    fprintf(pFile,"/C { setrgbcolor } bind def\n");
    fprintf(pFile,"/white { 1 1 1 C } bind def\n");
    fprintf(pFile,"/black { 0 0 0 C } bind def\n");

    // edge macro -> E
    fprintf(pFile,"\n%% Black stroke - x1 y1 x2 y2 E\n");
    fprintf(pFile,"/E {moveto lineto stroke} bind def\n");
    fprintf(pFile,"black\n\n");

    // for each halfedge
    for (pHalfedge = mesh->halfedges_begin();
         pHalfedge != mesh->halfedges_end();
         pHalfedge++)
    {
        double x1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().x();
        double y1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().y();
        double x2 = scale * mesh_adaptor.info(pHalfedge)->uv().x();
        double y2 = scale * mesh_adaptor.info(pHalfedge)->uv().y();
        fprintf(pFile,"%g %g %g %g E\n",x1,y1,x2,y2);
    }

    /* Emit EPS trailer. */
    fputs("grestore\n\n",pFile);
    fputs("showpage\n",pFile);

    fclose(pFile);
    return true;
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
    std::cerr << "  output: EPS" << std::endl;


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


    //***************************************
    // output
    //***************************************

    // Write Postscript file
    if (err == Parametizer::OK)
        write_file_eps(mesh_adaptor, output_filename);       

    fprintf(stderr, "\n");

    return (err == Parametizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}


