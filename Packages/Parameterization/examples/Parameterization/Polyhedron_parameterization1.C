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


#include <CGAL/basic.h>

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
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

#include <OpenNL/linear_solver.h>

#include "Mesh_adaptor_polyhedron_3.h"

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

// CGAL kernel
typedef CGAL::Cartesian<double>                         Kernel;

// Mesh true type and parameterization adaptors
typedef CGAL::Polyhedron_3<Kernel>                      Polyhedron;
typedef Mesh_adaptor_polyhedron_3<Polyhedron>           Mesh_adaptor_polyhedron;
// Parametizer for this kind of mesh
typedef CGAL::Parametizer_3<Mesh_adaptor_polyhedron>    Parametizer;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "\nPARAMETERIZATION" << std::endl;


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
    fprintf(stderr, "ok\n  fill mesh...");
    stream >> mesh;


    //***************************************
    // Create mesh adaptor
    // Note: parameterization methods support only 
    // meshes that are toplogical disks
    //***************************************

    // The parameterization package needs an adaptor to handle Polyhedron_3 meshes
    Mesh_adaptor_polyhedron mesh_adaptor(&mesh);


    //***************************************
    // Floater's mean value coordinates parameterization
    //***************************************

    Parametizer::ErrorCode err = parameterize(&mesh_adaptor);
    if (err != Parametizer::OK)
        fprintf(stderr, "\nFATAL ERROR: parameterization error # %d\n", (int)err);


    ////***************************************
    //// output
    ////***************************************

    //// Write Wavefront obj file
    //if (err == Parametizer::OK)
    //    mesh.write_file_obj(output_filename);       

    cerr << std::endl;

    return (err == Parametizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}


