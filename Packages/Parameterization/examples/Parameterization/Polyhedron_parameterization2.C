#include "short_names.h"                    // must be included first
#include <CGAL/basic.h>                     // must be included second
#include <CGAL/Cartesian.h>                 // must be included third

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Mesh_adaptor_polyhedron_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_authalic_parametizer_3.h>

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

// Discrete Authalic parametizer
typedef CGAL::Discrete_authalic_parametizer_3<Mesh_adaptor_polyhedron>
                                                        Parametizer;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "\nPARAMETERIZATION" << std::endl;
    std::cerr << "  Authalic parameterization" << std::endl;
    std::cerr << "  circle boundary" << std::endl;
    std::cerr << "  OpenNL solver" << std::endl;


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
    // Note: parameterization methods support only
    // meshes that are toplogical disks
    //***************************************

    // The parameterization package needs an adaptor to handle Polyhedron_3 meshes
    Mesh_adaptor_polyhedron mesh_adaptor(&mesh);


    //***************************************
    // Discrete Authalic Parameterization
    //***************************************

    Parametizer::Error_code err = CGAL::parameterize(&mesh_adaptor, Parametizer());
    if (err != Parametizer::OK)
        fprintf(stderr, "\nFATAL ERROR: parameterization error # %d\n", (int)err);


    fprintf(stderr, "\n");

    return (err == Parametizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}


