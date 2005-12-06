// Authalic_parameterization.C

#include "short_names.h"                    // must be included first

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>
#include <CGAL/Discrete_authalic_parameterizer_3.h>

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
typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron>     
                                                        Parameterization_polyhedron_adaptor;

// Discrete Authalic Parameterization
typedef CGAL::Discrete_authalic_parameterizer_3<Parameterization_polyhedron_adaptor>
                                                        Parameterizer;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "\nPARAMETERIZATION" << std::endl;
    std::cerr << "  Discrete Authalic Parameterization" << std::endl;
    std::cerr << "  circle border" << std::endl;
    std::cerr << "  OpenNL solver" << std::endl;


    //***************************************
    // decode parameters
    //***************************************

    if (argc-1 != 1)
    {
        std::cerr << "Usage: " << argv[0] << " input_file.off" << std::endl;
        return(EXIT_FAILURE);
    }

    // File name is:
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
    // meshes that are topological disks
    //***************************************

    // The parameterization package needs an adaptor to handle Polyhedron_3 meshes
    Parameterization_polyhedron_adaptor mesh_adaptor(&mesh);


    //***************************************
    // Discrete Authalic Parameterization
    //***************************************

    Parameterizer::Error_code err = CGAL::parameterize(&mesh_adaptor, Parameterizer());
    if (err != Parameterizer::OK)
        fprintf(stderr, "\nFATAL ERROR: parameterization error # %d\n", (int)err);


    fprintf(stderr, "\n");

    return (err == Parameterizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}


