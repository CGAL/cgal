// Polyhedron_parameterization4.C

#ifdef CGAL_USE_TAUCS


#include "short_names.h"                    // must be included first
#include <CGAL/basic.h>                     // must be included second
#include <CGAL/Cartesian.h>                 // must be included third

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Mesh_adaptor_polyhedron_3.h>
#include <CGAL/parameterize.h>

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

// Circular border parametizer (the default)
typedef CGAL::Circular_border_arc_length_parametizer_3<Mesh_adaptor_polyhedron>
                                                        Border_parametizer;
// TAUCS solver
typedef CGAL::Taucs_solver_traits<double>               Solver;

// Floater's mean value coordinates parametizer (circular border)
// with TAUCS solver
typedef CGAL::Mean_value_coordinates_parametizer_3<Mesh_adaptor_polyhedron,
                                                   Border_parametizer,
                                                   Solver>
                                                        Parametizer;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

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
#if _WIN32_WINNT >= 0x0400
    // Trick to be prompted by VisualC++ debugger when an assertion
    // fails even though we use NON debug runtime libraries
    // (the only ones compatible with TAUCS)
    if (IsDebuggerPresent())
        _set_error_mode(_OUT_TO_MSGBOX);
#endif

    std::cerr << "\nPARAMETERIZATION" << std::endl;
    std::cerr << "  Floater parameterization" << std::endl;
    std::cerr << "  circle boundary" << std::endl;
    std::cerr << "  TAUCS solver" << std::endl;
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
    // Note: parameterization methods support only
    // meshes that are topological disks
    //***************************************

    // The parameterization package needs an adaptor to handle Polyhedron_3 meshes
    Mesh_adaptor_polyhedron mesh_adaptor(&mesh);


    //***************************************
    // Floater's mean value coordinates parametizer (circular border)
    // with TAUCS solver
    //***************************************

    Parametizer::Error_code err = CGAL::parameterize(&mesh_adaptor, Parametizer());
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

