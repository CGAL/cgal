// Taucs_parameterization.C

#ifdef CGAL_USE_TAUCS


#include "short_names.h"                    // must be included first

#include <CGAL/Cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Parameterization_polyhedron_adaptor_3.h>
#include <CGAL/parameterize.h>

#include <CGAL/Taucs_solver_traits.h>

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <fstream>
#include <cassert>


// ----------------------------------------------------------------------------
// Private types
// ----------------------------------------------------------------------------

typedef CGAL::Cartesian<double>                         Kernel;
typedef CGAL::Polyhedron_3<Kernel>                      Polyhedron;

// Mesh adaptor
typedef CGAL::Parameterization_polyhedron_adaptor_3<Polyhedron>
                                                        Parameterization_polyhedron_adaptor;


// ----------------------------------------------------------------------------
// Private functions
// ----------------------------------------------------------------------------

// Dump parameterized mesh to an eps file
static bool write_file_eps(const Parameterization_polyhedron_adaptor& mesh_adaptor,
                           const char *pFilename,
                           double scale = 500.0)
{
    const Polyhedron* mesh = mesh_adaptor.get_adapted_mesh();
    assert(mesh != NULL);
    assert(pFilename != NULL);

    std::ofstream out(pFilename);
    if(!out) 
        return false;
    CGAL::set_ascii_mode(out);

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

    out << "%!PS-Adobe-2.0 EPSF-2.0" << std::endl;
    out << "%%BoundingBox: " << int(xmin+0.5) << " " 
                                << int(ymin+0.5) << " " 
                                << int(xmax+0.5) << " " 
                                << int(ymax+0.5) << std::endl;
    out << "%%HiResBoundingBox: " << xmin << " " 
                                    << ymin << " " 
                                    << xmax << " " 
                                    << ymax << std::endl;
    out << "%%EndComments" << std::endl;
    out << "gsave" << std::endl;
    out << "0.1 setlinewidth" << std::endl;

    // color macros
    out << std::endl;
    out << "% RGB color command - r g b C" << std::endl;
    out << "/C { setrgbcolor } bind def" << std::endl;
    out << "/white { 1 1 1 C } bind def" << std::endl;
    out << "/black { 0 0 0 C } bind def" << std::endl;

    // edge macro -> E
    out << std::endl;
    out << "% Black stroke - x1 y1 x2 y2 E" << std::endl;
    out << "/E {moveto lineto stroke} bind def" << std::endl;
    out << "black" << std::endl << std::endl;

    // for each halfedge
    for (pHalfedge = mesh->halfedges_begin();
         pHalfedge != mesh->halfedges_end();
         pHalfedge++)
    {
        double x1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().x();
        double y1 = scale * mesh_adaptor.info(pHalfedge->prev())->uv().y();
        double x2 = scale * mesh_adaptor.info(pHalfedge)->uv().x();
        double y2 = scale * mesh_adaptor.info(pHalfedge)->uv().y();
        out << x1 << " " << y1 << " " << x2 << " " << y2 << " E" << std::endl;
    }

    /* Emit EPS trailer. */
    out << "grestore" << std::endl;
    out << std::endl;
    out << "showpage" << std::endl;

    return true;
}


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "PARAMETERIZATION" << std::endl;
    std::cerr << "  Floater parameterization" << std::endl;
    std::cerr << "  Circle border" << std::endl;
    std::cerr << "  TAUCS solver" << std::endl;
    std::cerr << "  Output: EPS" << std::endl;

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
    // Read the mesh
    //***************************************

    // Read the mesh
    std::ifstream stream(input_filename);
    if(!stream) 
    {
        std::cerr << "FATAL ERROR: cannot open file " << input_filename << std::endl;
        return EXIT_FAILURE;
    }
    Polyhedron mesh;
    stream >> mesh;

    //***************************************
    // Create mesh adaptor
    // Note: parameterization methods support only
    // meshes that are topological disks
    //***************************************

    // The parameterization package needs an adaptor to handle Polyhedron_3 meshes
    Parameterization_polyhedron_adaptor mesh_adaptor(&mesh);

    //***************************************
    // Floater Mean Value Coordinates parameterizer (circular border)
    // with TAUCS solver
    //***************************************

    // Circular border parameterizer (the default)
    typedef CGAL::Circular_border_arc_length_parameterizer_3<Parameterization_polyhedron_adaptor>
                                                        Border_parameterizer;
    // TAUCS solver
    typedef CGAL::Taucs_solver_traits<double>           Solver;

    // Floater Mean Value Coordinates parameterizer (circular border)
    // with TAUCS solver
    typedef CGAL::Mean_value_coordinates_parameterizer_3<Parameterization_polyhedron_adaptor,
                                                         Border_parameterizer,
                                                         Solver>
                                                        Parameterizer;

    Parameterizer::Error_code err = CGAL::parameterize(&mesh_adaptor, Parameterizer());
    if (err != Parameterizer::OK)
        std::cerr << "FATAL ERROR: " << Parameterizer::get_error_message(err) << std::endl;

    //***************************************
    // Output
    //***************************************

    // Write Postscript file
    if (err == Parameterizer::OK)
    {
        if ( ! write_file_eps(mesh_adaptor, output_filename) )
        {
            std::cerr << "FATAL ERROR: cannot write file " << output_filename << std::endl;
            return EXIT_FAILURE;
        }   
    }

    return (err == Parameterizer::OK) ? EXIT_SUCCESS : EXIT_FAILURE;
}


#else // CGAL_USE_TAUCS


#include <stdio.h>
#include <stdlib.h>

// ----------------------------------------------------------------------------
// Empty main() if TAUCS is not installed
// ----------------------------------------------------------------------------

int main(int argc,char * argv[])
{
    std::cerr << "Skip test as TAUCS is not installed" << std::endl;
    return EXIT_SUCCESS;
}


#endif // CGAL_USE_TAUCS

