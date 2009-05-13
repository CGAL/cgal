// APSS_reconstruction.cpp

//----------------------------------------------------------
// APSS reconstruction method:
// Read a point set or a mesh's set of vertices, reconstruct a surface,
// and save the surface.
// Output format is .off.
//----------------------------------------------------------
// APSS_reconstruction file_in file_out [options]

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Implicit_surface_3.h>

// This package
#include <CGAL/APSS_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/point_set_property_map.h>
#include <CGAL/IO/read_xyz_point_set.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/polyhedron_connected_components.h>

#include "compute_normal.h"

#include <deque>
#include <cstdlib>
#include <fstream>


// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef Kernel::Sphere_3 Sphere;
typedef std::deque<Point_with_normal> PointList;

// polyhedron
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

// APSS implicit function
typedef CGAL::APSS_reconstruction_function<Kernel> APSS_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Implicit_surface_3<Kernel, APSS_reconstruction_function> Surface_3;


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "APSS reconstruction method." << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 < 2)
    {
      std::cerr << "Read a point set or a mesh's set of vertices, reconstruct a surface,\n";
      std::cerr << "and save the surface.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file_in file_out [options]\n";
      std::cerr << "Input file formats are .off (mesh) and .xyz or .pwn (point set).\n";
      std::cerr << "Output file format is .off.\n";
      std::cerr << "Options:\n";
      std::cerr << "  -sm_radius <float>     Radius upper bound (default=0.1 * point set radius)\n";
      std::cerr << "  -sm_distance <float>   Distance upper bound (default=0.003 * point set radius)\n";
      std::cerr << "  -smooth <float>        Smoothness factor (default = 2)\n";
      std::cerr << "                           - should be greater than 1,\n";
      std::cerr << "                           - high numbers lead to smoother surfaces.\n";
      return EXIT_FAILURE;
    }

    // APSS options
    FT sm_angle = 20.0; // Min triangle angle (degrees). 20 = fast, 30 guaranties convergence.
    FT sm_radius = 0.1; // Max triangle radius w.r.t. point set radius. 0.1 is fine.
    FT sm_distance = 0.003; // Approximation error w.r.t. p.s.r.. For APSS: 0.015 = fast, 0.003 = smooth.
    FT smoothness = 2; // smoothness factor (the large the smoother and slower)

    // decode parameters
    std::string input_filename  = argv[1];
    std::string output_filename = argv[2];
    for (int i=3; i+1<argc ; ++i)
    {
      if (std::string(argv[i])=="-sm_radius")
        sm_radius = atof(argv[++i]);
      else if (std::string(argv[i])=="-sm_distance")
        sm_distance = atof(argv[++i]);
      else if (std::string(argv[i])=="-smooth")
        smoothness = atof(argv[++i]);
      else
        std::cerr << "invalid option " << argv[i] << "\n";
    }

    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Load mesh/point set
    //***************************************

    PointList points;

    // If OFF file format
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;
    std::string extension = input_filename.substr(input_filename.find_last_of('.'));
    if (extension == ".off" || extension == ".OFF")
    {
      // Read the mesh file in a polyhedron
      std::ifstream stream(input_filename.c_str());
      Polyhedron input_mesh;
      CGAL::scan_OFF(stream, input_mesh, true /* verbose */);
      if(!stream || !input_mesh.is_valid() || input_mesh.empty())
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
      }

      // Convert Polyhedron vertices to point set.
      // Compute vertices' normals from connectivity.
      Polyhedron::Vertex_const_iterator v;
      for (v = input_mesh.vertices_begin(); v != input_mesh.vertices_end(); v++)
      {
        const Point& p = v->point();
        Vector n = compute_vertex_normal<Polyhedron::Vertex,Kernel>(*v);
        points.push_back(Point_with_normal(p,n));
      }
    }
    // If XYZ file format
    // FIXME fails to compile:
    //  CGAL/point_set_property_map.h:39: error: forming reference to void                                                                                                                                                                                 
    //  CGAL/point_set_property_map.h:42: error: forming reference to void                                                                                                                                                                                 
    //  CGAL/point_set_property_map.h:47: error: forming reference to void
    // 
//     else if (extension == ".xyz" || extension == ".XYZ" ||
//              extension == ".pwn" || extension == ".PWN")
//     {
//       // Reads the point set file in points[].
//       // Note: read_xyz_point_set() requires an iterator over points
//       //       + property maps to access each point's position and normal.
//       std::ifstream stream(input_filename.c_str());
//       if(!stream ||
//          !CGAL::read_xyz_point_set(stream,
//                                    std::back_inserter(points),
//                                    CGAL::make_dereference_property_map(std::back_inserter(points)),
//                                    CGAL::make_normal_vector_property_map(points.begin())))
//       {
//         std::cerr << "Error: cannot read file " << input_filename << std::endl;
//         return EXIT_FAILURE;
//       }
//     }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    // Print status
    int nb_vertices = points.size();
    std::cerr << "Read file " << input_filename << ": " << nb_vertices << " vertices, "
                                                        << task_timer.time() << " seconds"
                                                        << std::endl;
    task_timer.reset();

    //***************************************
    // Check requirements
    //***************************************

    if (nb_vertices == 0)
    {
      std::cerr << "Error: empty file" << std::endl;
      return EXIT_FAILURE;
    }

    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      return EXIT_FAILURE;
    }

    //***************************************
    // Compute implicit function
    //***************************************

    std::cerr << "Compute APSS implicit function (smoothness=" << smoothness << ")...\n";

    // Create implicit function and insert vertices.
    // Note: APSS_implicit_function() requires an iterator over points
    //       + property maps to access each point's position and normal.
    //       The position property map has a default value and is omitted here.
    APSS_reconstruction_function implicit_function(points.begin(), points.end(),
                                                   CGAL::make_normal_vector_property_map(points.begin()),
                                                   smoothness);

    // Recover memory used by points[]
    points.clear();

    // Print status
    std::cerr << "Compute implicit function: " << task_timer.time() << " seconds\n";
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    // Get point inside the implicit surface
    Point inner_point = implicit_function.get_inner_point();
    FT inner_point_value = implicit_function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
      return EXIT_FAILURE;
    }

    // Get implicit function's radius
    Sphere bounding_sphere = implicit_function.bounding_sphere();
    FT size = sqrt(bounding_sphere.squared_radius());

    // defining the implicit surface = implicit function + bounding sphere centered at inner_point
    Point sm_sphere_center = inner_point;
    FT    sm_sphere_radius = size + std::sqrt(CGAL::squared_distance(bounding_sphere.center(),inner_point));
    sm_sphere_radius *= 1.01; // <= the Surface Mesher fails if the sphere does not contain the surface
    Surface_3 surface(implicit_function,
                      Sphere(sm_sphere_center,sm_sphere_radius*sm_sphere_radius));

    // defining meshing criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*size,  // Max triangle radius
                                                        sm_distance*size); // Approximation error

    // meshing surface
    STr tr; // 3D-Delaunay triangulation for Surface Mesher
    C2t3 surface_mesher_c2t3 (tr); // 2D-complex in 3D-Delaunay triangulation
    CGAL::make_surface_mesh(surface_mesher_c2t3, surface, criteria, CGAL::Manifold_with_boundary_tag());

    // Print status
    std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
                                     << tr.number_of_vertices() << " output vertices"
                                     << std::endl;
    task_timer.reset();

    if(tr.number_of_vertices() == 0)
      return EXIT_FAILURE;

    // Convert to polyhedron
    Polyhedron output_mesh;
    CGAL::output_surface_facets_to_polyhedron(surface_mesher_c2t3, output_mesh);

    //***************************************
    // Erase small connected components
    //***************************************

    std::cerr << "Erase small connected components...\n";

    unsigned int nb_erased_components =
      CGAL::erase_small_polyhedron_connected_components(output_mesh);

    // Print status
    std::cerr << "Erase small connected components: " << task_timer.time() << " seconds, "
                                                      << nb_erased_components << " components erased"
                                                      << std::endl;
    task_timer.reset();

    //***************************************
    // save the mesh
    //***************************************

    std::cerr << "Write file " << output_filename << std::endl << std::endl;

    std::ofstream out(output_filename.c_str());
    out << output_mesh;

    return EXIT_SUCCESS;
}

