// poisson_reconstruction.cpp

//----------------------------------------------------------
// Poisson Delaunay Reconstruction method.
// Reads a point set or a mesh's set of vertices, reconstructs a surface using Poisson,
// and saves the surface.
// Output format is .off.
//----------------------------------------------------------
// poisson_reconstruction file_in file_out [options]

// CGAL
#include <CGAL/AABB_tree.h> // must be included before kernel
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Timer.h>
#include <CGAL/trace.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Poisson_implicit_surface_3.h>
#include <CGAL/IO/output_surface_facets_to_polyhedron.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/Point_with_normal_3.h>
#include <CGAL/IO/read_xyz_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <deque>
#include <cstdlib>
#include <fstream>
#include <math.h>
#include <boost/foreach.hpp>

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

// Poisson implicit function
typedef CGAL::Poisson_reconstruction_function<Kernel> Poisson_reconstruction_function;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;
typedef CGAL::Poisson_implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

// AABB tree
typedef CGAL::AABB_face_graph_triangle_primitive<Polyhedron> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> AABB_traits;
typedef CGAL::AABB_tree<AABB_traits> AABB_tree;


struct Counter {
  int i, N;
  Counter(int N)
    : i(0), N(N)
  {}

  void operator()()
  {
    i++;
    if(i == N){
      std::cerr << "Counter reached " << N << std::endl;
    }
  }
  
};

struct InsertVisitor {

  Counter& c;
  InsertVisitor(Counter& c)
    : c(c)
  {}

  void before_insertion()
  {
    c();
  }

};


// ----------------------------------------------------------------------------
// main()
// ----------------------------------------------------------------------------

int main(int argc, char * argv[])
{
    std::cerr << "Poisson Delaunay Reconstruction method" << std::endl;

    //***************************************
    // decode parameters
    //***************************************

    // usage
    if (argc-1 < 2)
    {
      std::cerr << "Reads a point set or a mesh's set of vertices, reconstructs a surface using Poisson,\n";
      std::cerr << "and saves the surface.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " file_in file_out [options]\n";
      std::cerr << "Input file formats are .off (mesh) and .xyz or .pwn (point set).\n";
      std::cerr << "Output file format is .off.\n";
      std::cerr << "Options:\n";
      std::cerr << "  -sm_radius <float>     Radius upper bound (default=100 * average spacing)\n";
      std::cerr << "  -sm_distance <float>   Distance upper bound (default=0.25 * average spacing)\n";
      
      return EXIT_FAILURE;
    }

    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle (degrees).
    FT sm_radius = 100; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.25; // Approximation error w.r.t. point set average spacing.
    std::string solver_name = "eigen"; // Sparse linear solver name.
    double approximation_ratio = 0.02;
    double average_spacing_ratio = 5;

    // decode parameters
    std::string input_filename  = argv[1];
    std::string output_filename = argv[2];
    for (int i=3; i+1<argc ; ++i)
    {
      if (std::string(argv[i])=="-sm_radius")
        sm_radius = atof(argv[++i]);
      else if (std::string(argv[i])=="-sm_distance")
        sm_distance = atof(argv[++i]);
      else if (std::string(argv[i])=="-solver")
        solver_name = argv[++i];
      else if (std::string(argv[i])=="-approx")
        approximation_ratio = atof(argv[++i]);
      else if (std::string(argv[i])=="-ratio")
        average_spacing_ratio = atof(argv[++i]);
      else {
        std::cerr << "Error: invalid option " << argv[i] << "\n";
        return EXIT_FAILURE;
      }
    }

    CGAL::Timer task_timer; task_timer.start();

    //***************************************
    // Loads mesh/point set
    //***************************************

    PointList points;

    // If OFF file format
    std::cerr << "Open " << input_filename << " for reading..." << std::endl;
    std::string extension = input_filename.substr(input_filename.find_last_of('.'));
    if (extension == ".off" || extension == ".OFF")
    {
      // Reads the mesh file in a polyhedron
      std::ifstream stream(input_filename.c_str());
      Polyhedron input_mesh;
      CGAL::scan_OFF(stream, input_mesh, true /* verbose */);
      if(!stream || !input_mesh.is_valid() || input_mesh.empty())
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
      }

      // Converts Polyhedron vertices to point set.
      // Computes vertices normal from connectivity.
      BOOST_FOREACH(boost::graph_traits<Polyhedron>::vertex_descriptor v,
                    vertices(input_mesh)){
        const Point& p = v->point();
        Vector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(v,input_mesh);
        points.push_back(Point_with_normal(p,n));
      }
    }
    // If XYZ file format
    else if (extension == ".xyz" || extension == ".XYZ" ||
             extension == ".pwn" || extension == ".PWN")
    {
      // Reads the point set file in points[].
      // Note: read_xyz_points_and_normals() requires an iterator over points
      // + property maps to access each point's position and normal.
      // The position property map can be omitted here as we use iterators over Point_3 elements.
      std::ifstream stream(input_filename.c_str());
      if (!stream ||
          !CGAL::read_xyz_points_and_normals(
                                stream,
                                std::back_inserter(points),
                                CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type())))
      {
        std::cerr << "Error: cannot read file " << input_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    // Prints status
    int nb_points = points.size();
    std::cerr << "Reads file " << input_filename << ": " << nb_points << " points, "
                                                        << task_timer.time() << " seconds"
                                                        << std::endl;
    task_timer.reset();

    //***************************************
    // Checks requirements
    //***************************************

    if (nb_points == 0)
    {
      std::cerr << "Error: empty point set" << std::endl;
      return EXIT_FAILURE;
    }

    bool points_have_normals = (points.begin()->normal() != CGAL::NULL_VECTOR);
    if ( ! points_have_normals )
    {
      std::cerr << "Input point set not supported: this reconstruction method requires oriented normals" << std::endl;
      return EXIT_FAILURE;
    }

    CGAL::Timer reconstruction_timer; reconstruction_timer.start();

    
    Counter counter(std::distance(points.begin(), points.end()));
    InsertVisitor visitor(counter) ;
    

    //***************************************
    // Computes implicit function
    //***************************************

    std::cerr << "Computes Poisson implicit function...\n";

    // Creates implicit function from the read points.
    // Note: this method requires an iterator over points
    // + property maps to access each point's position and normal.
    // The position property map can be omitted here as we use iterators over Point_3 elements.
    Poisson_reconstruction_function function(
                              points.begin(), points.end(),
                              CGAL::make_identity_property_map(PointList::value_type()),
                              CGAL::make_normal_of_point_with_normal_pmap(PointList::value_type()),
                              visitor);

    #ifdef CGAL_EIGEN3_ENABLED
    {
      if (solver_name == "eigen")
      {
        std::cerr << "Use Eigen 3\n";
        CGAL::Eigen_solver_traits<Eigen::ConjugateGradient<CGAL::Eigen_sparse_symmetric_matrix<double>::EigenType> > solver;
        if ( ! function.compute_implicit_function(solver, visitor, 
                                                approximation_ratio,
                                                average_spacing_ratio) )
        {
          std::cerr << "Error: cannot compute implicit function" << std::endl;
          return EXIT_FAILURE;
        }
      }    
      else
      {
        std::cerr << "Error: invalid solver " << solver_name << "\n";
        return EXIT_FAILURE;
      }
    }
    #else
    {
      std::cerr << "Error: invalid solver " << solver_name << "\n";
      return EXIT_FAILURE;
    }
    #endif


    // Prints status
    std::cerr << "Total implicit function (triangulation+refinement+solver): " << task_timer.time() << " seconds\n";
    task_timer.reset();

    //***************************************
    // Surface mesh generation
    //***************************************

    std::cerr << "Surface meshing...\n";

    // Computes average spacing
    FT average_spacing = CGAL::compute_average_spacing(points.begin(), points.end(),
                                                       6 /* knn = 1 ring */);

    // Gets one point inside the implicit surface
    Point inner_point = function.get_inner_point();
    FT inner_point_value = function(inner_point);
    if(inner_point_value >= 0.0)
    {
      std::cerr << "Error: unable to seed (" << inner_point_value << " at inner_point)" << std::endl;
      return EXIT_FAILURE;
    }

    // Gets implicit function's radius
    Sphere bsphere = function.bounding_sphere();
    FT radius = std::sqrt(bsphere.squared_radius());

    // Defines the implicit surface: requires defining a
  	// conservative bounding sphere centered at inner point.
    FT sm_sphere_radius = 5.0 * radius;
    FT sm_dichotomy_error = sm_distance*average_spacing/1000.0; // Dichotomy error must be << sm_distance
    Surface_3 surface(function,
                      Sphere(inner_point,sm_sphere_radius*sm_sphere_radius),
                      sm_dichotomy_error/sm_sphere_radius);

    // Defines surface mesh generation criteria
    CGAL::Surface_mesh_default_criteria_3<STr> criteria(sm_angle,  // Min triangle angle (degrees)
                                                        sm_radius*average_spacing,  // Max triangle size
                                                        sm_distance*average_spacing); // Approximation error

    CGAL_TRACE_STREAM << "  make_surface_mesh(sphere center=("<<inner_point << "),\n"
                      << "                    sphere radius="<<sm_sphere_radius<<",\n"
                      << "                    angle="<<sm_angle << " degrees,\n"
                      << "                    triangle size="<<sm_radius<<" * average spacing="<<sm_radius*average_spacing<<",\n"
                      << "                    distance="<<sm_distance<<" * average spacing="<<sm_distance*average_spacing<<",\n"
                      << "                    dichotomy error=distance/"<<sm_distance*average_spacing/sm_dichotomy_error<<",\n"
                      << "                    Manifold_with_boundary_tag)\n";

    // Generates surface mesh with manifold option
    STr tr; // 3D Delaunay triangulation for surface mesh generation
    C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
    CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                            surface,                              // implicit surface
                            criteria,                             // meshing criteria
                            CGAL::Manifold_with_boundary_tag());  // require manifold mesh

    // Prints status
    std::cerr << "Surface meshing: " << task_timer.time() << " seconds, "
                                     << tr.number_of_vertices() << " output vertices"
                                     << std::endl;
    task_timer.reset();

    if(tr.number_of_vertices() == 0)
      return EXIT_FAILURE;

    // Converts to polyhedron
    Polyhedron output_mesh;
    CGAL::output_surface_facets_to_polyhedron(c2t3, output_mesh);

    // Prints total reconstruction duration
    std::cerr << "Total reconstruction (implicit function + meshing): " << reconstruction_timer.time() << " seconds\n";

    //***************************************
    // Computes reconstruction error
    //***************************************

    // Constructs AABB tree and computes internal KD-tree
    // data structure to accelerate distance queries
    AABB_tree tree(faces(output_mesh).first, faces(output_mesh).second, output_mesh);
    tree.accelerate_distance_queries();

    // Computes distance from each input point to reconstructed mesh
    double max_distance = DBL_MIN;
    double avg_distance = 0;
    for (PointList::const_iterator p=points.begin(); p!=points.end(); p++)
    {
      double distance = std::sqrt(tree.squared_distance(*p));

      max_distance = (std::max)(max_distance, distance);
      avg_distance += distance;
    }
    avg_distance /= double(points.size());

    std::cerr << "Reconstruction error:\n"
              << "  max = " << max_distance << " = " << max_distance/average_spacing << " * average spacing\n"
              << "  avg = " << avg_distance << " = " << avg_distance/average_spacing << " * average spacing\n";

    //***************************************
    // Saves reconstructed surface mesh
    //***************************************

    std::cerr << "Write file " << output_filename << std::endl << std::endl;
    std::ofstream out(output_filename.c_str());
    out << output_mesh;

    return EXIT_SUCCESS;
}
