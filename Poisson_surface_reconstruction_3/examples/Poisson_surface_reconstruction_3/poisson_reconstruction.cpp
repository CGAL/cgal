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
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/Poisson_implicit_surface_3.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <deque>
#include <cstdlib>
#include <fstream>
#include <math.h>

// ----------------------------------------------------------------------------
// Types
// ----------------------------------------------------------------------------

// kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

// Simple geometric types
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef std::pair<Point, Vector> Point_with_normal;
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
  std::size_t i, N;
  Counter(std::size_t N)
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
    if(argc == 1)
    {
      std::cerr << "Reads a point set or a mesh's set of vertices, reconstructs a surface using Poisson,\n";
      std::cerr << "and saves the surface.\n";
      std::cerr << "\n";
      std::cerr << "Usage: " << argv[0] << " [file_in] [file_out] [options]\n";
      std::cerr << "Input file formats are .off (mesh) and .xyz or .pwn (point set).\n";
      std::cerr << "Output file format is .off.\n";
      std::cerr << "Options:\n";
      std::cerr << "  -sm_radius <float>     Radius upper bound (default=100 * average spacing)\n";
      std::cerr << "  -sm_distance <float>   Distance upper bound (default=0.25 * average spacing)\n";
      std::cerr << "Running " << argv[0] << "data/kitten.xyz kitten_poisson-20-100-0.5.off -sm_distance 0.5\n";
    }

    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle (degrees).
    FT sm_radius = 100; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.25; // Approximation error w.r.t. point set average spacing.
    std::string solver_name = "eigen"; // Sparse linear solver name.
    double approximation_ratio = 0.02;
    double average_spacing_ratio = 5;

    // decode parameters
    std::string input_filename  = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/kitten.xyz");
    std::string output_filename = (argc > 2) ? argv[2] : "kitten_poisson-20-100-0.5.off";
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

    if (argc == 1) sm_distance = 0.5;

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
      for(boost::graph_traits<Polyhedron>::vertex_descriptor v :
                    vertices(input_mesh)){
        const Point& p = v->point();
        Vector n = CGAL::Polygon_mesh_processing::compute_vertex_normal(v,input_mesh);
        points.push_back(std::make_pair(p,n));
      }
    }
    // If XYZ file format
    else if (extension == ".xyz" || extension == ".XYZ" ||
             extension == ".pwn" || extension == ".PWN")
    {
      // Reads the point set file in points[].
      // Note: read_points() requires an iterator over points
      // + property maps to access each point's position and normal.
      if (!CGAL::IO::read_points(input_filename.c_str(), std::back_inserter(points),
                                  CGAL::parameters::point_map(CGAL::make_first_of_pair_property_map(Point_with_normal()))
                                                    .normal_map(CGAL::make_second_of_pair_property_map(Point_with_normal()))))
      {
        std::cerr << "Error: cannot read input file!" << input_filename << std::endl;
        return EXIT_FAILURE;
      }
    }
    else
    {
      std::cerr << "Error: cannot read file " << input_filename << std::endl;
      return EXIT_FAILURE;
    }

    // Prints status
    std::size_t nb_points = points.size();
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

    bool points_have_normals = (points.begin()->second != CGAL::NULL_VECTOR);
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
    Poisson_reconstruction_function function(
                              points.begin(), points.end(),
                              CGAL::make_first_of_pair_property_map(Point_with_normal()),
                              CGAL::make_second_of_pair_property_map(Point_with_normal()),
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
    FT average_spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>
      (points, 6 /* knn = 1 ring */,
       CGAL::parameters::point_map (CGAL::make_first_of_pair_property_map(Point_with_normal())));

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

    std::cerr         << "  make_surface_mesh(sphere center=("<<inner_point << "),\n"
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
    CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, output_mesh);

    // Prints total reconstruction duration
    std::cerr << "Total reconstruction (implicit function + meshing): " << reconstruction_timer.time() << " seconds\n";

    //***************************************
    // Computes reconstruction error
    //***************************************

    // Constructs AABB tree and computes internal KD-tree
    // data structure to accelerate distance queries
    AABB_tree tree(faces(output_mesh).first, faces(output_mesh).second, output_mesh);

    // Computes distance from each input point to reconstructed mesh
    double max_distance = DBL_MIN;
    double avg_distance = 0;
    for (PointList::const_iterator p=points.begin(); p!=points.end(); p++)
    {
      double distance = std::sqrt(tree.squared_distance(p->first));

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
