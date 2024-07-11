#include <CGAL/Installation/internal/disable_deprecation_warnings_and_errors.h>

//----------------------------------------------------------
// Compares Poisson using Mesh_3
// VS Poisson using Surface_mesher and Poisson_implicit_surface_3
// see issue https://github.com/CGAL/cgal/issues/8266
//----------------------------------------------------------

// CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/facets_in_complex_3_to_triangle_mesh.h>

#include <CGAL/Surface_mesh_default_triangulation_3.h>
#include <CGAL/make_surface_mesh.h>
#include <CGAL/IO/facets_in_complex_2_to_triangle_mesh.h>

#include <CGAL/Poisson_implicit_surface_3.h>
#include <CGAL/Poisson_reconstruction_function.h>
#include <CGAL/IO/read_points.h>
#include <CGAL/compute_average_spacing.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>

#include <CGAL/Timer.h>

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
typedef CGAL::Poisson_implicit_surface_3<Kernel, Poisson_reconstruction_function> Surface_3;

// Surface mesher
typedef CGAL::Surface_mesh_default_triangulation_3 STr;
typedef CGAL::Surface_mesh_complex_2_in_triangulation_3<STr> C2t3;

// Mesh_3
typedef CGAL::Labeled_mesh_domain_3<Kernel> Mesh_domain;
typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;
typedef CGAL::Mesh_criteria_3<Tr> Mesh_criteria;


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
  CGAL::get_default_random() = CGAL::Random(0);
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
      std::cerr << "  -frac <float>          factor applied to sm_radius (default = 1.)\n";
      std::cerr << "Running " << argv[0] << "data/kitten.xyz kitten_poisson-20-100-0.5.off -sm_distance 0.5\n";
    }

    // Poisson options
    FT sm_angle = 20.0; // Min triangle angle (degrees).
    FT sm_radius = 100; // Max triangle size w.r.t. point set average spacing.
    FT sm_distance = 0.25; // Approximation error w.r.t. point set average spacing.
    std::string solver_name = "eigen"; // Sparse linear solver name.
    double approximation_ratio = 0.02;
    double average_spacing_ratio = 5;
    double frac = 1.;

    // decode parameters
    std::string input_filename  = (argc > 1) ? argv[1] : CGAL::data_file_path("points_3/kitten.xyz");
    std::string output_filename = (argc > 2) ? argv[2] : "kitten_poisson-20-100-0.5.off";
    for (int i=3; i+1<argc ; ++i)
    {
      if (std::string(argv[i]) == "-sm_radius")
        sm_radius = atof(argv[++i]);
      else if (std::string(argv[i]) == "-sm_distance")
        sm_distance = atof(argv[++i]);
      else if (std::string(argv[i]) == "-solver")
        solver_name = argv[++i];
      else if (std::string(argv[i]) == "-approx")
        approximation_ratio = atof(argv[++i]);
      else if (std::string(argv[i]) == "-ratio")
        average_spacing_ratio = atof(argv[++i]);
      else if (std::string(argv[i]) == "-frac")
        frac = atof(argv[++i]);
      else {
        std::cerr << "Error: invalid option " << argv[i] << "\n";
        return EXIT_FAILURE;
      }
    }

    if (argc == 1) sm_distance = 0.5;

    const std::size_t last_dot = output_filename.find_last_of(".");
    const std::string output_extension = output_filename.substr(last_dot);
    const std::string output_basename = output_filename.substr(0, last_dot);

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
             extension == ".pwn" || extension == ".PWN" ||
             extension == ".ply" || extension == ".PLY")
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

    std::cerr << "\nComputes Poisson implicit function...\n";

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
    std::cerr << std::endl << std::endl;

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

    // Meshing criteria
    const double fangle = sm_angle;
    const double fsize = frac * sm_radius * average_spacing;
    const double fdist = sm_distance * average_spacing;

    const double implicit_function_time = reconstruction_timer.time();
    reconstruction_timer.reset();

    // MESH_3
    {
      CGAL::Real_timer meshing_timer;
      meshing_timer.start();

      std::cout << "* Use Mesh_3 *" << std::endl;
      // Defines generation criteria
      Mesh_criteria criteria(CGAL::parameters::facet_angle = fangle,
                             CGAL::parameters::facet_size = fsize,
                             CGAL::parameters::facet_distance = fdist);

      // Defines mesh domain
      Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(function, bsphere,
          CGAL::parameters::relative_error_bound(sm_dichotomy_error / sm_sphere_radius));

      // Generates mesh with manifold option
      C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain, criteria,
                                          CGAL::parameters::no_exude().no_perturb()
                                          .manifold_with_boundary());
      meshing_timer.stop();

      const Tr& tr = c3t3.triangulation();
      // Prints status
      std::cerr << "Mesh_3 meshing: " << meshing_timer.time() << " seconds, "
                                      << tr.number_of_vertices() << " output vertices"
                                      << std::endl;

      if (tr.number_of_vertices() == 0)
        return EXIT_FAILURE;

      // Prints total reconstruction duration
      reconstruction_timer.stop();
      std::cerr << "Total reconstruction (implicit function + meshing): "
        << (implicit_function_time + reconstruction_timer.time()) << " seconds\n";
      reconstruction_timer.reset();

      // Converts to polyhedron
      Polyhedron output_mesh;
      CGAL::facets_in_complex_3_to_triangle_mesh(c3t3, output_mesh);

      std::ofstream out(output_basename + "_mesh_3.off");
      out << output_mesh;
      out.close();
    }

    // SURFACE_MESHER
    {
      CGAL::Real_timer meshing_timer;
      meshing_timer.start();
      reconstruction_timer.start();

      std::cout << "\n\n* Use Surface_mesher with Poisson_implicit_surface_3 *" << std::endl;
      Surface_3 surface(function,
                        Sphere(inner_point, sm_sphere_radius * sm_sphere_radius),
                        sm_dichotomy_error / sm_sphere_radius);

      // Defines surface mesh generation criteria
      CGAL::Surface_mesh_default_criteria_3<STr> criteria(fangle, fsize, fdist);

      // Generates surface mesh with manifold option
      STr tr; // 3D Delaunay triangulation for surface mesh generation
      C2t3 c2t3(tr); // 2D complex in 3D Delaunay triangulation
      CGAL::make_surface_mesh(c2t3,                                 // reconstructed mesh
                              surface,                              // implicit surface
                              criteria,                             // meshing criteria
                              CGAL::Manifold_with_boundary_tag());  // require manifold mesh
      meshing_timer.stop();

      // Prints status
      std::cerr << "Surface meshing: " << meshing_timer.time() << " seconds, "
                                       << tr.number_of_vertices() << " output vertices"
                                       << std::endl;

      if (tr.number_of_vertices() == 0)
        return EXIT_FAILURE;

      // Prints total reconstruction duration
      reconstruction_timer.stop();
      std::cerr << "Total reconstruction (implicit function + meshing): "
        << (implicit_function_time + reconstruction_timer.time()) << " seconds\n";

      Polyhedron output_mesh;
      CGAL::facets_in_complex_2_to_triangle_mesh(c2t3, output_mesh);

      std::ofstream out(output_basename + "_surface_mesher.off");
      out << output_mesh;
      out.close();
    }

    return EXIT_SUCCESS;
}
