#include "benchmark_config.h"
#include "benchmark_xml.h"
#include "mesh_quality.h"

std::string XML_perf_data::default_filename;

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Mesh_triangulation_3.h>
#include <CGAL/Mesh_complex_3_in_triangulation_3.h>
#include <CGAL/Mesh_criteria_3.h>
#include <CGAL/Labeled_mesh_domain_3.h>
#include <CGAL/Mesh_polyhedron_3.h>
#include <CGAL/Polyhedral_mesh_domain_3.h>
#include <CGAL/Polyhedral_mesh_domain_with_features_3.h>
#include <CGAL/Mesh_domain_with_polyline_features_3.h>
#include <CGAL/make_mesh_3.h>
#include <CGAL/refine_mesh_3.h>

#include <CGAL/boost/graph/IO/polygon_mesh_io.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_arena.h>
# define TBB_PREVIEW_GLOBAL_CONTROL 1
# include <tbb/global_control.h>
#endif


#include <filesystem>
#include <fstream>
#include <iostream>
#include <signal.h>
#include <string>
#include <sstream>

namespace PMP = CGAL::Polygon_mesh_processing;

// basic types from kernel
typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::FT FT;
typedef K::Point_3 Point;
typedef K::Sphere_3 Sphere;

// To avoid verbose function and named parameters call
using namespace CGAL::parameters;

#include "implicit_functions.h"

std::string get_output_filename(const std::string& input_name)
{
  std::string filename = std::string(input_name);
  filename = filename.substr(filename.find_last_of("/") + 1, filename.length() - 1);
  filename = filename.substr(0, filename.find_last_of("."));
  return filename;
}

std::string get_technique()
{
  std::string tech;
#ifdef CGAL_CONCURRENT_MESH_3

  tech += "Task-scheduler (auto";
# ifdef CGAL_MESH_3_LOAD_BASED_WORKSHARING
    tech += ", load-based worksharing";
#endif
  tech += ")";

#else // !CGAL_CONCURRENT_MESH_3

  tech += "Sequential ";
# if defined(CGAL_MESH_3_USE_LAZY_UNSORTED_REFINEMENT_QUEUE)
#   ifdef CGAL_MESH_3_IF_UNSORTED_QUEUE_JUST_SORT_AFTER_SCAN
  tech += "(sort after scan only";
#   else
  tech += "(unsorted";
#   endif
# elif defined(CGAL_MESH_3_USE_LAZY_SORTED_REFINEMENT_QUEUE)
  tech += "(sorted";
# else
  tech += "(NOT LAZY, sorted";
# endif

#ifdef CGAL_SEQUENTIAL_MESH_3_ADD_OUTSIDE_POINTS_ON_A_FAR_SPHERE
  tech += ", points on far sphere)";
#else
  tech += ")";
#endif

#endif // CGAL_CONCURRENT_MESH_3

  return tech;
}

void display_info(int num_threads)
{
  std::cout << get_technique() << std::endl;

#ifdef CGAL_CONCURRENT_MESH_3

  if(num_threads != -1)
    std::cout << "Num threads = " << num_threads << std::endl;
  else
    std::cout << "Num threads = AUTO" << std::endl;

#else // !CGAL_CONCURRENT_MESH_3

# ifdef CGAL_MESH_3_INITIAL_POINTS_NO_RANDOM_SHOOTING
  std::cout << "NO random shooting)" << std::endl;
# else
  std::cout << "WITH random shooting)" << std::endl;
# endif

#endif // CGAL_CONCURRENT_MESH_3
}

void xml_perf_set_technique()
{
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Technique", get_technique());
}

#ifdef CGAL_MESH_3_IMPLICIT_WITH_FEATURES
// To add a crease (feature) to some implicit function
typedef std::vector<Point> Crease;
typedef std::list<Crease> Creases;

void add_crease(const Point& a,
                const Point& b,
                Creases& creases)
{
  Crease crease;
  crease.push_back(a);
  crease.push_back(b);
  creases.push_back(crease);
}
#endif

enum Exit_code
{
  // Success
  VALID_OUTPUT = 0,

  // Failure
  INPUT_IS_INVALID = 1,
  OUTPUT_IS_INVALID = 2
};

// the call to Mesh_3 happens here
template <typename C3t3, typename Domain>
Exit_code make_mesh(const Domain& domain,
                    const CGAL::Mesh_criteria_3<typename C3t3::Triangulation>& criteria,
                    const std::string& output_filename)
{
#ifdef _DEBUG
  double timelimit = 10;
  double sliverbound = 2;
#else
  double timelimit = 0; // when making charts, we run this executable a timeout
  double sliverbound = 2;
#endif

  CGAL::Real_timer t;
  t.start();

  C3t3 c3t3 = CGAL::make_mesh_3<C3t3>(domain
                                   , criteria
# ifdef CGAL_MESH_3_BENCHMARK_LLOYD
                                   , lloyd(time_limit=timelimit)
# else
                                   , no_lloyd()
# endif
# ifdef CGAL_MESH_3_BENCHMARK_ODT
                                   , odt(time_limit=timelimit)
# else
                                   , no_odt()
#endif
# ifdef CGAL_MESH_3_BENCHMARK_PERTURB
                                   , perturb(time_limit = timelimit,
                                             sliver_bound = sliverbound)
# else
                                   , no_perturb()
#endif
#ifdef CGAL_MESH_3_BENCHMARK_EXUDE
                                   , exude(time_limit = timelimit,
                                           sliver_bound = sliverbound)
#else
                                   , no_exude()
#endif
#ifdef CGAL_MESH_3_MANIFOLD
                                   , manifold()
#else
                                   , non_manifold()
#endif
                                    );
  t.stop();

  CGAL_MESH_3_SET_PERFORMANCE_DATA("V", c3t3.triangulation().number_of_vertices());
  CGAL_MESH_3_SET_PERFORMANCE_DATA("F", c3t3.number_of_facets_in_complex());
  CGAL_MESH_3_SET_PERFORMANCE_DATA("C", c3t3.number_of_cells_in_complex());
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Mem", CGAL::Memory_sizer().virtual_size() >> 20);
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Total_time", t.time());

#ifdef CGAL_MESH_3_BENCHMARK_EXPORT_TO_MAYA
  std::cout << "Exporting to " << output_filename + ".maya (Maya)... ";
  std::ofstream out_maya(output_filename + ".maya");
  c3t3.output_to_maya(out_maya, true);
  std::cout << "done." << std::endl;
#endif

#ifdef CGAL_MESH_3_BENCHMARK_EXPORT_TO_MESH
  std::cout << "Exporting to " << output_filename + ".mesh (Medit)... ";
  // std::cout << "std::filesystem::current_path() = " << std::filesystem::current_path() << std::endl;
  std::ofstream out_medit(output_filename + ".mesh");
  c3t3.output_to_medit(out_medit, true);
  std::cout << "done." << std::endl;
#endif

  if(!c3t3.triangulation().is_valid() || !c3t3.is_valid())
  {
    std::cerr << "Error: invalid output" << std::endl;
    return OUTPUT_IS_INVALID;
  }

  if(c3t3.number_of_facets_in_complex() == 0)
  {
    std::cerr << "Error: no facets in output" << std::endl;
    return OUTPUT_IS_INVALID;
  }

  if(c3t3.number_of_cells_in_complex() == 0)
  {
    std::cerr << "Error: no cells in output" << std::endl;
    return OUTPUT_IS_INVALID;
  }

  generate_quality_metrics(c3t3);

  return VALID_OUTPUT;
}

template <typename C3t3, typename Domain>
Exit_code make_mesh(const Domain& domain,
                    double facet_sizing, double facet_approx, double facet_ang,
                    double cell_sizing, double cell_shape,
                    const std::string& output_filename)
{
  typedef CGAL::Mesh_criteria_3<typename C3t3::Triangulation> Mesh_criteria;

  const CGAL::Bbox_3 bbox = domain.bbox();

  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  if(USE_RELATIVE_CRITERIA_VALUES)
  {
    facet_sizing = diag_length / facet_sizing;
    facet_approx = diag_length / facet_approx;
    cell_sizing = diag_length / cell_sizing;
  }

  Mesh_criteria criteria(edge_size = facet_sizing,
                         facet_angle = facet_ang,
                         facet_size = facet_sizing,
                         facet_distance = facet_approx,
                         cell_size = cell_sizing,
                         cell_radius_edge_ratio = cell_shape);

  std::cout << " * edge max size: " << facet_sizing << std::endl
            << " * facet max size: " << facet_sizing << std::endl
            << " * facet approx error: " << facet_approx << std::endl
            << " * facet min angle: " << facet_ang << std::endl
            << " * cell max size: " << cell_sizing << std::endl
            << " * cell shape (radius-edge): " << cell_shape << std::endl;

  return make_mesh<C3t3>(domain, criteria, output_filename);
}

Exit_code make_mesh_polyhedron(const std::string& input_filename,
                               double facet_sizing, double facet_approx, double facet_angle,
                               double cell_sizing, double cell_shape)
{
  std::cout << "make_mesh_polyhedron(" << input_filename  << ")" << std::endl;

  // Domain
#ifdef CGAL_MESH_3_POLYHEDRON_WITH_FEATURES
  typedef CGAL::Mesh_polyhedron_3<K>::type Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_with_features_3<K> Mesh_domain;
#else
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::Polyhedral_mesh_domain_3<Polyhedron, K> Mesh_domain;
#endif

  // Triangulation
#ifdef CGAL_CONCURRENT_MESH_3
  typedef CGAL::Mesh_triangulation_3<Mesh_domain,
                                     CGAL::Kernel_traits<Mesh_domain>::Kernel,
                                     CGAL::Parallel_tag>::type Tr;
#else
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif

  // C3t3
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Create input polyhedron
  Polyhedron polyhedron;
  if(!CGAL::IO::read_polygon_mesh(input_filename, polyhedron))
  {
    std::cerr << "Error: Could not read '" << input_filename << "'" << std::endl;
    return INPUT_IS_INVALID;
  }

  if(is_empty(polyhedron) ||
     !is_triangle_mesh(polyhedron) ||
     !is_closed(polyhedron) ||
     has_degenerate_faces(polyhedron) ||
     PMP::does_self_intersect(polyhedron))
  {
    std::cerr << "Error: input has defects" << std::endl;
    return INPUT_IS_INVALID;
  }

  std::cout << "Building AABB tree... " << std::endl;
  // Create domain
  Mesh_domain domain(polyhedron);
  std::cout << "done." << std::endl;

#ifdef CGAL_MESH_3_POLYHEDRON_WITH_FEATURES
  std::cout << "Detecting features..." << std::endl;
  domain.detect_features();
  std::cout << "done." << std::endl;
#endif

  return make_mesh<C3t3>(domain,
                         facet_sizing, facet_approx, facet_angle,
                         cell_sizing, cell_shape,
                         get_output_filename(input_filename));
}

Exit_code make_mesh_3D_images(const std::string& input_filename,
                              double facet_sizing, double facet_approx, double facet_angle,
                              double cell_sizing, double cell_shape)
{
  std::cout << "make_mesh_3D_images(" << input_filename  << ")" << std::endl;

  // Domain
  typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;

  // Triangulation
#ifdef CGAL_CONCURRENT_MESH_3
  typedef CGAL::Mesh_triangulation_3<Mesh_domain,
                                     CGAL::Kernel_traits<Mesh_domain>::Kernel,
                                     CGAL::Parallel_tag>::type Tr;
#else
  typedef CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif

  // C3t3
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Load image
  CGAL::Image_3 image;
  image.read(input_filename.c_str());

  // Create domain
  Mesh_domain domain = Mesh_domain::create_labeled_image_mesh_domain(image);
  std::cout << "done." << std::endl;

  return make_mesh<C3t3>(domain,
                         facet_sizing, facet_approx, facet_angle,
                         cell_sizing, cell_shape,
                         get_output_filename(input_filename));
}

template <class ImplicitFunction>
Exit_code make_mesh_implicit(const std::string& function_name,
                             ImplicitFunction func,
                             double facet_sizing, double facet_approx, double facet_angle,
                             double cell_sizing, double cell_shape)
{
  std::cout << "make_mesh_implicit(" << function_name  << ")" << std::endl;

  // Domain
#ifdef CGAL_MESH_3_IMPLICIT_WITH_FEATURES
  typedef CGAL::Labeled_mesh_domain_3<K> Implicit_domain;
  typedef CGAL::Mesh_domain_with_polyline_features_3<Implicit_domain> Mesh_domain;
#else
  typedef CGAL::Labeled_mesh_domain_3<K> Mesh_domain;
#endif

  // Triangulation
#ifdef CGAL_CONCURRENT_MESH_3
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain,
                                              typename CGAL::Kernel_traits<Mesh_domain>::Kernel,
                                              CGAL::Parallel_tag>::type Tr;
#else
  typedef typename CGAL::Mesh_triangulation_3<Mesh_domain>::type Tr;
#endif

  // C3t3
  typedef CGAL::Mesh_complex_3_in_triangulation_3<Tr> C3t3;

  // Create domain
  Sphere bounding_sphere(CGAL::ORIGIN, 10.0 * 10.0);

  namespace p = CGAL::parameters;
  Mesh_domain domain = Mesh_domain::create_implicit_mesh_domain(p::function = func,
                                                                p::bounding_object = bounding_sphere
                                                                /*, p::relative_error_bound = 1e-7*/);

#ifdef CGAL_MESH_3_IMPLICIT_WITH_FEATURES
  // Add 12 feature creases
  Creases creases;
  Point p1(-1.0, -1.0, -1.0);
  Point p2(-1.0, -1.0,  1.0);
  Point p3(-1.0,  1.0,  1.0);
  Point p4(-1.0,  1.0, -1.0);
  Point p5( 1.0, -1.0, -1.0);
  Point p6( 1.0, -1.0,  1.0);
  Point p7( 1.0,  1.0,  1.0);
  Point p8( 1.0,  1.0, -1.0);

  add_crease(p1, p2, creases);
  add_crease(p2, p3, creases);
  add_crease(p3, p4, creases);
  add_crease(p4, p1, creases);

  add_crease(p5, p6, creases);
  add_crease(p6, p7, creases);
  add_crease(p7, p8, creases);
  add_crease(p8, p5, creases);

  add_crease(p5, p1, creases);
  add_crease(p6, p2, creases);
  add_crease(p7, p3, creases);
  add_crease(p8, p4, creases);

  domain.add_features(creases.begin(), creases.end());
#endif

  return make_mesh<C3t3>(domain,
                         facet_sizing, facet_approx, facet_angle,
                         cell_sizing, cell_shape,
                         function_name);
}

// logs the parameters, and dispatch to polyhedral, implicit or image domain
Exit_code run(const std::string& input,
              double facet_approx, double facet_sizing, double facet_angle,
              double cell_sizing, double cell_shape,
              int num_iteration)
{
  std::string domain = input;
  size_t slash_index = domain.find_last_of('/');
  if(slash_index == std::string::npos)
    slash_index = domain.find_last_of('\\');
  if(slash_index == std::string::npos)
    slash_index = 0;
  else
    ++slash_index;

  domain = domain.substr(slash_index, domain.find_last_of('.') - slash_index);

  CGAL_MESH_3_SET_PERFORMANCE_DATA("Domain", domain);
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Facet_size", facet_sizing);
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Facet_approx", facet_approx);
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Facet_angle", facet_angle);
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Cell_size", cell_sizing);
  CGAL_MESH_3_SET_PERFORMANCE_DATA("Cell_shape", cell_shape);
  xml_perf_set_technique();

  Exit_code res;
  for(int j=0; j<num_iteration; ++j)
  {
    std::cout << std::endl << "Refinement #" << j << "..." << std::endl;

    if(input == "Klein_function")
    {
      res = make_mesh_implicit(input, Klein_function(),
                               facet_sizing, facet_approx, facet_angle,
                               cell_sizing, cell_shape);
    }
    else if(input == "Tanglecube_function")
    {
      res = make_mesh_implicit(input, Tanglecube_function(),
                               facet_sizing, facet_approx, facet_angle,
                               cell_sizing, cell_shape);
    }
    else if(input == "Sphere_function")
    {
      Sphere_function f(1.);
      res = make_mesh_implicit(input, f,
                               facet_sizing, facet_approx, facet_angle,
                               cell_sizing, cell_shape);
    }
    else if(input == "Thin_cylinder_function")
    {
      Cylinder_function f(0.05, 3.);
      res = make_mesh_implicit(input, f,
                               facet_sizing, facet_approx, facet_angle,
                               cell_sizing, cell_shape);
    }
    else if(input == "Pancake_function")
    {
      Cylinder_function f(3., 0.1);
      res = make_mesh_implicit(input, f,
                               facet_sizing, facet_approx, facet_angle,
                               cell_sizing, cell_shape);
    }
    else // not a known implicit function
    {
      // Assume it's a CGAL data file...
      std::string full_path = CGAL::data_file_path(input);

      // ...and if it is not, then take it as it is
      if(!std::ifstream(full_path).good())
        full_path = input;

      const std::string extension = CGAL::IO::internal::get_file_extension(input);
      if(extension == "inr")
      {
        res = make_mesh_3D_images(full_path,
                                  facet_sizing, facet_approx, facet_angle,
                                  cell_sizing, cell_shape);
      }
      else // assume it's a polyhedron
      {
        res = make_mesh_polyhedron(full_path,
                                   facet_sizing, facet_approx, facet_angle,
                                   cell_sizing, cell_shape);
      }
    }

    std::cout << "Refinement #" << j << " done." << std::endl;
    std::cout << std::endl << "---------------------------------" << std::endl << std::endl;

    XML_perf_data::commit();
  }

  return res;
}

// reads filename & input parameters either from argv or from a script file
int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  // for the default xml filename, check XML_perf_data::build_filename()
  if(argc > 1)
    XML_perf_data::default_filename = argv[1];

#if defined(CHECK_MEMORY_LEAKS_ON_MSVC) && defined(_MSC_VER)
  _CrtSetDbgFlag ( _CRTDtbbBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF );
#endif

#ifdef CGAL_CONCURRENT_MESH_3
  Concurrent_mesher_config::load_config_file(CONCURRENT_MESHER_CONFIG_FILENAME, true);

  int max_nb_threads = Concurrent_mesher_config::get().num_threads;
  if(max_nb_threads == -1) // if not set in the config file, take the max available
    max_nb_threads = tbb::this_task_arena::max_concurrency();
#endif

  Exit_code res;

#ifdef CGAL_CONCURRENT_MESH_3
 #ifdef BENCHMARK_WITH_1_TO_MAX_THREADS
  for(int num_threads=1; num_threads<=max_nb_threads; ++num_threads)
 #else
  int num_threads = max_nb_threads;
 #endif // BENCHMARK_WITH_1_TO_MAX_THREADS
#endif // CGAL_CONCURRENT_MESH_3
  {
#ifdef CGAL_CONCURRENT_MESH_3
    std::cout << "-- Parallel Mesh_3 --" << std::endl;
    tbb::global_control control(tbb::global_control::max_allowed_parallelism, num_threads);
    display_info(num_threads);

    CGAL_MESH_3_SET_PERFORMANCE_DATA("Num_threads", num_threads);
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Lockgrid_size", Concurrent_mesher_config::get().locking_grid_num_cells_per_axis);
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Lock_radius", Concurrent_mesher_config::get().first_grid_lock_radius);
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Statgrid_size", Concurrent_mesher_config::get().work_stats_grid_num_cells_per_axis);
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Num_work_items_per_batch", Concurrent_mesher_config::get().num_work_items_per_batch);
#else
    std::cout << "-- Sequential Mesh_3 --" << std::endl;

    CGAL_MESH_3_SET_PERFORMANCE_DATA("Num_threads", "N/A");
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Lockgrid_size", "N/A");
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Lock_radius", "N/A");
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Statgrid_size", "N/A");
    CGAL_MESH_3_SET_PERFORMANCE_DATA("Num_work_items_per_batch", "N/A");
#endif

    // Script file format: each line gives
    //    - Filename (polyhedron and image) or "XXX_function" (implicit)
    //    - Facet and cell criteria
    //    - Number of iterations with these parameters
    std::ifstream script_file;
    script_file.open(BENCHMARK_INPUTS_FILENAME);
    if(script_file.is_open())
    {
      std::cout << "Found inputs file '" << BENCHMARK_INPUTS_FILENAME << "'" << std::endl;

      std::string line;
      while(std::getline(script_file, line))
      {
        if(line.empty() || line[0] == '#') // lines starting with '#' are ignored
          continue;

        std::cout << std::endl << std::endl;
        std::cout << "*****************************************" << std::endl;
        std::cout << "******* " << line << std::endl;
        std::cout << "*****************************************" << std::endl;

        std::stringstream sstr(line);

        std::string input;
        double facet_approx, facet_sizing, facet_angle;
        double cell_sizing, cell_shape;
        int num_iteration;
        if(!(sstr >> input
                  >> facet_approx
                  >> facet_sizing
                  >> facet_angle
                  >> cell_sizing
                  >> cell_shape
                  >> num_iteration))
        {
          std::cerr << "Error: failed to read input" << std::endl;
          return INPUT_IS_INVALID;
        }

        res = run(input,
                  facet_approx, facet_sizing, facet_angle,
                  cell_sizing, cell_shape,
                  num_iteration);
      }
    }
    else // no script
    {
      std::cout << "Inputs file '" << BENCHMARK_INPUTS_FILENAME << "' NOT found." << std::endl;

      // If the script is empty, use the command line arguments:
      // [this_program]
      // - filename
      // - facet_sizing
      // - facet_approx
      // - facet_angle
      // - cell_sizing
      // - cell_shape
      // - num_iteration

      std::string input = (argc > 2) ? argv[2] : "Sphere_function"; // @fixme klein assertion?
      double facet_sizing = (argc > 3) ? std::stod(argv[3]) : DEFAULT_FACE_SIZE;
      double facet_approx = (argc > 4) ? std::stod(argv[4]) : DEFAULT_FACE_APPROX;
      double facet_angle = (argc > 5) ? std::stod(argv[5]) : DEFAULT_FACE_ANGLE; // 25°
      double cell_sizing = (argc > 6) ? std::stod(argv[6]) : DEFAULT_CELL_SIZE;
      double cell_shape = (argc > 7) ? std::stod(argv[7]) : DEFAULT_CELL_RADIUS_RATIO; // 3

      if(facet_sizing <= 0)
        facet_sizing = DEFAULT_FACE_SIZE;
      if(facet_approx <= 0)
        facet_approx = DEFAULT_FACE_APPROX;
      if(facet_angle <= 0)
        facet_angle = DEFAULT_FACE_ANGLE;
      if(cell_sizing <= 0)
        cell_sizing = DEFAULT_CELL_SIZE;
      if(cell_shape <= 0)
        cell_shape = DEFAULT_CELL_RADIUS_RATIO;

      res = run(input,
                facet_approx, facet_sizing, facet_angle,
                cell_sizing, cell_shape,
                1 /*num_iteration*/);
    }
  }

  return res;
}
