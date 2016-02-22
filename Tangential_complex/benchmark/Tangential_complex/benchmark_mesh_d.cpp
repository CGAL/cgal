//#undef CGAL_LINKED_WITH_TBB // CJTODO TEMP

// Without TBB_USE_THREADING_TOOL Intel Inspector XE will report false positives in Intel TBB
// (http://software.intel.com/en-us/articles/compiler-settings-for-threading-error-analysis-in-intel-inspector-xe/)
#ifdef _DEBUG
# define TBB_USE_THREADING_TOOL
#endif

#include <cstddef>

const std::size_t ONLY_LOAD_THE_FIRST_N_POINTS = 100000000;

#include <CGAL/assertions_behaviour.h>
#include <CGAL/Mesh_d.h>
#include <CGAL/Epick_d.h>
#include <CGAL/Random.h>
#include <CGAL/Mesh_3/Profiling_tools.h>

#include "../../test/Tangential_complex/testing_utilities.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/range/adaptor/strided.hpp>

#include <cstdlib>
#include <ctime>
#include <fstream>
#include <math.h>

#ifdef CGAL_LINKED_WITH_TBB
# include <tbb/task_scheduler_init.h>
#endif

using namespace CGAL::Tangential_complex_;

const char * const BENCHMARK_SCRIPT_FILENAME = "benchmark_script.txt";

typedef CGAL::Epick_d<CGAL::Dynamic_dimension_tag>              Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_d                                         Point;
typedef CGAL::Mesh_d<Kernel, Kernel, CGAL::Parallel_tag>        Mesh;

//#define CHECK_IF_ALL_SIMPLICES_ARE_IN_THE_AMBIENT_DELAUNAY
//#define MESH_D_INPUT_STRIDES 10 // only take one point every MESH_D_INPUT_STRIDES points
//#define MESH_D_NO_EXPORT
//#define CGAL_MESH_D_USE_LINEAR_PROG_TO_COMPUTE_INTERSECTION
//#define CGAL_MESH_D_FILTER_BY_TESTING_ALL_VERTICES_TANGENT_PLANES


template <typename Mesh, typename Indexed_simplex_range = void>
bool export_to_off(
  Mesh const& mesh,
  std::string const& input_name_stripped,
  std::string const& suffix,
  Indexed_simplex_range const *p_simpl_to_color_in_red = NULL,
  Indexed_simplex_range const *p_simpl_to_color_in_green = NULL,
  Indexed_simplex_range const *p_simpl_to_color_in_blue = NULL)
{
#ifdef MESH_D_NO_EXPORT
  return true;
#endif
  if (mesh.intrinsic_dimension() <= 3)
  {
    std::stringstream output_filename;
    output_filename << "output/" << input_name_stripped << "_"
      << mesh.intrinsic_dimension() << "_in_R"
      << mesh.ambient_dimension() << "_"
      << mesh.number_of_vertices() << "v"
      << suffix << ".off";

    std::ofstream off_stream(output_filename.str().c_str());

      mesh.export_to_off(
        off_stream,
        p_simpl_to_color_in_red,
        p_simpl_to_color_in_green,
        p_simpl_to_color_in_blue);
    return true;
  }

  return false;
}


void make_mesh(
  std::vector<Point> &points, 
  int intrinsic_dim,
  bool sparsify = true,
  double sparsity = 0.01,
  const char *input_name = "mesh")
{
  Kernel k;
  Kernel lk; // local kernel (intrinsic dim)

  //===========================================================================
  // Init
  //===========================================================================
  Wall_clock_timer t;

  // Get input_name_stripped
  std::string input_name_stripped(input_name);
  size_t slash_index = input_name_stripped.find_last_of('/');
  if (slash_index == std::string::npos)
    slash_index = input_name_stripped.find_last_of('\\');
  if (slash_index == std::string::npos)
    slash_index = 0;
  else
    ++slash_index;
  input_name_stripped = input_name_stripped.substr(
    slash_index, input_name_stripped.find_last_of('.') - slash_index);

  int ambient_dim = k.point_dimension_d_object()(*points.begin());

#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  std::vector<Point> points_not_sparse = points;
#endif

  //===========================================================================
  // Sparsify point set if requested
  //===========================================================================
  if (sparsify)
  {
    std::size_t num_points_before = points.size();
    points = sparsify_point_set(k, points, sparsity*sparsity);
    std::cerr << "Number of points before/after sparsification: "
      << num_points_before << " / " << points.size() << "\n";
  }

  //===========================================================================
  // Compute Tangential Complex
  //===========================================================================

#ifdef CGAL_TC_USE_ANOTHER_POINT_SET_FOR_TANGENT_SPACE_ESTIM
  Mesh mesh(points.begin(), points.end(), sparsity, intrinsic_dim,
    points_not_sparse.begin(), points_not_sparse.end(), k, lk);
#else
  Mesh mesh(points.begin(), points.end(), sparsity, intrinsic_dim, k, lk);
#endif

  double init_time = t.elapsed(); t.reset();

  std::vector<std::set<std::size_t> > uncertain_simplices;
  mesh.compute_mesh(&uncertain_simplices);
  double computation_time = t.elapsed(); t.reset();

  mesh.display_stats();

  //=========================================================================
  // Export to OFF
  //=========================================================================
  t.reset();
  bool exported = export_to_off(
    mesh, input_name_stripped, "", &uncertain_simplices);
  double export_time = (exported ? t.elapsed() : -1);
  t.reset();

  //===========================================================================
  // Is the result a pure pseudomanifold?
  //===========================================================================
  std::size_t num_wrong_dim_simplices, 
              num_wrong_number_of_cofaces, 
              num_unconnected_stars;
  std::set<std::set<std::size_t> > wrong_dim_simplices;
  std::set<std::set<std::size_t> > wrong_number_of_cofaces_simplices;
  std::set<std::set<std::size_t> > unconnected_stars_simplices;
  bool is_pure_pseudomanifold = mesh.complex().is_pure_pseudomanifold(
    intrinsic_dim, mesh.number_of_vertices(), false, false, 1,
    &num_wrong_dim_simplices, &num_wrong_number_of_cofaces, 
    &num_unconnected_stars,
    &wrong_dim_simplices, &wrong_number_of_cofaces_simplices, 
    &unconnected_stars_simplices);

  // Stats about the simplices
  mesh.complex().display_stats();
  std::size_t num_edges = mesh.complex().num_K_simplices<1>();
  std::size_t num_triangles = mesh.complex().num_K_simplices<2>();
  std::cerr << "Euler caract.: V - E + F = "
    << mesh.number_of_vertices()
    << " - " << (std::ptrdiff_t) num_edges
    << " + " << (std::ptrdiff_t) num_triangles
    << " = "
    << yellow
    << (std::ptrdiff_t) mesh.number_of_vertices()
    - (std::ptrdiff_t) num_edges
    + (std::ptrdiff_t) num_triangles
    << white << "\n";

  //===========================================================================
  // Display info
  //===========================================================================

  std::cerr
    << "\n================================================\n"
    << "Number of vertices: " << mesh.number_of_vertices() << "\n"
    << "Pure pseudomanifold: " << yellow 
    << (is_pure_pseudomanifold ? "YES" : "NO") << white << "\n"
    << "Computation times (seconds): \n"
    << "  * Mesh: " << init_time + computation_time << "\n"
    << "    - Init + kd-tree = " << init_time << "\n"
    << "    - Mesh computation = " << computation_time << "\n"
    //<< "  * Export to OFF : " << export_before_time << "\n"
    << "================================================\n\n";
}

int main()
{
  CGAL::set_error_behaviour(CGAL::ABORT);

#ifdef CGAL_LINKED_WITH_TBB
# ifdef _DEBUG
  int num_threads = 1;
# else
  int num_threads = 8;
# endif
#endif

  unsigned int seed = static_cast<unsigned int>(time(NULL));
  CGAL::default_random = CGAL::Random(seed);
  std::cerr << "Random seed = " << seed << "\n";

  std::ifstream script_file;
  script_file.open(BENCHMARK_SCRIPT_FILENAME);
  // Script?
  // Script file format: each line gives
  //    - Filename (point set) or "generate_XXX" (point set generation)
  //    - Ambient dim
  //    - Intrinsic dim
  //    - Number of iterations with these parameters
  if (script_file.is_open())
  {
    int i = 1;
#ifdef CGAL_LINKED_WITH_TBB
# ifdef BENCHMARK_WITH_1_TO_MAX_THREADS
    for(num_threads = 1 ;
          num_threads <= tbb::task_scheduler_init::default_num_threads() ;
          ++num_threads)
# endif
#endif
    {
#ifdef CGAL_LINKED_WITH_TBB
      tbb::task_scheduler_init init(
        num_threads > 0 ? num_threads : tbb::task_scheduler_init::automatic);
#endif

      std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME 
        << "' found.\n";
      script_file.seekg(0);
      while (script_file.good())
      {
        std::string line;
        std::getline(script_file, line);
        if (line.size() > 1 && line[0] != '#')
        {
          boost::replace_all(line, "\t", " ");
          boost::trim_all(line);
          std::cerr << "\n\n";
          std::cerr << "*****************************************\n";
          std::cerr << "******* " << line << "\n";
          std::cerr << "*****************************************\n";
          std::stringstream sstr(line);

          std::string input;
          std::string param1;
          std::string param2;
          std::string param3;
          std::size_t num_points;
          int ambient_dim;
          int intrinsic_dim;
          char sparsify;
          double sparsity;
          char perturb, add_high_dim_simpl, collapse;
          double time_limit_for_perturb;
          int num_iteration;
          sstr >> input;
          sstr >> param1;
          sstr >> param2;
          sstr >> param3;
          sstr >> num_points;
          sstr >> ambient_dim;
          sstr >> intrinsic_dim;
          sstr >> sparsify;
          sstr >> sparsity;
          sstr >> perturb;
          sstr >> add_high_dim_simpl;
          sstr >> collapse;
          sstr >> time_limit_for_perturb;
          sstr >> num_iteration;

          for (int j = 0 ; j < num_iteration ; ++j)
          {
            std::string input_stripped = input;
            size_t slash_index = input_stripped.find_last_of('/');
            if (slash_index == std::string::npos)
              slash_index = input_stripped.find_last_of('\\');
            if (slash_index == std::string::npos)
              slash_index = 0;
            else
              ++slash_index;
            input_stripped = input_stripped.substr(
              slash_index, input_stripped.find_last_of('.') - slash_index);

            std::cerr << "\nMesh #" << i << "...\n";
          
#ifdef CGAL_MESH_D_PROFILING
            Wall_clock_timer t_gen;
#endif

            std::vector<Point> points;
            Mesh::TS_container tangent_spaces;

            if (input == "generate_moment_curve")
            {
              points = generate_points_on_moment_curve<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_plane")
            {
              points = generate_points_on_plane<Kernel>(
                num_points, intrinsic_dim, ambient_dim);
            }
            else if (input == "generate_sphere_d")
            {
              points = generate_points_on_sphere_d<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()),  // radius
                std::atof(param2.c_str())); // radius_noise_percentage
            }
            else if (input == "generate_two_spheres_d")
            {
              points = generate_points_on_two_spheres_d<Kernel>(
                num_points, ambient_dim,
                std::atof(param1.c_str()),
                std::atof(param2.c_str()),
                std::atof(param3.c_str()));
            }
            else if (input == "generate_3sphere_and_circle_d")
            {
              CGAL_assertion(intrinsic_dim == 3);
              CGAL_assertion(ambient_dim == 5);
              points = generate_points_on_3sphere_and_circle<Kernel>(
                num_points,
                std::atof(param1.c_str()));
            }
            else if (input == "generate_torus_3D")
            {
              points = generate_points_on_torus_3D<Kernel>(
                num_points,
                std::atof(param1.c_str()),
                std::atof(param2.c_str()),
                param3 == "Y");
            }
            else if (input == "generate_torus_d")
            {
              points = generate_points_on_torus_d<Kernel>(
                num_points, 
                intrinsic_dim,
                param1 == "Y", // uniform
                std::atof(param2.c_str())); // radius_noise_percentage
            }
            else if (input == "generate_klein_bottle_3D")
            {
              points = generate_points_on_klein_bottle_3D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_klein_bottle_4D")
            {
              points = generate_points_on_klein_bottle_4D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else if (input == "generate_klein_bottle_variant_5D")
            {
              points = generate_points_on_klein_bottle_variant_5D<Kernel>(
                num_points,
                std::atof(param1.c_str()), std::atof(param2.c_str()));
            }
            else
            {
              load_points_from_file<Kernel, typename Mesh::Tangent_space_basis>(
                input, std::back_inserter(points),
                std::back_inserter(tangent_spaces),
                ONLY_LOAD_THE_FIRST_N_POINTS);
            }

#ifdef CGAL_MESH_D_PROFILING
            std::cerr << "Point set generated/loaded in " << t_gen.elapsed()
                      << " seconds.\n";
#endif

            if (!points.empty())
            {
#if defined(MESH_D_INPUT_STRIDES) && MESH_D_INPUT_STRIDES > 1
              auto p = points | boost::adaptors::strided(MESH_D_INPUT_STRIDES); // CJTODO C++11 (auto)
              std::vector<Point> points(p.begin(), p.end());
              std::cerr << "****************************************\n"
                << "WARNING: taking 1 point every " << MESH_D_INPUT_STRIDES
                << " points.\n"
                << "****************************************\n";
#endif
              make_mesh(
                points, intrinsic_dim, sparsify == 'Y', sparsity, 
                input.c_str());

              std::cerr << "Mesh #" << i++ << " done.\n";
              std::cerr << "\n---------------------------------\n\n";
            }
            else
            {
              std::cerr << "Mesh #" << i++ << ": no points loaded.\n";
            }
          }
        }
      }
      script_file.seekg(0);
      script_file.clear();
    }

    script_file.close();
  }
  // Or not script?
  else
  {
    std::cerr << "Script file '" << BENCHMARK_SCRIPT_FILENAME << "' NOT found.\n";
  }

  system("pause");
  return 0;
}
