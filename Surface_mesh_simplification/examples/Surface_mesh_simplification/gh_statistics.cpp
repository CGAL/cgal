#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_probabilistic_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_probabilistic_tri_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_triangle_policies.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>

#include <iostream>
#include <fstream>

#include <chrono>
#include <vector>
#include <array>
#include <functional>

#include <boost/format.hpp> 
#include <boost/filesystem.hpp>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

/**
 * this is a largely undocumented file that was used to gather most statistics
 * for my GSoC 2021 project on mesh decimation, maybe it is useful in the future  
 */

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

typedef typename boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
typedef typename boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = boost::filesystem;

typedef SMS::GarlandHeckbert_policies<Surface_mesh, Kernel> Classic_plane;
typedef SMS::GarlandHeckbert_probabilistic_policies<Surface_mesh, Kernel> Prob_plane;
typedef SMS::GarlandHeckbert_triangle_policies<Surface_mesh, Kernel> Classic_tri;
typedef SMS::GarlandHeckbert_probabilistic_tri_policies<Surface_mesh, Kernel> Prob_tri;

// settings for benchmarking - throw away the first n_burns results and keep the n_samples
// samples
constexpr int n_burns = 1;
constexpr int n_samples = 10;

constexpr std::size_t classic_plane_index = 0;
constexpr std::size_t prob_plane_index = 1;
constexpr std::size_t classic_tri_index = 2;
constexpr std::size_t prob_tri_index = 3;

template<typename Policy>
Surface_mesh edge_collapse(Surface_mesh& mesh, double ratio = 0.2)
{
  typedef typename Policy::Get_cost Cost; 
  typedef typename Policy::Get_placement Placement;
  
  typedef SMS::Bounded_normal_change_placement<Placement> Bounded_placement;
  
  const Policy p {mesh, 100};
  
  const Cost& cost = p.get_cost();
  const Placement& unbounded_placement = p.get_placement();

  Bounded_placement placement(unbounded_placement);

  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);
  
  // collapse edges ignoring result code
  SMS::edge_collapse(mesh, stop,
      CGAL::parameters::get_cost(cost).get_placement(placement));
  
  return mesh;
}

template<typename Policy>
void time_all_vector(const std::vector<Surface_mesh>& meshes, const fs::path& output_file)
{
  namespace time = std::chrono;

  fs::ofstream out_stream {output_file};
  
  // always decimate meshes in this vector to avoid timing the
  // copying of the meshes
  std::vector<Surface_mesh> temp_meshes = meshes;
  
  for(int i = 0; i < n_burns; ++i)
  {
    for(Surface_mesh& mesh : temp_meshes)
    {
      edge_collapse<Policy>(mesh);
    }
  }
  temp_meshes = meshes;

  // measure each run 
  for(int i = 0; i < n_samples; ++i)
  {
    // measure time taken by the edge_collapse function over each mesh individually
    time::time_point start_time = time::steady_clock::now();

    for(Surface_mesh& mesh : temp_meshes)
    {
      edge_collapse<Policy>(mesh);
    }

    time::time_point end_time = time::steady_clock::now();

    // get elapsed time in nanoseconds
    unsigned long elapsed_ns = time::duration_cast<time::nanoseconds>
      (end_time - start_time).count();

    // add the elapsed time as data point
    out_stream << std::to_string(elapsed_ns) << '\n';   

    temp_meshes = meshes;
  }
}

bool read_from_file(Surface_mesh& mesh, const fs::path& p)
{
  fs::ifstream in {p};

  // make sure the mesh is empty and ready for the next set of data
  mesh.clear();

  if (!in || !(in >> mesh))
  {
    return false;
  }

  return true;
}

std::vector<std::pair<Surface_mesh, std::string>> get_all_meshes(const fs::path& dir)
{
  // vector for storing results 
  std::vector<std::pair<Surface_mesh, std::string>> meshes { };

  // read all meshes in the given directory into the vector
  if (fs::is_directory(dir))
  {
    fs::directory_iterator end_iter;
    Surface_mesh mesh;

    // iterate through all files in the given directory
    for(fs::directory_iterator dir_iter(dir); dir_iter != end_iter; ++dir_iter)
    {
      // look for .off files
      if (fs::is_regular_file(dir_iter->status())
          && dir_iter->path().extension() == ".off")
      {
        const fs::path& path = dir_iter->path();

        // try to read the file into the surface mesh, upon failure, we continue looping 
        // iterating through the directory
        if (!read_from_file(mesh, path))
        {
          std::cerr << "Failed to read input mesh " << path 
            << " in directory " << dir << "." << std::endl;
        }
        else if (!CGAL::is_triangle_mesh(mesh))
        {
          std::cerr << "Input geometry is not triangulated." << std::endl;
        }
        else
        {
          // we can move this mesh since we don't need it anymore - the next one will be read into 
          // the same variable
          meshes.emplace_back(std::move(mesh), path.filename().string());
        }
      }
    }
  }

  return meshes;
}

  template<typename Policy>
void time_policy(const fs::path& dir, const fs::path& output_file)
{
  auto vec = get_all_meshes(dir);
  std::vector<Surface_mesh> meshes { };

  for(const auto& p : vec) {
    std::cout << p.second << '\n';
  }

  // keep the first values from the pair (we don't need the names)(we don't need the names)
  std::transform(vec.begin(), vec.end(), std::back_inserter(meshes),
      [] (const std::pair<Surface_mesh, std::string>& p) { return p.first; }
      );

  return time_all_vector<Policy>(meshes, output_file);
}

  template<typename Policy>
double hausdorff_error(const Surface_mesh& mesh, double ratio = 0.2)
{
  // an arbitrarily chosen small value 
  constexpr double error_bound = 0.00001;

  // make a copy of the mesh so we can compare later
  Surface_mesh temp = mesh;

  edge_collapse<Policy>(temp, ratio);
  return PMP::bounded_error_symmetric_Hausdorff_distance<CGAL::Sequential_tag>
    (mesh, temp, error_bound);
}

// calculate approximate Hausdorff errors for all different policies at 
// the same decimation ratio
std::array<FT, 4> hausdorff_errors(const Surface_mesh& mesh, double ratio)
{
  std::array<FT, 4> ret { {0, 0, 0, 0} }; 

  ret[classic_plane_index] = hausdorff_error<Classic_plane>(mesh, ratio);
  ret[prob_plane_index] = hausdorff_error<Prob_plane>(mesh, ratio);
  ret[classic_tri_index] = hausdorff_error<Classic_tri>(mesh, ratio);
  ret[prob_tri_index] = hausdorff_error<Prob_tri>(mesh, ratio);

  return ret;
}

  template<typename InputIt>
void hausdorff_errors_mesh(const Surface_mesh& mesh, const fs::path& out, InputIt begin, InputIt end)
{
  fs::ofstream out_stream {out};

  for(InputIt it = begin; it != end; ++it)
  {
    std::array<double, 4> errs = hausdorff_errors(mesh, *it);

    out_stream << "ratio: " << *it << '\n';

    for(double e : errs)
    {
      out_stream << std::to_string(e) << '\n';
    }
  }
}

  template<typename InputIt>
void hausdorff_errors_mesh(const fs::path& in, const fs::path& out, InputIt begin, InputIt end)
{
  fs::ifstream is {in};
  Surface_mesh mesh;

  if (!read_from_file(mesh, in))
  {
    std::cerr << "Failed to read input mesh " << in << '\n'; 
  }
  else
  {
    hausdorff_errors_mesh(mesh, out, begin, end);
  }
} 

void hausdorff_errors_dir(const fs::path& dir, const fs::path& out, double ratio = 0.2)
{
  // get a vector of pairs of a mesh and its name by loading all .off files in a directory
  const auto meshes = get_all_meshes(dir); 

  fs::ofstream out_stream {out};

  for(auto& p : meshes)
  {
    // calculate the errors and output to file
    std::array<double, 4> errors = hausdorff_errors(p.first, ratio);

    out_stream << p.second << '\n';
    for(double err : errors)
    {
      out_stream << std::to_string(err) << '\n';
    }
  }
}

void write_aspect_ratios(fs::ofstream& out, const Surface_mesh& mesh) 
{
  for(auto face : faces(mesh))
  {
    out << std::to_string(PMP::face_aspect_ratio(face, mesh)) << '\n'; 
  }
}

void gather_face_aspect_ratio(const fs::path& file, const fs::path& out)
{
  Surface_mesh cp { }; 
  Surface_mesh ct { }; 
  Surface_mesh pp { }; 
  Surface_mesh pt { }; 

  read_from_file(cp, file);
  read_from_file(ct, file);
  read_from_file(pp, file);
  read_from_file(pt, file);

  edge_collapse<Classic_plane>(cp);
  edge_collapse<Prob_plane>(pp);
  edge_collapse<Classic_tri>(ct);
  edge_collapse<Prob_tri>(pt);

  fs::ofstream out_stream { out };

  out_stream << "classic_plane\n";
  write_aspect_ratios(out_stream, cp);
  out_stream << "prob_plane\n";
  write_aspect_ratios(out_stream, pp);
  out_stream << "classic_tri\n";
  write_aspect_ratios(out_stream, ct);
  out_stream << "prob_tri\n";
  write_aspect_ratios(out_stream, pt);
}

enum class Mode { time, hausdorff_ratio, hausdorff_mesh, run, apect_ratio };

void time_policy_string(const fs::path& in, const fs::path& out, const std::string& policy) {
  if (policy == "classic_plane")
  {
    time_policy<Classic_plane>(in, out);
  }
  else if (policy == "prob_plane")
  {
    time_policy<Prob_plane>(in, out);
  }
  else if (policy == "classic_tri")
  {
    time_policy<Classic_tri>(in, out);
  }
  else if (policy == "prob_tri")
  {
    time_policy<Prob_tri>(in, out);
  }
}

int main(int argc, char* argv[])
{

  // default is time
  Mode m = Mode::time;

  if (argc > 1)
  {
    std::string input { argv[1] };

    if (input == "time")
    {
      m = Mode::time;

      if (argc != 5)
      {
        std::cerr << "Expected 3 arugments to time: [input_dir] [output_file] [policy].\n";
        return EXIT_FAILURE;
      }

      time_policy_string(argv[2], argv[3], argv[4]);
    }
    else if (input == "ratio")
    {
      m = Mode::hausdorff_ratio;
    }
    else if (input == "mesh")
    {
      m = Mode::hausdorff_mesh;
    }
    else if (input == "run")
    {
      m = Mode::run;
      if (argc != 3 && argc != 4)
      {
        std::cerr << "Expected 1 (2) argument(s) to run: [input_file] [ratio].\n";
        return EXIT_FAILURE;
      }
      // default init empty mesh
      Surface_mesh m { };

      FT ratio = 0.2;
      if (argc == 4) {
        ratio = std::stod(argv[3]);
      }

      read_from_file(m, argv[2]);

      Surface_mesh cp = m;
      Surface_mesh pp = m;
      Surface_mesh ct = m;
      Surface_mesh pt = m;

      edge_collapse<Classic_plane>(cp, ratio);
      edge_collapse<Prob_plane>(pp, ratio);
      edge_collapse<Classic_tri>(ct, ratio);
      edge_collapse<Prob_tri>(pt, ratio);

      fs::ofstream classic_plane_stream {"run_classic_plane.off"};
      fs::ofstream prob_plane_stream {"run_prob_plane.off"};
      fs::ofstream classic_tri_stream {"run_classic_tri.off"};
      fs::ofstream prob_tri_stream {"run_prob_tri.off"};

      classic_plane_stream << cp;
      prob_plane_stream << pp;
      classic_tri_stream << ct;
      prob_tri_stream << pt;
    }
    else if (input == "aspect")
    {
      gather_face_aspect_ratio(argv[2], argv[3]);
    } 
    else 
    {
      std::cerr << "Didn't recognise command line argument " << input << 
        "(options are time, ratio or mesh)\n"; 

      return EXIT_FAILURE;
    }
  }

  if (m == Mode::hausdorff_mesh)
  {
    if (argc != 3)
    {
      std::cerr << "Expected 1 argument to mesh: [input_file].\n";
      return EXIT_FAILURE;
    }
    else 
    {
      std::array<double, 10> range {{ 0.7, 0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15 }};
      hausdorff_errors_mesh(argv[2], "mesh_hausdorff_err", range.begin(), range.end());
    }
  }
}
