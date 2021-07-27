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

#include <boost/histogram.hpp>
#include <boost/format.hpp>
#include <boost/filesystem.hpp>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/distance.h>


typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

typedef typename boost::graph_traits<Surface_mesh>::edge_descriptor edge_descriptor;
typedef typename boost::graph_traits<Surface_mesh>::face_descriptor face_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;
namespace fs = boost::filesystem;

typedef SMS::GarlandHeckbert_probabilistic_policies<Surface_mesh, Kernel> Prob_plane;
typedef SMS::GarlandHeckbert_policies<Surface_mesh, Kernel> Classic_plane;
typedef SMS::GarlandHeckbert_probabilistic_tri_policies<Surface_mesh, Kernel> Prob_tri;
typedef SMS::GarlandHeckbert_triangle_policies<Surface_mesh, Kernel> Classic_tri;

using hist = boost::histogram::histogram<std::tuple<boost::histogram::axis::regular<double, 
      boost::use_default, boost::use_default, boost::use_default>>>; 

// settings for benchmarking - throw away the first n_burns results and keep the n_samples
// samples
constexpr int n_burns = 1;
constexpr int n_samples = 10;

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
  
  for (int i = 0; i < n_burns; ++i)
  {
    for (Surface_mesh& mesh : temp_meshes)
    {
      edge_collapse<Policy>(mesh);
    }
  }
  temp_meshes = meshes;

  // measure each run 
  for (int i = 0; i < n_samples; ++i)
  {
    // measure time taken by the edge_collapse function over all meshes
    time::time_point start_time = time::steady_clock::now();
  
    for (Surface_mesh& mesh : temp_meshes)
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
    for (fs::directory_iterator dir_iter(dir); dir_iter != end_iter; ++dir_iter)
    {
      // look for .off files
      if (fs::is_regular_file(dir_iter->status())
          && dir_iter->path().extension() == ".off")
      {
        const fs::path& path = dir_iter->path();

        // input stream from the file
        fs::ifstream is(path);

        // try to read the file into the surface mesh, upon failure, we continue looping 
        // iterating through the directory
        if (!is || !(is >> mesh))
        {
          std::cerr << "Failed to read input mesh " << path 
            << " in directory " << dir << std::endl;
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
  std::for_each(std::begin(vec), std::end(vec), [] (const std::pair<Surface_mesh, std::string>&
        p) { return p.first; });
  return time_all_vector<Policy>(vec, output_file);
}

// taken from a boost example
template<typename histo>
void print_hist(histo h)
{
  using namespace boost::histogram;
  
  std::ostringstream os;
  for (auto&& x : indexed(h, coverage::all)) 
  {
    // alternative output formatting
    //os << boost::format("bin %2i [%4.1f, %4.1f): %i\n")
    //      % x.index() % x.bin().lower() % x.bin().upper() % *x;
    os << boost::format("%i, ")
          % *x;
  }

  std::cout << os.str() << std::endl;
}

hist generate_edge_statistics(Surface_mesh mesh, hist histo, std::function<FT(edge_descriptor, 
      const Surface_mesh&)> f)
{
  for (auto e : mesh.edges())
  {
    FT value = f(e, mesh); 
    histo(value);
  }

  return histo;
}

hist generate_face_statistics(Surface_mesh mesh, hist histo, std::function<FT(face_descriptor, 
      const Surface_mesh&)> f)
{
  for (auto face : mesh.faces())
  {
    FT value = f(face, mesh); 
    histo(value);
  }

  return histo;
}

template<typename Policy>
double hausdorff_error(const Surface_mesh& mesh, double ratio = 0.2)
{
  // an arbitrarily chosen small value 
  constexpr double error_bound = 0.00001;
  
  // make a copy of the mesh so we can compare later
  Surface_mesh temp = mesh;
 
  edge_collapse<Policy>(temp);
  return PMP::bounded_error_symmetric_Hausdorff_distance<CGAL::Sequential_tag>
    (mesh, temp, error_bound);
}
// calculate approximate Hausdorff errors for all different policies at the same
// decimation ratio
std::array<double, 4> hausdorff_errors(const Surface_mesh& mesh, double ratio = 0.2)
{
  std::array<double, 4> ret { {0, 0, 0, 0} }; 

  ret[0] = hausdorff_error<Classic_plane>(mesh, ratio);
  ret[1] = hausdorff_error<Prob_plane>(mesh, ratio);
  ret[2] = hausdorff_error<Classic_tri>(mesh, ratio);
  ret[3] = hausdorff_error<Prob_tri>(mesh, ratio);
    
  return ret;
}

void hausdorff_errors_dir(const fs::path& dir, const fs::path& out, double ratio = 0.2)
{
  // get a vector of pairs of a mesh and its name by loading all .off files in a directory
  const auto meshes = get_all_meshes(dir); 
  
  fs::ofstream out_stream {out};
  
  for (auto& p : meshes)
  {
    // calculate the errors and output to file
    std::array<double, 4> errors = hausdorff_errors(p.first, ratio);

    out_stream << p.second << '\n';
    for (double err : errors)
    {
      out_stream << std::to_string(err) << '\n';
    }
  }
}

int main(int argc, char** argv)
{
  using namespace boost::histogram;
  hist edge_histo = make_histogram(axis::regular<>(14, 0.1, 4.0, "length"));
  
  auto f = [] (edge_descriptor e, const Surface_mesh& mesh) { 
    return PMP::edge_length(mesh.halfedge(e), mesh);
  };
  
  auto g = [] (face_descriptor face, const Surface_mesh& mesh) { 
    return PMP::face_aspect_ratio(face, mesh);
  };
  
    
  const fs::path output {"hausdorff_errors"};

  if (fs::exists(output))
  {
    // clear file
    fs::resize_file(output, 0); 
  }

  hausdorff_errors_dir("data/", output);
  
  //time_policy<Classic_tri>("../data/", output);
  
  /*
  Surface_mesh probabilistic = edge_collapse
    <SMS::GarlandHeckbert_probabilistic_policies<Surface_mesh, Kernel>>(surface_mesh);
  
  Surface_mesh classic = edge_collapse
    <SMS::GarlandHeckbert_policies<Surface_mesh, Kernel>>(surface_mesh);
  
  Surface_mesh classic_tri = edge_collapse
    <SMS::GarlandHeckbert_triangle_policies<Surface_mesh, Kernel>>(surface_mesh);
 
  Surface_mesh probabilistic_tri = edge_collapse
    <SMS::GarlandHeckbert_probabilistic_tri_policies<Surface_mesh, Kernel>>(surface_mesh);
  
  std::cout << "Original mesh histogram:\n";
  print_hist(generate_edge_statistics(surface_mesh, edge_histo, f));

  std::cout << "Decimated classic mesh histogram:\n";
  print_hist(generate_edge_statistics(classic, edge_histo, f));
 
  std::cout << "Decimated probabilistic mesh histogram:\n";
  print_hist(generate_edge_statistics(probabilistic, edge_histo, f));

  std::cout << "Decimated classic tri mesh histogram:\n";
  print_hist(generate_edge_statistics(classic_tri, edge_histo, f));

  std::cout << "Decimated probabilistic tri mesh histogram:\n";
  print_hist(generate_edge_statistics(probabilistic_tri, edge_histo, f));

  hist face_histo = make_histogram(axis::regular<>(10, 1.0, 5.5));
  
  std::cout << "Original aspect ratio:\n";
  print_hist(generate_face_statistics(surface_mesh, face_histo, g));

  std::cout << "Classic aspect ratio:\n";
  print_hist(generate_face_statistics(classic, face_histo, g));
  
  std::cout << "Probabilistic aspect ratio:\n";
  print_hist(generate_face_statistics(probabilistic, face_histo, g));
  
  std::cout << "Classic tri aspect ratio:\n";
  print_hist(generate_face_statistics(classic_tri, face_histo, g));
  
  std::cout << "Probabilistic tri aspect ratio:\n";
  print_hist(generate_face_statistics(probabilistic_tri, face_histo, g));
  */

  
  return EXIT_SUCCESS;
}
