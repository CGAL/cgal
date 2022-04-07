#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>

#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <array>
#include <chrono>
#include <iostream>
#include <fstream>
#include <functional>
#include <vector>

typedef CGAL::Simple_cartesian<double>                                             Kernel;
typedef Kernel::FT                                                                 FT;
typedef Kernel::Point_3                                                            Point_3;
typedef CGAL::Surface_mesh<Point_3>                                                Surface_mesh;

typedef typename boost::graph_traits<Surface_mesh>::edge_descriptor                edge_descriptor;
typedef typename boost::graph_traits<Surface_mesh>::face_descriptor                face_descriptor;

namespace SMS = CGAL::Surface_mesh_simplification;
namespace PMP = CGAL::Polygon_mesh_processing;

typedef SMS::GarlandHeckbert_plane_policies<Surface_mesh, Kernel>                  Classic_plane;
typedef SMS::GarlandHeckbert_probabilistic_plane_policies<Surface_mesh, Kernel>    Prob_plane;
typedef SMS::GarlandHeckbert_triangle_policies<Surface_mesh, Kernel>               Classic_tri;
typedef SMS::GarlandHeckbert_probabilistic_triangle_policies<Surface_mesh, Kernel> Prob_tri;

// settings for benchmarking - throw away the first n_burns results and keep the n_samples samples
constexpr int n_burns = 1;
constexpr int n_samples = 3;

constexpr std::size_t classic_plane_index = 0;
constexpr std::size_t prob_plane_index = 1;
constexpr std::size_t classic_tri_index = 2;
constexpr std::size_t prob_tri_index = 3;

// =================================================================================================
// =================================================================================================
// =================================================================================================

bool read_from_file(const std::string& filename,
                    Surface_mesh& mesh)
{
  // make sure the mesh is empty and ready for the next set of data
  clear(mesh);

  if(!PMP::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Error: failed to read input: " << filename << std::endl;
    return false;
  }

  return true;
}

std::vector<std::pair<Surface_mesh, std::string> >
get_all_meshes(const std::vector<std::string>& filenames)
{
  // vector for storing results
  std::vector<std::pair<Surface_mesh, std::string> > meshes;

  for(const std::string& filename : filenames)
  {
    Surface_mesh mesh;
    if(!read_from_file(filename, mesh))
      std::cerr << "Error: Failed to read input mesh " << filename << std::endl;
    else if(!CGAL::is_triangle_mesh(mesh))
      std::cerr << "Error: Input geometry is not triangulated." << std::endl;
    else
      // we can move this mesh since we don't need it anymore
      // the next one will be read into the same variable
      meshes.emplace_back(std::move(mesh), filename);
  }

  return meshes;
}

// =================================================================================================
// =================================================================================================
// =================================================================================================

template <typename Policy>
Surface_mesh edge_collapse(Surface_mesh& mesh,
                           const double ratio = 0.2)
{
  typedef typename Policy::Get_cost Cost;
  typedef typename Policy::Get_placement Placement;
  typedef SMS::Bounded_normal_change_placement<Placement> Bounded_placement;

  std::cout << "Edge collapse mesh of " << num_edges(mesh) << " edges. Policy: " << typeid(Policy).name() << std::endl;

  const Policy p { mesh, 100 };

  const Cost& cost = p.get_cost();
  const Placement& unbounded_placement = p.get_placement();
  Bounded_placement bounded_placement(unbounded_placement);
  SMS::Count_ratio_stop_predicate<Surface_mesh> stop(ratio);

  std::chrono::time_point<std::chrono::steady_clock> start_time = std::chrono::steady_clock::now();

  SMS::edge_collapse(mesh, stop, CGAL::parameters::get_cost(cost)
                                                  .get_placement(unbounded_placement));

  std::chrono::time_point<std::chrono::steady_clock> end_time = std::chrono::steady_clock::now();

  auto elapsed_ns = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();
  std::cout << edges(mesh).size() << " edges. Elapsed: " << std::to_string(elapsed_ns) << " (ms)\n";

  return mesh;
}

// =================================================================================================
// =================================================================================================
// =================================================================================================

// always decimate meshes in this vector to avoid timing the copying of the meshes
template <typename Policy, typename TriangleMesh>
void time_mesh(const TriangleMesh mesh,
               std::ostream& out)
{
  std::vector<TriangleMesh> meshes (n_burns + n_samples, mesh);

  for(int i=0; i<n_burns; ++i)
    edge_collapse<Policy>(meshes[i]);

  std::chrono::time_point<std::chrono::steady_clock> start_time = std::chrono::steady_clock::now();
  for(int i=0; i<n_samples; ++i)
  {
    // measure time taken by the edge_collapse function over each mesh individually
    edge_collapse<Policy>(meshes[n_burns + i]);
  }
  std::chrono::time_point<std::chrono::steady_clock> end_time = std::chrono::steady_clock::now();
  auto elapsed_ns = std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count();

  out << "Policy: " << typeid(Policy).name() << "\n"
      << "Elapsed: " << elapsed_ns << " (ms)\n"
      << "Average: " << elapsed_ns / n_samples << " (ms)" << std::endl;
}

template <typename TriangleMesh>
void time_policy(const TriangleMesh& mesh,
                 std::ostream& out,
                 const std::string& policy)
{
  if(policy == "classic_plane")
    time_mesh<Classic_plane>(mesh, out);
  else if(policy == "classic_tri")
    time_mesh<Classic_tri>(mesh, out);
  else if(policy == "prob_plane")
    time_mesh<Prob_plane>(mesh, out);
  else if(policy == "prob_tri")
    time_mesh<Prob_tri>(mesh, out);
}

template <typename TriangleMesh>
void time_all_policies(const TriangleMesh& mesh,
                       std::ostream& out)
{
  std::cout << "   ### Time" << std::endl;

  time_mesh<Classic_plane>(mesh, out);
  time_mesh<Classic_tri>(mesh, out);
  time_mesh<Prob_plane>(mesh, out);
  time_mesh<Prob_tri>(mesh, out);
}

// =================================================================================================
// =================================================================================================
// =================================================================================================

template <typename Policy, typename TriangleMesh>
double hausdorff_error(const TriangleMesh& mesh,
                       double ratio = 0.2)
{
  // make a copy of the mesh so we can compare later
  TriangleMesh tmp = mesh;
  edge_collapse<Policy>(tmp, ratio);

  // arbitrary error bound
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(mesh);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));
  const double error_bound = 1e-4 * diag_length;

  return PMP::bounded_error_symmetric_Hausdorff_distance<CGAL::Parallel_if_available_tag>(mesh, tmp, error_bound);
}

// calculate approximate Hausdorff errors for all different policies at the same decimation ratio
template <typename TriangleMesh>
std::array<FT, 4> hausdorff_errors(const TriangleMesh& mesh,
                                   double ratio)
{
  std::array<FT, 4> ret { {0, 0, 0, 0} };

  ret[classic_plane_index] = hausdorff_error<Classic_plane>(mesh, ratio);
  ret[prob_plane_index] = hausdorff_error<Prob_plane>(mesh, ratio);
  ret[classic_tri_index] = hausdorff_error<Classic_tri>(mesh, ratio);
  ret[prob_tri_index] = hausdorff_error<Prob_tri>(mesh, ratio);

  return ret;
}

template <typename TriangleMesh, typename InputIt>
void hausdorff_errors(const TriangleMesh& mesh,
                      std::ostream& out,
                      InputIt begin, InputIt end)
{
  std::cout << "   ### Hausdorff error" << std::endl;

  for(InputIt it=begin; it!=end; ++it)
  {
    std::array<double, 4> errs = hausdorff_errors(mesh, *it);
    out << " -- Hausdorff error for collapse with ratio: " << *it << '\n';

    out << "classic plane: " << errs[classic_plane_index] << std::endl;
    out << "prob plane   : " << errs[prob_plane_index] << std::endl;
    out << "classic tri  : " << errs[classic_tri_index] << std::endl;
    out << "prob tri     : " << errs[prob_tri_index] << std::endl;
  }
}

// =================================================================================================
// =================================================================================================
// =================================================================================================

template <typename TriangleMesh, typename OutStream>
OutStream& write_aspect_ratios(const TriangleMesh& mesh,
                               OutStream& out)
{
  for(auto face : faces(mesh))
    out << std::to_string(PMP::face_aspect_ratio(face, mesh)) << '\n';

  return out;
}

template <typename TriangleMesh>
void gather_face_aspect_ratio(const TriangleMesh& mesh,
                              std::ostream& out)
{
  std::cout << "   ### Face aspect ratio" << std::endl;

  Surface_mesh cp = mesh;
  Surface_mesh pp = mesh;
  Surface_mesh ct = mesh;
  Surface_mesh pt = mesh;

  edge_collapse<Classic_plane>(cp);
  edge_collapse<Prob_plane>(pp);
  edge_collapse<Classic_tri>(ct);
  edge_collapse<Prob_tri>(pt);

  out << "Face aspect-ratio: classic plane\n";
  write_aspect_ratios(cp, out);
  out << "Face aspect-ratio: prob plane\n";
  write_aspect_ratios(pp, out);
  out << "Face aspect-ratio: classic triangle\n";
  write_aspect_ratios(ct, out);
  out << "Face aspect-ratio: prob triangle\n";
  write_aspect_ratios(pt, out);
}

template <typename TriangleMesh>
void run(const std::pair<TriangleMesh, std::string>& input)
{
  std::cout << std::endl;
  std::cout << " ====================================================================== " << std::endl;
  std::cout << " ====================================================================== " << std::endl;
  std::cout << "Run on " << input.second << std::endl;

  std::ostream& out = std::cout;

  time_all_policies(input.first, out);

//  std::array<double, 10> range {{ 0.15 }};
//  std::array<double, 10> range {{ 0.7, 0.6, 0.5, 0.45, 0.4, 0.35, 0.3, 0.25, 0.2, 0.15 }};
//  hausdorff_errors(input.first, out, range.begin(), range.end());

  gather_face_aspect_ratio(input.first, out);
}

int main(int argc, char** argv)
{
  std::vector<std::string> default_data = { "data/helmet.off",
                                            "data/femur.off",
                                            "data/oni.off",
                                            "data/genus1_null_edges.off" };

  std::vector<std::string> data;
  if(argc > 1)
    data = { std::string(argv[1]) };
  else
    data = default_data;

  std::vector<std::pair<Surface_mesh, std::string> > named_meshes = get_all_meshes(data);
  for(const auto& e : named_meshes)
    run(e);

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
