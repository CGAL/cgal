#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Edge_count_ratio_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Bounded_normal_change_placement.h>
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/GarlandHeckbert_policies.h>

#include <CGAL/Polygon_mesh_processing/autorefinement.h>

#include <CGAL/Polygon_mesh_processing/manifoldness.h>
#include <CGAL/Polygon_mesh_processing/distance.h>

#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <CGAL/IO/polygon_mesh_io.h>

#include <chrono>
#include <fstream>
#include <iostream>
#include <filesystem>
#include <vector>

#define TAG CGAL::Parallel_if_available_tag

typedef CGAL::Simple_cartesian<double>                          Kernel;
typedef Kernel::FT                                              FT;
typedef Kernel::Point_3                                         Point_3;
typedef CGAL::Surface_mesh<Point_3>                             Surface_mesh;

namespace SMS = CGAL::Surface_mesh_simplification;


namespace PMP = CGAL::Polygon_mesh_processing;

typedef SMS::GarlandHeckbert_plane_policies<Surface_mesh, Kernel>                  Classic_plane;
typedef SMS::GarlandHeckbert_probabilistic_plane_policies<Surface_mesh, Kernel>    Prob_plane;
typedef SMS::GarlandHeckbert_triangle_policies<Surface_mesh, Kernel>               Classic_tri;
typedef SMS::GarlandHeckbert_probabilistic_triangle_policies<Surface_mesh, Kernel> Prob_tri;
typedef SMS::GarlandHeckbert_policies<Surface_mesh, Kernel>                        Classic_plane_plus_line;

// Old policies composed with line policies
typedef SMS::internal::GarlandHeckbert_line_policies<Surface_mesh, Kernel> Line_quadric;
typedef SMS::internal::GarlandHeckbert_composed_policies<Surface_mesh, Kernel, Prob_plane, Line_quadric> Proba_plane_plus_line;
typedef SMS::internal::GarlandHeckbert_composed_policies<Surface_mesh, Kernel, Classic_tri, Line_quadric> Classic_tri_plus_line;
typedef SMS::internal::GarlandHeckbert_composed_policies<Surface_mesh, Kernel, Prob_tri, Line_quadric> Proba_tri_plus_line;

double mean_aspect_ratio(Surface_mesh& mesh){
  double total_aspect_ratio=0;
  double total_area=0;
  for(auto f: mesh.faces()){
    double a=PMP::face_area(f, mesh);
    total_aspect_ratio+=a*CGAL::to_double(PMP::face_aspect_ratio(f, mesh));
    total_area+=a;
  }
  return total_aspect_ratio/total_area;
}

template <typename GHPolicies>
double collapse_gh(Surface_mesh& mesh,
                 const double ratio,
                 GHPolicies gh_policies)
{
  std::chrono::steady_clock::time_point start_time = std::chrono::steady_clock::now();

  SMS::Edge_count_ratio_stop_predicate<Surface_mesh> stop(ratio);

  // Garland&Heckbert simplification policies

  typedef typename GHPolicies::Get_cost                                        GH_cost;
  typedef typename GHPolicies::Get_placement                                   GH_placement;
  typedef SMS::Bounded_normal_change_placement<GH_placement>                   Bounded_GH_placement;

  // GHPolicies gh_policies(mesh);
  const GH_cost& gh_cost = gh_policies.get_cost();
  const GH_placement& gh_placement = gh_policies.get_placement();
  Bounded_GH_placement placement(gh_placement);

  int r = SMS::edge_collapse(mesh, stop,
                             CGAL::parameters::get_cost(gh_cost)
                                              .get_placement(placement));

  std::chrono::steady_clock::time_point end_time = std::chrono::steady_clock::now();

  return (std::chrono::duration_cast<std::chrono::milliseconds>(end_time - start_time).count()) / 1000.;
}

std::string getFileNameWithoutExtension(const std::string& filePath) {
    std::filesystem::path path(filePath);
    return path.stem().string();  // 'stem' gets the filename without extension
}

// Usage:
// ./command [input] [ratio] [policy] [output]
// policy can be "cp" (classic plane), "ct" (classic triangle), "pp" (probabilistic plane), "pt" (probabilistic triangle)
int main(int argc, char** argv)
{
  Surface_mesh mesh;
  const std::string filename = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/cube-meshed.off");
  PMP::duplicate_non_manifold_vertices(mesh);
  PMP::autorefine(mesh);

  if(!CGAL::IO::read_polygon_mesh(filename, mesh))
  {
    std::cerr << "Failed to read input mesh: " << filename << std::endl;
    return EXIT_FAILURE;
  }

  if(!CGAL::is_triangle_mesh(mesh))
  {
    std::cerr << "Input geometry is not triangulated." << std::endl;
    return EXIT_FAILURE;
  }

  const double ratio = (argc > 2) ? std::stod(argv[2]) : 0.2;
  std::cout << "\n\nPolicy Approx_Hausdorff Dist_points_to_input mean_aspect_ratio running_time" << std::endl;
  auto collapse=[&](std::string policy){
    Surface_mesh tmesh(mesh);
    double time;
    if(policy == "Classic_plane")
      time=collapse_gh(tmesh, ratio, Classic_plane(tmesh));
    else if(policy == "Classic_triangle")
      time=collapse_gh(tmesh, ratio, Classic_tri(tmesh));
    else if(policy == "Prob_plane")
      time=collapse_gh(tmesh, ratio, Prob_plane(tmesh));
    else if(policy == "Prob_triangle")
      time=collapse_gh(tmesh, ratio, Prob_tri(tmesh));
    else if(policy == "Plane_plus_line_0.1")
      time=collapse_gh(tmesh, ratio, Classic_plane_plus_line(tmesh, 100, 0.1));
    else if(policy == "Plane_plus_line_0.01")
      time=collapse_gh(tmesh, ratio, Classic_plane_plus_line(tmesh, 100, 0.01));
    else if(policy == "Plane_plus_line_0.001")
      time=collapse_gh(tmesh, ratio, Classic_plane_plus_line(tmesh, 100, 0.001));
    else if(policy == "Classic_tri_line")
      time=collapse_gh(tmesh, ratio, Classic_tri_plus_line(tmesh, 100));
    else if(policy == "Proba_plane_line")
      time=collapse_gh(tmesh, ratio, Proba_plane_plus_line(tmesh, Prob_plane(tmesh), Line_quadric(tmesh), 100));
    else if(policy == "Proba_tri_line")
      time=collapse_gh(tmesh, ratio, Proba_tri_plus_line(tmesh, Prob_tri(tmesh), Line_quadric(tmesh), 100));

    std::cout << policy << " " << PMP::approximate_Hausdorff_distance<TAG>(tmesh, mesh)
                        << " " << PMP::max_distance_to_triangle_mesh<TAG>(tmesh.points(),mesh)
                        << " " << mean_aspect_ratio(tmesh)
                        << " " << time << std::endl;
    CGAL::IO::write_polygon_mesh("out_"+policy+".off", tmesh, CGAL::parameters::stream_precision(17));
  };

  // collapse("Classic_plane");
  // collapse("Classic_triangle");
  collapse("Prob_plane");
  collapse("Prob_triangle");
  // collapse("Plane_plus_line_0.1");
  collapse("Plane_plus_line_0.01");
  collapse("Classic_tri_line");
  collapse("Proba_plane_line");
  collapse("Proba_tri_line");
  // collapse("Plane_plus_line_0.001");

  return EXIT_SUCCESS;
}
