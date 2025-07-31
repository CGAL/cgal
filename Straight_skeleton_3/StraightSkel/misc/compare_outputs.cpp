#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/distance.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/shape_predicates.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <iostream>
#include <fstream>

namespace PMP = ::CGAL::Polygon_mesh_processing;

using EPECK = CGAL::Exact_predicates_inexact_constructions_kernel;
using K = EPECK;
using FT = K::FT;
using Point = K::Point_3;
using Vector = K::Vector_3;

using Mesh = CGAL::Surface_mesh<Point>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

template <typename Mesh>
bool has_degenerated_faces(const Mesh& mesh)
{
  for(auto f : faces(mesh))
    if(PMP::is_degenerate_triangle_face(f, mesh))
      return true;

  return false;
}

bool compare(const Mesh& ours,
             const Mesh& theirs)
{
  bool are_similar = true;

  // volume
  const FT vol_ours = PMP::volume(ours);
  const FT vol_theirs = PMP::volume(theirs);
  std::cout << "vol_ours = " << vol_ours << std::endl;
  std::cout << "vol_theirs = " << vol_theirs << std::endl;
  std::cout << "volume ratio = " << vol_ours / vol_theirs << std::endl;

  if(CGAL::abs(vol_ours - vol_theirs) > 0.001 * vol_ours) {
    std::cerr << "Error: difference in volume" << std::endl;
    if(CGAL::abs(vol_ours - vol_theirs) > 0.01 * vol_ours) {
      std::cerr << "Error: significant difference in volume" << std::endl;
      if(CGAL::abs(vol_ours - vol_theirs) > 0.1 * vol_ours) {
        std::cerr << "Error: REALLY significant difference in volume" << std::endl;
      }
    }
    are_similar = false;
  }

  // area
  const FT area_ours = PMP::area(ours);
  const FT area_theirs = PMP::area(theirs);
  std::cout << "area_ours = " << area_ours << std::endl;
  std::cout << "area_theirs = " << area_theirs << std::endl;
  std::cout << "area ratio = " << area_ours / area_theirs << std::endl;

  if(CGAL::abs(area_ours - area_theirs) > 0.01 * area_ours) {
    std::cerr << "Error: significant difference in area" << std::endl;
    are_similar = false;
  }

  // Hausdorff
  CGAL::Bbox_3 bbox = CGAL::Polygon_mesh_processing::bbox(ours);
  const double diag_length = std::sqrt(CGAL::square(bbox.xmax() - bbox.xmin()) +
                                       CGAL::square(bbox.ymax() - bbox.ymin()) +
                                       CGAL::square(bbox.zmax() - bbox.zmin()));

  FT hd = PMP::approximate_symmetric_Hausdorff_distance<CGAL::Parallel_if_available_tag>(
                ours, theirs,
                CGAL::parameters::number_of_points_on_edges(10)
                                 .number_of_points_on_faces(100)
                                 .random_seed(0));
  std::cout << "symmetric Hausdorff distance: " << hd << std::endl;
  std::cout << "symmetric Hausdorff ratio: " << hd / diag_length << std::endl;
  if(hd > 0.001 * diag_length) {
    std::cerr << "Error: large Hausdorff distance" << std::endl;
    are_similar = false;
  }

  return are_similar;
}

int main(int argc, char** argv)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  if(argc < 2)
  {
    std::cerr << "Usage: "
              << argv[0] << "\n"
              << "\tours_filename\n"
              << "\ttheirs_filename\n"
              << std::endl;
    std::cerr << "Input format: any format readable with CGAL::IO::read_polygon_mesh()" << std::endl;
    return EXIT_FAILURE;
  }

  const char* ours_filename = argv[1];
  const char* theirs_filename = argv[2];

  Mesh ours;
  if(!CGAL::IO::read_polygon_mesh(ours_filename, ours)) {
    std::cerr << "Error: failed to read input (ours)" << std::endl;
    return EXIT_FAILURE;
  }

  Mesh theirs;
  if(!CGAL::IO::read_polygon_mesh(theirs_filename, theirs)) {
    std::cerr << "Error: failed to read input (theirs)" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Ours: " << num_vertices(ours) << " NV " << num_faces(theirs) << " NF" << std::endl;
  std::cout << "Theirs: " << num_vertices(theirs) << " NV " << num_faces(theirs) << " NF" << std::endl;

  const bool closed_ours = CGAL::is_closed(ours);
  const bool closed_theirs = CGAL::is_closed(theirs);
  const bool has_degen_ours = has_degenerated_faces(ours);
  const bool has_degen_theirs = has_degenerated_faces(theirs);
  const bool SI_ours = PMP::does_self_intersect(ours);
  const bool SI_theirs = PMP::does_self_intersect(theirs);

  std::cout << "Health checks (ours) " << closed_ours << " " << !has_degen_ours << " " << !SI_ours << std::endl;
  std::cout << "Health checks (theirs) " << closed_theirs << " " << !has_degen_theirs << " " << !SI_theirs << std::endl;

  if(!closed_ours || has_degen_ours || SI_ours) {
    std::cerr << "Error: 'ours' is BAD" << std::endl;
    return EXIT_FAILURE;
  }

  if(!closed_theirs || has_degen_theirs || SI_theirs) {
    std::cerr << "Error: 'theirs' is BAD" << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "Valid inputs" << std::endl;

  if(!compare(ours, theirs)) {
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;
}
