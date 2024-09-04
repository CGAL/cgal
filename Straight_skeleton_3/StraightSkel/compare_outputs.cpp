#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/self_intersections.h>

#include <iostream>
#include <fstream>

namespace PMP = ::CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using FT = K::FT;
using Point = K::Point_3;
using Vector = K::Vector_3;

using Mesh = CGAL::Surface_mesh<Point>;
using face_descriptor = boost::graph_traits<Mesh>::face_descriptor;

bool compare(const Mesh& ours,
             const Mesh& theirs)
{
  std::cout << "Ours: " << num_vertices(ours) << " NV" << num_faces(ours) << " NF" << std::endl;
  std::cout << "Theirs: " << num_vertices(theirs) << " NV" << num_faces(theirs) << " NF" << std::endl;

  const bool SI_ours = PMP::does_self_intersect(ours);
  const bool SI_theirs = PMP::does_self_intersect(theirs);
  const bool closed_ours = CGAL::is_closed(ours);
  const bool closed_theirs = CGAL::is_closed(theirs);

  if(SI_theirs || !closed_theirs) {
    std::cerr << "Error: result to compare to is BAD (self-intersections / not closed)" << std::endl;
    return false;
  }

  if(SI_ours || !closed_ours) {
    std::cerr << "Error: our result is BAD (self-intersections / not closed)" << std::endl;
    return false;
  }

  // volume
  const FT vol_ours = PMP::volume(ours);
  const FT vol_theirs = PMP::volume(theirs);
  std::cout << "vol_ours = " << vol_ours << std::endl;
  std::cout << "vol_theirs = " << vol_theirs << std::endl;

  if(CGAL::abs(vol_ours - vol_theirs) > 0.001 * vol_ours) {
    std::cerr << "Error: significant difference in volume" << std::endl;
    return false;
  }

  // area
  const FT area_ours = PMP::area(ours);
  const FT area_theirs = PMP::area(theirs);
  std::cout << "area_ours = " << area_ours << std::endl;
  std::cout << "area_theirs = " << area_theirs << std::endl;

  if(CGAL::abs(area_ours - area_theirs) > 0.01 * area_ours) {
    std::cerr << "Error: significant difference in area" << std::endl;
    return false;
  }

  // @todo add:
  // - Hausdorff

  return true;
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

  if(!compare(ours, theirs))
    return EXIT_FAILURE;

  return EXIT_SUCCESS;
}
