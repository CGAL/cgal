#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>

#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
typedef CGAL::Surface_mesh<K::Point_3>                        Mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main(int argc, char* argv[])
{
  const std::string filename1 = (argc > 1) ? argv[1] : CGAL::data_file_path("meshes/blobby.off");
  const std::string filename2 = (argc > 2) ? argv[2] : CGAL::data_file_path("meshes/eight.off");

  Mesh mesh1, mesh2;
  if(!PMP::IO::read_polygon_mesh(filename1, mesh1) || !PMP::IO::read_polygon_mesh(filename2, mesh2))
  {
    std::cerr << "Invalid input." << std::endl;
    return 1;
  }


  PMP::Corefinement::Visitor_for_non_manifold_output<K, Mesh> visitor(mesh1, mesh2);

  std::array<Mesh,4> out_meshes;
  std::array<std::optional<Mesh*>, 4> output = {&out_meshes[0], &out_meshes[1], &out_meshes[2], &out_meshes[3]};

  std::array<bool,4> res = PMP::corefine_and_compute_boolean_operations(mesh1, mesh2, output, CGAL::parameters::visitor(visitor));

  // UNION
  {
  std::vector<K::Point_3> points;
  std::vector< std::array<std::size_t, 3> > polygons;
  visitor.extract_union(points, polygons);
  CGAL::IO::write_polygon_soup("union.off", points, polygons, CGAL::parameters::stream_precision(17));
  assert(res[PMP::Corefinement::UNION] == PMP::is_polygon_soup_a_polygon_mesh(polygons));
  }
  // INTERSECTION
  {
  std::vector<K::Point_3> points;
  std::vector< std::array<std::size_t, 3> > polygons;
  visitor.extract_intersection(points, polygons);
  CGAL::IO::write_polygon_soup("inter.off", points, polygons, CGAL::parameters::stream_precision(17));
  assert(res[PMP::Corefinement::INTERSECTION] == PMP::is_polygon_soup_a_polygon_mesh(polygons));
  }
  // TM1_MINUS_TM2
  {
  std::vector<K::Point_3> points;
  std::vector< std::array<std::size_t, 3> > polygons;
  visitor.extract_tm1_minus_tm2(points, polygons);
  CGAL::IO::write_polygon_soup("tm1_minus_tm2.off", points, polygons, CGAL::parameters::stream_precision(17));
  assert(res[PMP::Corefinement::TM1_MINUS_TM2] == PMP::is_polygon_soup_a_polygon_mesh(polygons));
  }
  // TM2_MINUS_TM1
  {
  std::vector<K::Point_3> points;
  std::vector< std::array<std::size_t, 3> > polygons;
  visitor.extract_tm2_minus_tm1(points, polygons);
  CGAL::IO::write_polygon_soup("tm2_minus_tm1.off", points, polygons, CGAL::parameters::stream_precision(17));
  assert(res[PMP::Corefinement::TM2_MINUS_TM1] == PMP::is_polygon_soup_a_polygon_mesh(polygons));
  }

  //TODO common faces

  return 0;
}
