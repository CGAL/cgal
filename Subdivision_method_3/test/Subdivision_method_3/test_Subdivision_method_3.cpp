#include <CGAL/subdivision_method_3.h>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <iostream>
#include <fstream>
#include <cassert>

namespace PMP = ::CGAL::Polygon_mesh_processing;
namespace SM = ::CGAL::Subdivision_method_3;
namespace params = SM::parameters;

#define TESTMESH_QUAD      CGAL::data_file_path("meshes/corner.off")
#define TESTMESH_QUAD_OPEN CGAL::data_file_path("meshes/corner_with_hole.off")

#define TESTMESH_TRI       CGAL::data_file_path("meshes/quint_tris.off")
#define TESTMESH_TRI_OPEN  CGAL::data_file_path("meshes/nefertiti.off")

template <typename Mesh>
void test_Subdivision_surface(const int depth = 3,
                              const bool do_not_modify_geometry = false)
{
  auto check_volume_unchanged = [](Mesh& P)
  {
    if (!is_triangle_mesh(P))
      PMP::triangulate_faces(P);

    double area_before = PMP::area(P);
    double area_after = PMP::area(P);
    assert(std::abs(area_before - area_after) < 1e-6);

    if (is_closed(P)) {
      double volume_before = PMP::volume(P);
      double volume_after = PMP::volume(P);
      assert(std::abs(volume_before - volume_after) < 1e-6);
    }
  };

  auto write_output = [](const Mesh& P, const std::string& test_name, int depth, bool do_not_modify_geometry) {
    std::string filename = "out_" + test_name + "_" + std::to_string(depth)
                             + (do_not_modify_geometry ? "_fixed." : ".") + "off";
    if (!CGAL::IO::write_polygon_mesh(filename, P, CGAL::parameters::stream_precision(17))) {
      std::cerr << "Error: Could not write to file " << filename << std::endl;
    }
  };

  // Test Catmull-Clark subdivision on quad mesh
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_QUAD, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_QUAD << std::endl;
      return;
    }

    SM::CatmullClark_subdivision(P, params::vertex_point_map(get(CGAL::vertex_point, P))
                                           .number_of_iterations(depth));
    assert(CGAL::is_valid_polygon_mesh(P));

    write_output(P, "CatmullClark_quad", depth, false /*do not modify*/);
  }

  // Test Catmull-Clark subdivision on 'opened' quad mesh
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_QUAD_OPEN, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_QUAD_OPEN << std::endl;
      return;
    }

    SM::CatmullClark_subdivision(P, params::vertex_point_map(get(CGAL::vertex_point, P))
                                           .number_of_iterations(depth)
                                           .do_not_modify_geometry(do_not_modify_geometry));
    assert(CGAL::is_valid_polygon_mesh(P));
    if (do_not_modify_geometry) {
      check_volume_unchanged(P);
    }

    write_output(P, "CatmullClark_quad_open", depth, do_not_modify_geometry);
  }

  // Test Catmull-Clark subdivision on 'opened' tri mesh
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_TRI_OPEN, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_TRI_OPEN << std::endl;
      return;
    }

    SM::CatmullClark_subdivision(P, params::vertex_point_map(get(CGAL::vertex_point, P))
                                           .number_of_iterations(depth)
                                           .do_not_modify_geometry(do_not_modify_geometry));
    assert(CGAL::is_valid_polygon_mesh(P));
    if (do_not_modify_geometry) {
      check_volume_unchanged(P);
    }

    write_output(P, "CatmullClark_tri_open", depth, do_not_modify_geometry);
  }

  // Test Loop subdivision on tri mesh
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_TRI, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_TRI << std::endl;
      return;
    }

    SM::Loop_subdivision(P, params::vertex_point_map(get(CGAL::vertex_point, P))
                                   .number_of_iterations(depth)
                                   .do_not_modify_geometry(do_not_modify_geometry));
    assert(CGAL::is_valid_polygon_mesh(P));
    if (do_not_modify_geometry) {
      check_volume_unchanged(P);
    }

    write_output(P, "Loop_tri", depth, do_not_modify_geometry);
  }

  // Test Loop subdivision on 'opened' tri mesh
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_TRI_OPEN, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_TRI_OPEN << std::endl;
      return;
    }

    SM::Loop_subdivision(P, params::vertex_point_map(get(CGAL::vertex_point, P))
                                   .number_of_iterations(depth)
                                   .do_not_modify_geometry(do_not_modify_geometry));
    assert(CGAL::is_valid_polygon_mesh(P));
    if (do_not_modify_geometry) {
      check_volume_unchanged(P);
    }

    write_output(P, "Loop_tri_open", depth, do_not_modify_geometry);
  }

  // Test Doo-Sabin subdivision on general mesh
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_TRI_OPEN, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_TRI_OPEN << std::endl;
      return;
    }

    SM::DooSabin_subdivision(P, params::vertex_point_map(get(CGAL::vertex_point, P))
                                       .number_of_iterations(depth));
    assert(CGAL::is_valid_polygon_mesh(P));

    write_output(P, "DooSabin_general", depth, false);
  }

  // Test Sqrt-3 subdivision on tri mesh
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_TRI, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_TRI << std::endl;
      return;
    }

    SM::Sqrt3_subdivision(P, params::vertex_point_map(get(CGAL::vertex_point, P))
                                    .number_of_iterations(depth));
    assert(CGAL::is_valid_polygon_mesh(P));
    check_volume_unchanged(P);

    write_output(P, "Sqrt3_tri", depth, false);
  }

  // Test Sqrt-3 subdivision on 'opened' tri mesh
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_TRI_OPEN, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_TRI_OPEN << std::endl;
      return;
    }

    SM::Sqrt3_subdivision(P, params::vertex_point_map(get(CGAL::vertex_point, P))
                                    .number_of_iterations(depth));
    assert(CGAL::is_valid_polygon_mesh(P));

    write_output(P, "Sqrt3_tri_open", depth, false);
  }

  // Test Sqrt-3 subdivision on 'opened' tri mesh with external property map
  {
    Mesh P;
    if (!CGAL::IO::read_polygon_mesh(TESTMESH_TRI_OPEN, P) || is_empty(P)) {
      std::cerr << "Error: failed to read or empty mesh: " << TESTMESH_TRI_OPEN << std::endl;
      return;
    }

    typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPM;
    typedef typename boost::property_traits<VPM>::value_type Point;
    typedef typename boost::property_traits<VPM>::reference PointRef;
    typedef typename CGAL::Kernel_traits<Point>::type::Vector_3 Vector;

    std::unordered_map<typename boost::graph_traits<Mesh>::vertex_descriptor, Point> um;
    boost::associative_property_map<std::unordered_map<typename boost::graph_traits<Mesh>::vertex_descriptor, Point>> apm(um);

    VPM vpm = get(CGAL::vertex_point, P);

    // Assign arbitrary new coordinates
    for (auto vd : vertices(P)) {
      PointRef pt = get(vpm, vd);
      Vector v = pt - Point(0., 0., -3.);
      put(apm, vd, pt + 0.5 * v);
    }

    SM::Sqrt3_subdivision(P, params::vertex_point_map(apm)
                                    .number_of_iterations(depth));
    assert(CGAL::is_valid_polygon_mesh(P));

    write_output(P, "Sqrt3_tri_open_external_map", depth, false);
  }
}

int main(int argc, char* argv[])
{
  const int depth = (argc > 1) ? std::stoi(argv[1]) : 3;

  typedef CGAL::Simple_cartesian<double> Kernel;
  typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
  typedef CGAL::Surface_mesh<Kernel::Point_3> SurfaceMesh;

  test_Subdivision_surface<Polyhedron>(depth);
  test_Subdivision_surface<SurfaceMesh>(depth);
  test_Subdivision_surface<SurfaceMesh>(depth, true /*do_not_modify_geometry*/);

  std::cerr << "Done" << std::endl;

  return 0;
}
