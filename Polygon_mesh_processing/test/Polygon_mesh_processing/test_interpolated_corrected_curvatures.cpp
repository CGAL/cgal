#include <CGAL/Polygon_mesh_processing/interpolated_corrected_curvatures.h>
#include <CGAL/Polygon_mesh_processing/random_perturbation.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/property_map.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/graph/graph_traits.hpp>

#include <iostream>
#include <unordered_map>

#define ABS_ERROR 1e-6

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epic_kernel;
typedef CGAL::Surface_mesh<Epic_kernel::Point_3> SMesh;
typedef CGAL::Polyhedron_3<Epic_kernel> Polyhedron;

struct Average_test_info {
  Epic_kernel::FT expansion_radius = -1;
  Epic_kernel::FT mean_curvature_avg;
  Epic_kernel::FT gaussian_curvature_avg;
  Epic_kernel::FT principal_curvature_avg;
  Epic_kernel::FT tolerance = 0.9;

  Average_test_info(
    Epic_kernel::FT mean_curvature_avg,
    Epic_kernel::FT gaussian_curvature_avg,
    Epic_kernel::FT principal_curvature_avg,
    Epic_kernel::FT expansion_radius = -1,
    Epic_kernel::FT tolerance = 0.9
  ) :
    expansion_radius(expansion_radius),
    mean_curvature_avg(mean_curvature_avg),
    gaussian_curvature_avg(gaussian_curvature_avg),
    principal_curvature_avg(principal_curvature_avg),
    tolerance(tolerance)
  {
  }

};

bool passes_comparison(Epic_kernel::FT result, Epic_kernel::FT expected, Epic_kernel::FT tolerance) {
  if (abs(expected) < ABS_ERROR && abs(result) < ABS_ERROR)
    return true; // expected 0, got 0
  else if (abs(expected) < ABS_ERROR)
    return false; // expected 0, got non-0

  return (std::min)(result, expected) / (std::max)(result, expected) > tolerance;
}

template <typename PolygonMesh>
void test_average_curvatures(std::string mesh_path,
  Average_test_info test_info,
  bool compare_single_vertex = false,
  int scale_factor_exponent = 0
) {
  typedef typename boost::graph_traits<PolygonMesh>::vertex_descriptor vertex_descriptor;

  PolygonMesh pmesh;
  const std::string filename = CGAL::data_file_path(mesh_path);

  if (!CGAL::IO::read_polygon_mesh(filename, pmesh) || faces(pmesh).size() == 0)
  {
    std::cerr << "Invalid input file." << std::endl;
  }

  // The following part is used to scale the given mesh and expected curvatures by a constant factor
  // this is used to test the stability of the implementation across different scales
  if (scale_factor_exponent) {
    Epic_kernel::FT factor = pow(10, scale_factor_exponent);

    test_info.expansion_radius *= factor;
    test_info.mean_curvature_avg /= factor;
    test_info.gaussian_curvature_avg /= factor * factor;
    test_info.principal_curvature_avg /= factor;

    auto vpm = get(CGAL::vertex_point, pmesh);

    for (vertex_descriptor vi : vertices(pmesh)) {
      Epic_kernel::Point_3 pi = get(vpm, vi);
      Epic_kernel::Point_3 pi_new(pi.x() * factor, pi.y() * factor, pi.z() * factor);
      put(vpm, vi, pi_new);
    }
  }


  typename boost::property_map<PolygonMesh, CGAL::dynamic_vertex_property_t<Epic_kernel::FT>>::type
    mean_curvature_map = get(CGAL::dynamic_vertex_property_t<Epic_kernel::FT>(), pmesh),
    gaussian_curvature_map = get(CGAL::dynamic_vertex_property_t<Epic_kernel::FT>(), pmesh);
  typename boost::property_map
    <PolygonMesh, CGAL::dynamic_vertex_property_t<PMP::Principal_curvatures_and_directions<Epic_kernel>>>::type
      principal_curvatures_and_directions_map =
        get(CGAL::dynamic_vertex_property_t<PMP::Principal_curvatures_and_directions<Epic_kernel>>(), pmesh);

  PMP::interpolated_corrected_curvatures(
    pmesh,
    CGAL::parameters::ball_radius(test_info.expansion_radius)
    .vertex_mean_curvature_map(mean_curvature_map)
    .vertex_Gaussian_curvature_map(gaussian_curvature_map)
    .vertex_principal_curvatures_and_directions_map(principal_curvatures_and_directions_map)
  );

  Epic_kernel::FT mean_curvature_avg = 0, gaussian_curvature_avg = 0, principal_curvature_avg = 0;

  for (vertex_descriptor v : vertices(pmesh)) {
    mean_curvature_avg += get(mean_curvature_map, v);
    gaussian_curvature_avg += get(gaussian_curvature_map, v);
    principal_curvature_avg += get(principal_curvatures_and_directions_map, v).min_curvature
      + get(principal_curvatures_and_directions_map, v).max_curvature;
  }

  mean_curvature_avg /= vertices(pmesh).size();
  gaussian_curvature_avg /= vertices(pmesh).size();
  principal_curvature_avg /= vertices(pmesh).size() * 2;

  // are average curvatures equal to expected?
  assert(passes_comparison(mean_curvature_avg, test_info.mean_curvature_avg, test_info.tolerance));
  assert(passes_comparison(gaussian_curvature_avg, test_info.gaussian_curvature_avg, test_info.tolerance));
  assert(passes_comparison(principal_curvature_avg, test_info.principal_curvature_avg, test_info.tolerance));

  Epic_kernel::FT new_mean_curvature_avg = 0, new_Gaussian_curvature_avg = 0, new_principal_curvature_avg = 0;

  for (vertex_descriptor v : vertices(pmesh)) {
    new_mean_curvature_avg += get(mean_curvature_map, v);
    new_Gaussian_curvature_avg += get(gaussian_curvature_map, v);
    new_principal_curvature_avg += get(principal_curvatures_and_directions_map, v).min_curvature
      + get(principal_curvatures_and_directions_map, v).max_curvature;
  }

  new_mean_curvature_avg /= vertices(pmesh).size();
  new_Gaussian_curvature_avg /= vertices(pmesh).size();
  new_principal_curvature_avg /= vertices(pmesh).size() * 2;

  // are average curvatures computed from interpolated_corrected_curvatures()
  // equal to average curvatures each computed on its own?
  assert(passes_comparison(mean_curvature_avg, new_mean_curvature_avg, 0.99));
  assert(passes_comparison(gaussian_curvature_avg, new_Gaussian_curvature_avg, 0.99));
  assert(passes_comparison(principal_curvature_avg, new_principal_curvature_avg, 0.99));

  if (compare_single_vertex) {
    // computing curvatures together from interpolated_corrected_curvatures()

    Epic_kernel::FT single_vertex_mean_curvature_avg = 0,
      single_vertex_Gaussian_curvature_avg = 0,
      single_vertex_principal_curvature_avg = 0;

    Epic_kernel::FT h, g;
    PMP::Principal_curvatures_and_directions<Epic_kernel> p;

    for (vertex_descriptor v : vertices(pmesh)) {
      PMP::interpolated_corrected_curvatures(
        v,
        pmesh,
        CGAL::parameters::vertex_Gaussian_curvature(std::ref(g))
        .vertex_mean_curvature(std::ref(h))
        .vertex_principal_curvatures_and_directions(std::ref(p))
        .ball_radius(test_info.expansion_radius)
      );

      single_vertex_mean_curvature_avg += h;
      single_vertex_Gaussian_curvature_avg += g;
      single_vertex_principal_curvature_avg += p.min_curvature + p.max_curvature;
    }

    single_vertex_mean_curvature_avg /= vertices(pmesh).size();
    single_vertex_Gaussian_curvature_avg /= vertices(pmesh).size();
    single_vertex_principal_curvature_avg /= vertices(pmesh).size() * 2;

    assert(passes_comparison(mean_curvature_avg, single_vertex_mean_curvature_avg, 0.99));
    assert(passes_comparison(gaussian_curvature_avg, single_vertex_Gaussian_curvature_avg, 0.99));
    assert(passes_comparison(principal_curvature_avg, single_vertex_principal_curvature_avg, 0.99));
  }

}

int main()
{
  // testing on a simple sphere(r = 0.5), on both Polyhedron & SurfaceMesh:
  // For this mesh, ina addition to the whole mesh functions, we also compare against the single vertex
  // curvature functions to make sure the produce the same results
  // Expected: Mean Curvature = 2, Gaussian Curvature = 4, Principal Curvatures = 2 & 2 so 2 on avg.
  test_average_curvatures<Polyhedron>("meshes/sphere.off", Average_test_info(2, 4, 2), true);
  test_average_curvatures<SMesh>("meshes/sphere.off", Average_test_info(2, 4, 2), true);

  // Same mesh but with specified expansion radii of 0 and 0.25 (half radius of sphere)
  test_average_curvatures<SMesh>("meshes/sphere.off", Average_test_info(2, 4, 2, 0), true);
  test_average_curvatures<SMesh>("meshes/sphere.off", Average_test_info(2, 4, 2, 0.25), true);

  // testing on a simple sphere(r = 10), on both Polyhedron & SurfaceMesh:
  // Expected: Mean Curvature = 0.1, Gaussian Curvature = 0.01, Principal Curvatures = 0.1 & 0.1 so 0.1 on avg.
  test_average_curvatures<Polyhedron>("meshes/sphere966.off", Average_test_info(0.1, 0.01, 0.1));
  test_average_curvatures<SMesh>("meshes/sphere966.off", Average_test_info(0.1, 0.01, 0.1));

  // Same mesh but with specified expansion radii of 0 and 5 (half radius of sphere)
  test_average_curvatures<SMesh>("meshes/sphere966.off", Average_test_info(0.1, 0.01, 0.1, 0));
  test_average_curvatures<SMesh>("meshes/sphere966.off", Average_test_info(0.1, 0.01, 0.1, 5));

  // testing on a simple half cylinder(r = 1), on both Polyhedron & SurfaceMesh:
  // Expected: Mean Curvature = 0.5, Gaussian Curvature = 0, Principal Curvatures = 0 & 1 so 0.5 on avg.
  test_average_curvatures<Polyhedron>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5));
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5));

  // Same mesh but with specified expansion radii of 0 and 0.5 (half radius of cylinder)
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0));
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0.5));

  // Same tests as last one, but with a scaling on the mesh with different values to check for scale stability
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0), false, -6);
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0.5), false, -6);

  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0), false, -3);
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0.5), false, -3);

  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0), false, -1);
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0.5), false, -1);

  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0), false, 1);
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0.5), false, 1);

  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0), false, 3);
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0.5), false, 3);

  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0), false, 6);
  test_average_curvatures<SMesh>("meshes/cylinder.off", Average_test_info(0.5, 0, 0.5, 0.5), false, 6);
}
