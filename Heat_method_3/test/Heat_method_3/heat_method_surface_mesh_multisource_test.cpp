// Multi-source accuracy test for Surface_mesh_geodesic_distances_3.
//
// Mesh: data/flat_star_disk_irregular_rim.off — a FLAT star-shaped disk (z = 0)
// with irregular, variable-density boundary sampling; every boundary vertex is
// used as a source (508 sources, 3000 vertices). On a flat, simply covered mesh
// the exact geodesic distance to the source set equals, in the near-boundary
// band, the euclidean distance to the nearest source (verified to ~1e-15
// against exact polyhedral geodesics when the fixture was generated), and the
// euclidean distance is a hard LOWER bound everywhere. The test therefore
// needs no external oracle.
//
// The irregular rim sampling makes the Poisson potential phi spread across the
// source set, which exposes the source-set normalization of
// value_at_source_set(): mapping phi through  min_s |phi(i) - phi(s)|  (the
// distance from the VALUE phi(i) to the SET of source values) folds every
// vertex whose phi lies inside the source-value spread toward 0 and shifts the
// far field by max_s phi(s) instead of a single constant. Measured on this
// mesh (band = 2h..6h of the euclidean distance to the source set, h = mean
// edge length): mean d/d_true = 0.54 and 83% of the band vertices read BELOW
// the euclidean lower bound — impossible for a distance. With the mean-shift
// normalization the ratio is ~1.1 and lower-bound violations are marginal.

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Heat_method_3/Surface_mesh_geodesic_distances_3.h>

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

typedef CGAL::Simple_cartesian<double>                         Kernel;
typedef Kernel::Point_3                                        Point_3;
typedef CGAL::Surface_mesh<Point_3>                            Surface_mesh;
typedef boost::graph_traits<Surface_mesh>::vertex_descriptor   vertex_descriptor;
typedef Surface_mesh::Property_map<vertex_descriptor, double>  Vertex_distance_map;

template <typename Mode>
void test_mode(const char* mode_name, Surface_mesh& sm)
{
  typedef CGAL::Heat_method_3::Surface_mesh_geodesic_distances_3<Surface_mesh, Mode> Heat_method;

  Vertex_distance_map vertex_distance =
      sm.add_property_map<vertex_descriptor, double>("v:distance", 0).first;

  // the source set: every boundary vertex
  std::vector<vertex_descriptor> sources;
  for(vertex_descriptor vd : vertices(sm)) {
    if(sm.is_border(vd)) {
      sources.push_back(vd);
    }
  }
  assert(sources.size() == 508); // fixture integrity

  Heat_method hm(sm);
  hm.add_sources(sources);
  hm.estimate_geodesic_distances(vertex_distance);

  // mean edge length (band scale)
  double h = 0.;
  for(auto ed : edges(sm)) {
    vertex_descriptor s = source(halfedge(ed, sm), sm);
    vertex_descriptor t = target(halfedge(ed, sm), sm);
    h += std::sqrt(CGAL::squared_distance(sm.point(s), sm.point(t)));
  }
  h /= static_cast<double>(num_edges(sm));

  // euclidean distance to the nearest source: exact truth in the near band on
  // this flat mesh, and a hard lower bound everywhere
  std::size_t n_band = 0, n_below_bound = 0;
  double ratio_sum = 0.;
  for(vertex_descriptor vd : vertices(sm)) {
    double eu2 = (std::numeric_limits<double>::max)();
    for(vertex_descriptor sv : sources) {
      eu2 = (std::min)(eu2, CGAL::to_double(CGAL::squared_distance(sm.point(vd), sm.point(sv))));
    }
    const double eu = std::sqrt(eu2);
    if(eu <= 2. * h || eu >= 6. * h) {
      continue;
    }
    const double d = get(vertex_distance, vd);
    n_band++;
    ratio_sum += d / eu;
    if(d < 0.75 * eu) {
      n_below_bound++;
    }
  }
  assert(n_band > 1000);

  const double mean_ratio = ratio_sum / static_cast<double>(n_band);
  const double below_fraction = static_cast<double>(n_below_bound) / static_cast<double>(n_band);
  std::cout << mode_name << ": band vertices = " << n_band
            << ", mean d/d_true = " << mean_ratio
            << ", below euclidean lower bound = " << 100. * below_fraction << "%"
            << std::endl;

  // before the fix: mean_ratio ~ 0.54 and ~83% below the lower bound
  assert(mean_ratio > 0.85);
  assert(mean_ratio < 1.30);
  assert(below_fraction < 0.10);
}

int main(int argc, char* argv[])
{
  Surface_mesh sm;
  const std::string filename = (argc > 1) ? argv[1] : "data/flat_star_disk_irregular_rim.off";
  std::ifstream in(filename);
  in >> sm;
  if(!in || sm.number_of_vertices() == 0) {
    std::cerr << "could not read " << filename << std::endl;
    return EXIT_FAILURE;
  }
  assert(sm.number_of_vertices() == 3000);

  test_mode<CGAL::Heat_method_3::Direct>("Direct", sm);
  test_mode<CGAL::Heat_method_3::Intrinsic_Delaunay>("Intrinsic_Delaunay", sm);

  std::cout << "done" << std::endl;
  return EXIT_SUCCESS;
}
