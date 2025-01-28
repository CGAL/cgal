#include <iostream>
#include <fstream>
#include <ctime>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Bbox_3.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>

#include <CGAL/Variational_shape_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;

typedef CGAL::Surface_mesh<Kernel::Point_3> Mesh;
typedef boost::graph_traits<Mesh>::face_descriptor face_descriptor;
typedef boost::graph_traits<Mesh>::vertex_descriptor vertex_descriptor;
typedef boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;

typedef boost::property_map<Mesh, boost::vertex_point_t>::type Vertex_point_map;
typedef Mesh::Property_map<face_descriptor, std::size_t> Face_proxy_map;

typedef CGAL::Variational_shape_approximation<Mesh, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;
typedef L21_approx::Proxy Plane_proxies;

#define CGAL_VSA_TEST_TOLERANCE 1e-8

bool check_strict_ordering(const std::vector<FT> &error)
{
  if (error.empty()) {
    std::cout << "Empty error sequence." << std::endl;
    return true;
  }

  FT pre = error.front();
  for(const FT e : error) {
    if (pre < e)
      return false;
  }

  return true;
}

/**
 * This file tests the teleportation functionality on a plane-sphere shape.
 * We seed randomly first, then verify that teleporting all (most) planes
 * from planar part to the spherical one and lower the error.
 */
int main()
{
  Mesh mesh;
  std::ifstream input("data/plane-sphere-high.off");
  if (!input || !(input >> mesh) || !CGAL::is_triangle_mesh(mesh)) {
    std::cerr << "Invalid input file." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Teleportation test." << std::endl;

  // algorithm instance
  Vertex_point_map vpmap = get(boost::vertex_point, const_cast<Mesh &>(mesh));
  L21_metric error_metric(mesh, vpmap);
  L21_approx approx(mesh, vpmap, error_metric);

  std::cout << "Random seeding by number." << std::endl;
  std::srand(static_cast<unsigned int>(std::time(nullptr)));

  std::size_t count = 0;
  while (!count) {
    approx.initialize_seeds(CGAL::parameters::seeding_method(CGAL::Surface_mesh_approximation::RANDOM)
      .max_number_of_proxies(50));
    if (approx.number_of_proxies() != 50)
      return EXIT_FAILURE;
    approx.run(10);

    std::cout << "Teleport until merge test failed." << std::endl;
    std::vector<FT> error;
    count = 0;
    while(approx.teleport_proxies(1) == 1) {
      FT sum_err(0);
      approx.run(10);
      sum_err += approx.compute_total_error();
      error.push_back(sum_err / FT(10.0));
      ++count;
    }
    std::cout << "#teleported " << count << std::endl;

    if (!check_strict_ordering(error)) {
      std::cout << "Failed: teleportation error decrease inconsistent." << std::endl;
      return EXIT_FAILURE;
    }
  }

  // test partition placement
  std::cout << "Test partition placement." << std::endl;
  Face_proxy_map fproxymap =
    mesh.add_property_map<face_descriptor, std::size_t>("f:porxy_id", 0).first;
  approx.proxy_map(fproxymap);
  std::vector<Plane_proxies> proxies;
  approx.proxies(std::back_inserter(proxies));


  CGAL::Bbox_3 bbox;
  for(const vertex_descriptor& v : vertices(mesh))
    bbox += vpmap[v].bbox();
  const FT ymin = bbox.ymin(), ymax = bbox.ymax(), yrange = ymax - ymin;
  std::cout << "Range along y axis: [" << ymin << ", " << ymax << "]" << std::endl;

  // test if all faces on the planar part are in the same proxy
  std::size_t planar_pxidx = static_cast<std::size_t>(-1);
  std::size_t num_planar_faces = 0;
  bool first = true;
  for(const face_descriptor& f : faces(mesh)) {
    const halfedge_descriptor he = halfedge(f, mesh);
    const Point_3 &p0 = vpmap[source(he, mesh)];
    const Point_3 &p1 = vpmap[target(he, mesh)];
    const Point_3 &p2 = vpmap[target(next(he, mesh), mesh)];
    const Point_3 fcenter = CGAL::centroid(p0, p1, p2);
    const Vector_3 fnormal = CGAL::unit_normal(p0, p1, p2);

    // check the face center and normal to see if it is on the planar part of the geometry
    double dis_var = std::abs(CGAL::to_double((fcenter.y() - ymin) / yrange));
    double dir_var = std::abs(CGAL::to_double(fnormal.y()) - 1.0);
    if (dis_var < CGAL_VSA_TEST_TOLERANCE && dir_var < CGAL_VSA_TEST_TOLERANCE) {
      ++num_planar_faces;
      const std::size_t pxidx = fproxymap[f];
      if (first) {
        first = false;
        planar_pxidx = pxidx;
      }

      // only one proxy on the planar part of the geometry
      if (pxidx != planar_pxidx) {
        std::cout << "Failed: more than one proxy on the planar part of the model." << std::endl;
        return EXIT_FAILURE;
      }
    }
  }
  std::cout << "#faces on planar part " << num_planar_faces << std::endl;
  if (num_planar_faces != 922) {
    std::cout << "Failed: the plane-sphere model have 922 faces on the planar part." << std::endl;
    return EXIT_FAILURE;
  }

  // proxy of the planar part should facing straight towards the y positive.
  double px_dir_var = std::abs(CGAL::to_double(proxies[planar_pxidx].y()) - 1.0);
  if (px_dir_var > CGAL_VSA_TEST_TOLERANCE) {
    std::cout << "Failed: the proxy of planar part is incorrect." << std::endl;
    return EXIT_FAILURE;
  }

  // force teleportation test
  if ( approx.find_best_merge(true) != std::nullopt )
  {
    std::cout << "Failed: should be no possible merge with test." << std::endl;
    return EXIT_FAILURE;
  }
  std::cout << "Test forced teleportation." << std::endl;
  if (approx.teleport_proxies(1, 5, true) != 1) {
    std::cout << "Failed: forced teleportation failed." << std::endl;
  }

  std::cout << "Succeeded." << std::endl;
  return EXIT_SUCCESS;
}
