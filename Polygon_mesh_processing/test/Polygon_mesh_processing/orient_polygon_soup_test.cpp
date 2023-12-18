#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/Cartesian_converter.h>
#include <CGAL/IO/polygon_mesh_io.h>
#include <CGAL/IO/polygon_soup_io.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orientation.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup_extension.h>

#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Epick;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epeck;

typedef Epeck::Point_3 Exact_point;

namespace PMP = CGAL::Polygon_mesh_processing;

template <typename Polygons>
void shuffle_soup(Polygons& polygons)
{
  // reverse orientation randomly
  for(std::size_t i=0; i<polygons.size(); ++i)
    if(std::rand() % 2 == 0)
      std::swap(polygons[i][0], polygons[i][1]);
}

template <typename K>
bool test_orient(const bool save_oriented)
{
  std::cout << "test_orient() with K = " << typeid(K).name() << std::endl;

  typedef typename K::Point_3 Point_3;
  typedef CGAL::Polyhedron_3<K> Polyhedron;
  typedef CGAL::Surface_mesh<Point_3> Surface_mesh;

  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > polygons;
  if(!CGAL::IO::read_polygon_soup(CGAL::data_file_path("meshes/elephant.off"), points, polygons))
  {
    std::cerr << "Error " << __LINE__ << ": failed to read polygon soup.\n";
    return false;
  }

  shuffle_soup(polygons);
  std::cout << "Is the soup a mesh? " << PMP::is_polygon_soup_a_polygon_mesh(polygons) << std::endl;

  bool oriented = PMP::orient_polygon_soup(points, polygons);
  std::cerr << (oriented ? "Oriented." : "Not orientabled.") << std::endl;
  if(!oriented || !PMP::is_polygon_soup_a_polygon_mesh(polygons))
    return false;

  if(oriented)
  {
    Surface_mesh mesh;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, mesh);
    if(!is_valid_polygon_mesh(mesh))
      return false;

    Polyhedron poly;
    PMP::polygon_soup_to_polygon_mesh(points, polygons, poly);

    if(save_oriented)
      CGAL::IO::write_polygon_mesh("elephant-oriented.stl", poly, CGAL::parameters::use_binary_mode(false));

    if(!is_valid_polygon_mesh(poly))
      return false;
  }

  return true;
}

template <class K, class Tag>
bool test_pipeline()
{
  std::cout << "test_pipeline() with K = " << typeid(K).name() << std::endl;

  typedef typename K::Point_3 Point_3;
  typedef CGAL::Cartesian_converter<K, Epeck> K2E;

  typedef CGAL::Polyhedron_3<K> Polyhedron;

  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > polygons;
  if(!CGAL::IO::read_polygon_soup(CGAL::data_file_path("meshes/elephant.off"), points, polygons))
  {
    std::cerr << "Error " << __LINE__ << ": failed to read polygon soup.\n";
    return false;
  }

  std::map<Point_3, Exact_point> k2e;
  boost::associative_property_map<std::map<Point_3, Exact_point>> cpm(k2e);
  for(const auto& p : points)
    put(cpm, p, K2E()(p));

  shuffle_soup(polygons);
  std::cout << "Is the soup a mesh? " << PMP::is_polygon_soup_a_polygon_mesh(polygons) << std::endl;

  Polyhedron ref1;
  if(!CGAL::IO::read_polygon_mesh(CGAL::data_file_path("meshes/elephant.off"), ref1))
  {
    std::cerr << "Error " << __LINE__ << ": failed to read reference mesh.\n";
    return false;
  }

  typedef CGAL::dynamic_vertex_property_t<Exact_point>                      Point_property;
  typedef typename boost::property_map<Polyhedron, Point_property>::type    Custom_VPM;
  Custom_VPM cvpm = get(Point_property(), ref1);
  for(auto v : vertices(ref1))
    put(cvpm, v, K2E()(get(CGAL::vertex_point, ref1, v)));

  if(PMP::is_outward_oriented(ref1, CGAL::parameters::vertex_point_map(cvpm)))
    PMP::reverse_face_orientations(ref1);

  PMP::orient_triangle_soup_with_reference_triangle_mesh<Tag>(ref1, points, polygons,
                                                              CGAL::parameters::vertex_point_map(cvpm),
                                                              CGAL::parameters::point_map(cpm));

  PMP::duplicate_non_manifold_edges_in_polygon_soup(points, polygons);

  Polyhedron poly;
  PMP::polygon_soup_to_polygon_mesh(points, polygons, poly);

  typedef typename boost::property_map<Polyhedron, CGAL::dynamic_face_property_t<std::size_t> >::type Fccmap;
  Fccmap fim = get(CGAL::dynamic_face_property_t<std::size_t>(), poly);
  std::size_t id =0;
  for(auto f : faces(poly))
    put(fim, f, id++);

  PMP::merge_reversible_connected_components(poly, CGAL::parameters::face_index_map(fim));

  Fccmap fccmap = get(CGAL::dynamic_face_property_t<std::size_t>(), poly);
  if(PMP::connected_components(poly, fccmap, CGAL::parameters::face_index_map(fim)) != 1)
    return false;

  if(PMP::is_outward_oriented(poly))
    return false;

  PMP::reverse_face_orientations(ref1);

  std::vector<Point_3> ref_points;
  std::vector<std::vector<std::size_t> > ref_polygons;
  PMP::polygon_mesh_to_polygon_soup(ref1, ref_points, ref_polygons);

  std::map<Point_3, Exact_point> ref_k2e;
  boost::associative_property_map<std::map<Point_3, Exact_point> > ref_cpm(ref_k2e);
  for(const auto& p : ref_points)
    put(ref_cpm, p, K2E()(p));

  shuffle_soup(polygons);
  std::cout << "Is the soup a mesh? " << PMP::is_polygon_soup_a_polygon_mesh(polygons) << std::endl;

  PMP::orient_triangle_soup_with_reference_triangle_soup(ref_points, ref_polygons, points, polygons,
                                                         CGAL::parameters::point_map(ref_cpm),
                                                         CGAL::parameters::point_map(cpm));
  if(!PMP::is_polygon_soup_a_polygon_mesh(polygons))
  {
    std::cerr << "Error: Orient_TS_with_ref_TS failed" << std::endl;
    return false;
  }

  CGAL::clear(poly);
  PMP::polygon_soup_to_polygon_mesh(points, polygons, poly);
  if(!is_valid_polygon_mesh(poly) || !PMP::is_outward_oriented(poly))
  {
    std::cerr << "Error: result is not invalid or not outward oriented" << std::endl;
    return false;
  }

  return true;
}

int main()
{
  bool res = test_orient<Epick>(false /*save_oriented*/);
  assert(res);
  res = test_orient<Epeck>(false /*save_oriented*/);
  assert(res);

  res = test_pipeline<Epick, CGAL::Sequential_tag>();
  assert(res);
  res = test_pipeline<Epeck, CGAL::Sequential_tag>();
  assert(res);

#if defined(CGAL_LINKED_WITH_TBB)
  res = test_pipeline<Epick, CGAL::Parallel_tag>();
  assert(res);

  res = test_pipeline<Epeck, CGAL::Parallel_tag>();
  assert(res);
#endif

  return 0;
}
