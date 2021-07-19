#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>

#include <CGAL/boost/graph/helpers.h>
#include <CGAL/boost/graph/property_maps.h>
#include <CGAL/Dynamic_property_map.h>
#include <CGAL/Polygon_mesh_processing/polygon_soup_to_polygon_mesh.h>
#include <CGAL/Polygon_mesh_processing/polygon_mesh_to_polygon_soup.h>
#include <CGAL/Polygon_mesh_processing/orient_polygon_soup.h>
#include <CGAL/IO/OFF.h>

#include <array>
#include <string>
#include <fstream>
#include <iostream>

typedef CGAL::Simple_cartesian<double> SC;
typedef CGAL::Exact_predicates_exact_constructions_kernel Epec;


template <typename K>
void test_polygon_soup(std::string fname, bool expected)
{
  typedef typename K::Point_3                                                 Point;

  typedef CGAL::Polyhedron_3<K>                                               Polyhedron;
  typedef typename boost::graph_traits<Polyhedron>::vertex_descriptor         vertex_descriptor;

  std::vector<Point> points;
  std::vector< std::vector<std::size_t> > polygons;
  std::ifstream input(fname.c_str());

  if(!input)
  {
    std::cerr << "Error opening file " << fname << "\n";
    exit(EXIT_FAILURE);
  }

  if(!CGAL::IO::read_OFF(input, points, polygons))
  {
    std::cerr << "Error parsing the OFF file " << fname << "\n";
    exit(EXIT_FAILURE);
  }

  bool is_mesh = CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons);
  std::cout << "is_polygon_soup_a_polygon_mesh(" << fname << ") == "
            << std::boolalpha << is_mesh << ";" << std::endl;
  assert(is_mesh == expected);

  if(is_mesh)
  {
    Polyhedron p;

    // just to test the named paramers
    typedef std::pair<Point, bool>                                            Point_with_Boolean;
    std::vector<Point_with_Boolean> points_with_pairs;
    for(const Point& pt : points)
      points_with_pairs.emplace_back(pt, false);

    typedef CGAL::dynamic_vertex_property_t<Point>                            Point_property;
    typedef typename boost::property_map<Polyhedron, Point_property>::type    Custom_VPM;

    Custom_VPM vpm = get(Point_property(), p);

    CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(
          points_with_pairs, polygons, p,
          CGAL::parameters::point_map(CGAL::First_of_pair_property_map<Point_with_Boolean>()),
          CGAL::parameters::vertex_point_map(vpm));

    std::cout << num_vertices(p) << " nv and " << num_faces(p) << " nf" << std::endl;
    assert(!CGAL::is_empty(p) && CGAL::is_valid_polygon_mesh(p));

    std::set<Point> ppts;
    for( vertex_descriptor v : vertices(p))
      ppts.insert(get(vpm, v));

    assert(ppts.size() == num_vertices(p));

    // twice to check if adds correctly
    std::deque<Point> soup_points;
    std::vector<std::array<std::size_t, 3> > soup_polygons;

    CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(p, soup_points, soup_polygons);
    CGAL::Polygon_mesh_processing::polygon_mesh_to_polygon_soup(p, soup_points, soup_polygons);

    std::size_t nv = static_cast<std::size_t>(num_vertices(p));
    std::size_t nf = static_cast<std::size_t>(num_faces(p));

    assert(soup_points.size() == 2 * nv);
    assert(soup_polygons.size() == 2 * nf);

    // check sanity of the polygons
    for(std::size_t fi=0; fi<nf; ++fi)
    {
      for(const std::size_t pi : soup_polygons[fi]) {
        assert(pi < nv);
      }
    }

    for(std::size_t fi=nf; fi<2*nf; ++fi)
    {
      for(const std::size_t pi : soup_polygons[fi]) {
        assert(nv <= pi && pi < 2 * nv);
      }
    }
  }

  if(!expected)
  {
    CGAL::Polygon_mesh_processing::orient_polygon_soup(points, polygons);
    bool is_mesh = CGAL::Polygon_mesh_processing::is_polygon_soup_a_polygon_mesh(polygons);
    std::cout << "After orientation: is_polygon_soup_a_polygon_mesh(" << fname << ") == "
              << std::boolalpha << is_mesh << ";" << std::endl;
    if(is_mesh)
    {
      Polyhedron p;
      CGAL::Polygon_mesh_processing::polygon_soup_to_polygon_mesh(points, polygons, p);
      std::cout << num_vertices(p) << " nv and " << num_faces(p) << " nf" << std::endl;
      assert(!CGAL::is_empty(p) && CGAL::is_valid_polygon_mesh(p));
    }
  }

  std::cout << fname << " OK\n\n\n";
}

int main()
{
  test_polygon_soup<SC>("data_polygon_soup/bad_cube.off", false);
  test_polygon_soup<Epec>("data_polygon_soup/bad_cube.off", false);
  test_polygon_soup<SC>("data_polygon_soup/isolated_singular_vertex_one_cc.off", false);

  test_polygon_soup<SC>("data_polygon_soup/isolated_vertices.off", false);

  test_polygon_soup<SC>("data_polygon_soup/nm_vertex_and_edge.off", false);
  test_polygon_soup<SC>("data_polygon_soup/one_duplicated_edge.off", false);
  test_polygon_soup<SC>("data_polygon_soup/one_duplicated_edge_sharing_vertex.off", false);
  test_polygon_soup<SC>("data_polygon_soup/partial_overlap.off", false);
  test_polygon_soup<SC>("data_polygon_soup/incompatible_orientation.off", false);

  test_polygon_soup<SC>("data/blobby_3cc.off", true);
  test_polygon_soup<SC>("data/elephant.off", true);
  test_polygon_soup<SC>("data/joint_refined.off", true);
  test_polygon_soup<SC>("data/mech-holes-shark.off", true);
  test_polygon_soup<SC>("data/non_manifold_vertex.off", false);
  test_polygon_soup<SC>("data/two_tris_collinear.off", true);
  test_polygon_soup<SC>("data/U.off", true);

  return EXIT_SUCCESS;
}
