#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
#include <CGAL/Polyhedral_envelope.h>
#include <CGAL/Polygon_mesh_processing/repair_self_intersections.h>

namespace PMP = CGAL::Polygon_mesh_processing;

typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;
typedef CGAL::Exact_predicates_exact_constructions_kernel EPEC;

void test_API()
{
  std::cout << "---- test_API() ----\n";
  std::vector<EPIC::Point_3> points;
  std::vector<std::array<int, 3> > triangles;
  std::ifstream in("data/eight.off");
  CGAL::IO::read_OFF(in, points, triangles);
  CGAL::Surface_mesh<EPIC::Point_3> sm;
  CGAL::Polyhedron_3<EPIC> poly;
  PMP::polygon_soup_to_polygon_mesh(points, triangles, sm);
  PMP::polygon_soup_to_polygon_mesh(points, triangles, poly);
  auto epsilon_map = sm.add_property_map<CGAL::Surface_mesh<EPIC::Point_3>::Face_index, double>("f:em", 0.0001).first;

  // build from the same kernel
  {
    CGAL::Polyhedral_envelope<EPIC> envelope;
    assert(envelope.is_empty());
    envelope = CGAL::Polyhedral_envelope<EPIC>(sm,0.0001);
    assert(!envelope.is_empty());
    assert(envelope(sm));
    assert(envelope(poly));
    assert(envelope(points, triangles));
    envelope = CGAL::Polyhedral_envelope<EPIC>(sm,0, CGAL::parameters::face_epsilon_map(epsilon_map));
    assert(!envelope.is_empty());
    assert(envelope(sm));
    assert(envelope(poly));
    assert(envelope(points, triangles));

    envelope=CGAL::Polyhedral_envelope<EPIC>(points,triangles,0.0001);
    assert(!envelope.is_empty());
    assert(envelope(sm));
    assert(envelope(poly));
    assert(envelope(points, triangles));
    std::vector<double> eps(triangles.size(), 0.0001);
    envelope=CGAL::Polyhedral_envelope<EPIC>(points,triangles,0,CGAL::parameters::face_epsilon_map(CGAL::make_property_map(eps)));
    assert(!envelope.is_empty());
    assert(envelope(sm));
    assert(envelope(poly));
    assert(envelope(points, triangles));

    std::vector<decltype(sm)::Face_index> subfaces(faces(sm).begin(), std::next(faces(sm).begin(), 150));
    envelope=CGAL::Polyhedral_envelope<EPIC>(subfaces,sm,0.0001);
    assert(!envelope.is_empty());
    assert(!envelope(sm));
    assert(!envelope(poly));
    assert(!envelope(points, triangles));
    envelope=CGAL::Polyhedral_envelope<EPIC>(subfaces,sm,0, CGAL::parameters::face_epsilon_map(epsilon_map));
    assert(!envelope.is_empty());
    assert(!envelope(sm));
    assert(!envelope(poly));
    assert(!envelope(points, triangles));
  }

  // build from different kernels
  {
    auto sm_np = CGAL::parameters::vertex_point_map(
      CGAL::make_cartesian_converter_property_map<EPEC::Point_3>(
        get(boost::vertex_point, sm)));
    auto poly_np = CGAL::parameters::vertex_point_map(
      CGAL::make_cartesian_converter_property_map<EPEC::Point_3>(
        get(boost::vertex_point, poly)));
    auto soup_np =  CGAL::parameters::point_map(
      CGAL::make_cartesian_converter_property_map<EPEC::Point_3>(
        CGAL::Identity_property_map<EPIC::Point_3>()));

    CGAL::Polyhedral_envelope<EPEC> envelope(poly,0.0001,poly_np);
    assert(envelope(sm, sm_np));
    assert(envelope(poly, poly_np));
    assert(envelope(points, triangles, soup_np));

    envelope = CGAL::Polyhedral_envelope<EPEC>(points,triangles,0.0001, soup_np);
    assert(envelope(sm, sm_np));
    assert(envelope(poly, poly_np));
    assert(envelope(points, triangles, soup_np));

    std::vector<decltype(sm)::Face_index> subfaces(faces(sm).begin(), std::next(faces(sm).begin(), 150));
    envelope=CGAL::Polyhedral_envelope<EPEC>(subfaces,sm,0.0001,sm_np);
    assert(!envelope.is_empty());
    assert(!envelope(sm, sm_np));
    assert(!envelope(poly, poly_np));
    assert(!envelope(points, triangles, soup_np));
  }
}

void test_remove_si()
{
  std::cout << "---- test_remove_si() ----\n";
  CGAL::Surface_mesh<EPIC::Point_3> tm;
  std::ifstream in("data/pig.off");
  in >> tm;
  assert(tm.vertices().size()!=0);

  // try with a small bound --> reject fix
  PMP::experimental::remove_self_intersections(tm, CGAL::parameters::polyhedral_envelope_epsilon(0.0001));
  assert(PMP::does_self_intersect(tm));
  // try with a larger bound --> accept fix
  PMP::experimental::remove_self_intersections(tm, CGAL::parameters::polyhedral_envelope_epsilon(0.001));
  assert(!PMP::does_self_intersect(tm));
}

void cube_test()
{
  std::cout << "---- cube_test() ----\n";
  CGAL::Surface_mesh<EPIC::Point_3> tm;
  std::ifstream in("data-coref/cube_meshed.off");
  in >> tm;
  assert(tm.vertices().size()!=0);

  typedef EPIC::Point_3 P;
  CGAL::Polyhedral_envelope<EPIC> envelope(tm, 0.1);
  assert(envelope(P(0,0,0), P(1,0,0), P(1,1,0)));
  assert(envelope(P(0,0,0.05), P(1,0,0.05), P(1,1,0.05)));
  // assert(envelope(P(0,0,0.1), P(1,0,0.1), P(1,1,0.1))); // NOT WORKING
  assert(!envelope(P(0,0,0.2), P(1,0,0.2), P(1,1,0.2)));
}

int main()
{
  test_remove_si();
  test_API();
  cube_test();
}
