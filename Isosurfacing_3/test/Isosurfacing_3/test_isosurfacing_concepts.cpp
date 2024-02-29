#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Isosurfacing_3/marching_cubes_3.h>
#include <CGAL/Isosurfacing_3/dual_contouring_3.h>

#include <CGAL/assertions.h>
#include <CGAL/Origin.h>

#include <array>
#include <vector>

namespace IS = CGAL::Isosurfacing;

template <typename K>
struct Traits
{
  using FT = typename K::FT;
  using Point_3 = typename K::Point_3;
  using Vector_3 = typename K::Vector_3;
  using Iso_cuboid_3 = typename K::Iso_cuboid_3;

  struct Compute_x_3
  {
    FT operator()(const Point_3&) const { return 0; }
    FT operator()(const Vector_3&) const { return 0; }
  };

  struct Compute_y_3
  {
    FT operator()(const Point_3&) const { return 0; }
    FT operator()(const Vector_3&) const { return 0; }
  };

  struct Compute_z_3
  {
    FT operator()(const Point_3&) const { return 0; }
    FT operator()(const Vector_3&) const { return 0; }
  };

  struct Construct_point_3
  {
    Point_3 operator()(FT, FT, FT) const { return CGAL::ORIGIN; }
  };

  struct Construct_vector_3
  {
    Vector_3 operator()(FT, FT, FT) const { return CGAL::NULL_VECTOR; }
  };

  struct Construct_iso_cuboid_3
  {
    Iso_cuboid_3 operator()(FT, FT, FT, FT, FT, FT) const { return Iso_cuboid_3(); }
  };

  struct Construct_vertex_3
  {
    Point_3 operator()(const Iso_cuboid_3&, int) const { return CGAL::ORIGIN; }
  };

  Compute_x_3 compute_x_3_object() { return Compute_x_3(); }
  Compute_y_3 compute_y_3_object() { return Compute_y_3();}
  Compute_z_3 compute_z_3_object() { return Compute_z_3(); }
  Construct_point_3 construct_point_3_object() { return Construct_point_3();}
  Construct_vector_3 construct_vector_3_object() { return Construct_vector_3(); }
  Construct_iso_cuboid_3 construct_iso_cuboid_3_object() { return Construct_iso_cuboid_3(); }
  Construct_vertex_3 construct_vertex_3_object() { return Construct_vertex_3(); }
};

template <typename K>
struct MC_Domain
{
  typedef Traits<K> Geom_traits;
  typedef typename Geom_traits::FT FT;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef std::size_t vertex_descriptor;
  typedef std::size_t edge_descriptor;
  typedef std::size_t cell_descriptor;
  typedef std::vector<std::size_t> Vertices_incident_to_edge;
  typedef std::vector<std::size_t> Cells_incident_to_edge;
  typedef std::vector<std::size_t> Cell_vertices;
  typedef std::vector<std::size_t> Cell_edges;

  Geom_traits geom_traits() const { return Geom_traits(); }

  Point_3 point(const vertex_descriptor& /*v*/) const { return CGAL::ORIGIN; }

  FT value(const Point_3& /*p*/) const { return 0; }
  FT value(const vertex_descriptor& /*v*/) const { return 0; }

  Vertices_incident_to_edge incident_vertices(const edge_descriptor& /*e*/) const { return {}; }
  Cells_incident_to_edge incident_cells(const edge_descriptor& /*e*/) const { return {}; }
  Cell_vertices cell_vertices(const cell_descriptor& /*c*/) const { return {}; }
  Cell_edges cell_edges(const cell_descriptor& /*c*/) const { return {}; }

  template <typename ConcurrencyTag, typename Functor>
  void for_each_vertex(Functor& /*f*/) const { }
  template <typename ConcurrencyTag, typename Functor>
  void for_each_edge(Functor& /*f*/) const { }
  template <typename ConcurrencyTag, typename Functor>
  void for_each_cell(Functor& /*f*/) const { }

  bool construct_intersection(const Point_3&, const Point_3& /*p_1*/,
                              const FT /*v_0*/, const FT /*v_1*/,
                              const FT /*isovalue*/,
                              Point_3& /*intersection*/) const { return false; }
};

template <typename K>
struct DC_Domain
  : public MC_Domain<K>
{
  typedef Traits<K> Geom_traits;
  typedef typename Geom_traits::Point_3 Point_3;
  typedef typename Geom_traits::Vector_3 Vector_3;

  Vector_3 gradient(const Point_3& /*p*/) const { return CGAL::NULL_VECTOR; }
};

template <typename K>
void test()
{
  using Point_3 = typename K::Point_3;

  MC_Domain<K> mc_domain;
  DC_Domain<K> dc_domain;

  std::vector<Point_3> points;
  std::vector<std::vector<std::size_t> > faces;

  IS::marching_cubes(mc_domain, 0, points, faces);

  IS::dual_contouring(dc_domain, 0, points, faces);
}

int main(int, char**)
{
  test<CGAL::Simple_cartesian<double> >();
  test<CGAL::Exact_predicates_inexact_constructions_kernel>();

  std::cout << "Done" << std::endl;

  return EXIT_SUCCESS;
}
