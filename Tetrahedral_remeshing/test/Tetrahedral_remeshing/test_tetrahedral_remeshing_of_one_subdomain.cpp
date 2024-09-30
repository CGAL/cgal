//#define CGAL_TETRAHEDRAL_REMESHING_GENERATE_INPUT_FILES

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/tetrahedral_remeshing_io.h>

#include <CGAL/Random.h>

#include <iostream>
#include <fstream>
#include <cassert>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

void generate_input_two_subdomains(const std::size_t nbv, Remeshing_triangulation& tr)
{
  CGAL::Random rng;

  typedef Remeshing_triangulation::Point Point;
  typedef Remeshing_triangulation::Cell_handle Cell_handle;

  while (tr.number_of_vertices() < nbv)
    tr.insert(Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));

  const Remeshing_triangulation::Geom_traits::Plane_3
    plane(Point(0, 0, 0), Point(0, 1, 0), Point(0, 0, 1));

  for (Cell_handle c : tr.finite_cell_handles())
  {
    if (plane.has_on_positive_side(
      CGAL::centroid(c->vertex(0)->point(), c->vertex(1)->point(),
                     c->vertex(2)->point(), c->vertex(3)->point())))
      c->set_subdomain_index(1);
    else
      c->set_subdomain_index(2);
  }
  assert(tr.is_valid(true));

#ifdef CGAL_TETRAHEDRAL_REMESHING_GENERATE_INPUT_FILES
  std::ofstream os("data/triangulation_two_subdomains.binary.cgal",
                   std::ios_base::out | std::ios_base::binary);
  CGAL::save_binary_triangulation(os, tr);
  os.close();
#endif
}

template<typename Tr>
struct Cells_of_subdomain_pmap
{
private:
  using Cell_handle = typename Tr::Cell_handle;

  const int m_subdomain;

public:
  using key_type = Cell_handle;
  using value_type = bool;
  using reference = bool;
  using category = boost::read_write_property_map_tag;

  Cells_of_subdomain_pmap(const int& subdomain)
    : m_subdomain(subdomain)
  {}

  friend value_type get(
    const Cells_of_subdomain_pmap& map, const key_type& c)
  {
    return (map.m_subdomain == c->subdomain_index());
  }
  friend void put(
    Cells_of_subdomain_pmap&, const key_type&, const value_type)
  {} //nothing to do : subdomain indices are updated in remeshing

};

int main(int argc, char* argv[])
{
  std::cout << "CGAL Random seed = "
    << CGAL::get_default_random().get_seed() << std::endl;

  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;

  Remeshing_triangulation tr;
  generate_input_two_subdomains(1000, tr);

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length,
      CGAL::parameters::cell_is_selected_map(
        Cells_of_subdomain_pmap<Remeshing_triangulation>(2)));

  return EXIT_SUCCESS;
}
