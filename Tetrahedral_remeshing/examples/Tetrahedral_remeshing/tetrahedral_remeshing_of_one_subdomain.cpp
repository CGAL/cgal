#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include "tetrahedral_remeshing_generate_input.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;


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

  friend value_type get(const Cells_of_subdomain_pmap& map,
                        const key_type& c)
  {
    return (map.m_subdomain == c->subdomain_index());
  }
  friend void put(Cells_of_subdomain_pmap&,
                  const key_type&,
                  const value_type)
  {
    ; //nothing to do : subdomain indices are updated in remeshing
  }
};

int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.1;
  const std::size_t nbv = (argc > 2) ? atoi(argv[2]) : 1000;

  Remeshing_triangulation tr;
  CGAL::Tetrahedral_remeshing::generate_input_two_subdomains(nbv, tr);

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length,
      CGAL::parameters::cell_is_selected_map(
        Cells_of_subdomain_pmap<Remeshing_triangulation>(2)));

  return EXIT_SUCCESS;
}

