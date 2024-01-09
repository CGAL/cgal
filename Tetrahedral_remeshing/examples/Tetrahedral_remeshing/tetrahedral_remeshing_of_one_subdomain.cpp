#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include "tetrahedral_remeshing_generate_input.h"

using K = CGAL::Exact_predicates_inexact_constructions_kernel;

using Remeshing_triangulation = CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K>;
using Subdomain_index = Remeshing_triangulation::Tds::Cell::Subdomain_index;

template<typename Tr>
struct Cells_of_subdomain_pmap
{
  const int m_subdomain;

public:
  using key_type = typename Tr::Cell_handle;
  using value_type = bool;
  using reference = bool;
  using category = boost::read_write_property_map_tag;

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

Subdomain_index subdomain_on_side_of_plane(const K::Point_3& p,
                                           const K::Plane_3& plane)
{
  if (plane.has_on_positive_side(p))
    return 1;
  else
    return 2;
}

K::Point_3 centroid(const Remeshing_triangulation::Cell_handle c)
{
  return CGAL::centroid(c->vertex(0)->point(),
                        c->vertex(1)->point(),
                        c->vertex(2)->point(),
                        c->vertex(3)->point());
}

int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.05;
  const std::size_t nbv = (argc > 2) ? atoi(argv[2]) : 1000;
  const std::size_t nbv_on_plane = (argc > 3) ? atoi(argv[3]) : 100;

  /// Create a randomly generated triangulation of a sphere
  Remeshing_triangulation tr;
  CGAL::Tetrahedral_remeshing::insert_random_points_in_cube(nbv, tr);

  /// Use vertical plane to split the mesh into two subdomains
  const K::Plane_3 plane(0, 0, 1, 0);
  CGAL::Tetrahedral_remeshing::insert_points_on_plane(plane, nbv_on_plane, tr);

  /// A subdomain index 0 is considered outside and is not remeshed
  /// so we set finite cells to a non-zero `Subdomain_index`
  /// (depending on the side of the plane they are on)
  for (auto cell : tr.finite_cell_handles())
  {
    const K::Point_3 cc = centroid(cell);
    cell->set_subdomain_index(subdomain_on_side_of_plane(cc, plane));
  }

  /// Remesh only the cells of subdomain 2
  CGAL::tetrahedral_isotropic_remeshing(tr,
      target_edge_length,
      CGAL::parameters::cell_is_selected_map(
        Cells_of_subdomain_pmap<Remeshing_triangulation>{2}));

  std::ofstream os("out.mesh");
  CGAL::IO::write_MEDIT(os, tr);

  return EXIT_SUCCESS;
}

