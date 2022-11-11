#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/property_map.h>

#include <unordered_set>
#include <iostream>
#include <utility>
#include <cassert>

#include "tetrahedral_remeshing_generate_input.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K> Remeshing_triangulation;

typedef Remeshing_triangulation::Point         Point;
typedef Remeshing_triangulation::Vertex_handle Vertex_handle;
typedef Remeshing_triangulation::Cell_handle   Cell_handle;
typedef Remeshing_triangulation::Edge          Edge;

class Constrained_edges_property_map
{
public:
  typedef bool                               value_type;
  typedef bool                               reference;
  typedef std::pair<Vertex_handle, Vertex_handle> key_type;
  typedef boost::read_write_property_map_tag category;

private:
  std::unordered_set<key_type, boost::hash<key_type>>* m_set_ptr;

public:
  Constrained_edges_property_map()
    : m_set_ptr(nullptr)
  {}
  Constrained_edges_property_map(std::unordered_set<key_type, boost::hash<key_type>>* set_)
    : m_set_ptr(set_)
  {}

public:
  friend void put(Constrained_edges_property_map& map,
                  const key_type& k,
                  const bool b)
  {
    assert(map.m_set_ptr != nullptr);
    assert(k.first < k.second);
    if (b)  map.m_set_ptr->insert(k);
    else    map.m_set_ptr->erase(k);
  }

  friend value_type get(const Constrained_edges_property_map& map,
                        const key_type& k)
  {
    assert(map.m_set_ptr != nullptr);
    assert(k.first < k.second);
    return (map.m_set_ptr->count(k) > 0);
  }
};

void set_subdomain(Remeshing_triangulation& tr, const int index)
{
  for (auto cit : tr.finite_cell_handles())
    cit->set_subdomain_index(index);
}

int main(int argc, char* argv[])
{
  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.02;
  const int nb_iter = (argc > 2) ? atoi(argv[2]) : 1;
  const int nbv = (argc > 3) ? atoi(argv[3]) : 500;

  Remeshing_triangulation t3;
  typedef std::pair<Vertex_handle, Vertex_handle> Vertex_pair;
  std::unordered_set<Vertex_pair, boost::hash<Vertex_pair>> constraints;

  CGAL::Tetrahedral_remeshing::generate_input_cube(nbv, t3, constraints);
  make_constraints_from_cube_edges(t3, constraints);

  CGAL::tetrahedral_isotropic_remeshing(t3, target_edge_length,
    CGAL::parameters::edge_is_constrained_map(
      Constrained_edges_property_map(&constraints))
    .number_of_iterations(nb_iter));

  return EXIT_SUCCESS;
}
