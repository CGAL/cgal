#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#define CGAL_DUMP_REMESHING_STEPS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

#include <CGAL/Random.h>
#include <CGAL/property_map.h>

#include <boost/unordered_set.hpp>

#include <iostream>
#include <fstream>
#include <utility>

#include "tetrahedral_remeshing_io.h"

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
  boost::unordered_set<key_type>* m_set_ptr;

public:
  Constrained_edges_property_map()
    : m_set_ptr(NULL)
  {}
  Constrained_edges_property_map(boost::unordered_set<key_type>* set_)
    : m_set_ptr(set_)
  {}

public:
  friend void put(Constrained_edges_property_map& map,
                  const key_type& k,
                  const bool b)
  {
    CGAL_assertion(map.m_set_ptr != NULL);
    CGAL_assertion(k.first < k.second);
    if (b)  map.m_set_ptr->insert(k);
    else    map.m_set_ptr->erase(k);
  }

  friend value_type get(const Constrained_edges_property_map& map,
                              const key_type& k)
  {
    CGAL_assertion(map.m_set_ptr != NULL);
    CGAL_assertion(k.first < k.second);
    return map.m_set_ptr->count(k);
  }
};

void add_edge(Vertex_handle v1,
              Vertex_handle v2,
              const Remeshing_triangulation& tr,
              boost::unordered_set<std::pair<Vertex_handle, Vertex_handle> >& constraints)
{
  Cell_handle c;
  int i, j;
  if(tr.is_edge(v1, v2, c, i, j))
    constraints.insert(std::make_pair(v1, v2));
}

void make_constraints_from_cube_edges(
              Remeshing_triangulation& tr,
              boost::unordered_set<std::pair<Vertex_handle, Vertex_handle> >& constraints)
{
  Remeshing_triangulation::Locate_type lt;
  int li, lj;

  Cell_handle c = tr.locate(Point(-2., -2., -2.), lt, li, lj);
  Vertex_handle v0 = c->vertex(li);
  c = tr.locate(Point(-2., -2.,  2.));
  Vertex_handle v1 = c->vertex(li);

  c = tr.locate(Point( 2., -2., -2.));
  Vertex_handle v2 = c->vertex(li);
  c = tr.locate(Point( 2., -2.,  2.));
  Vertex_handle v3 = c->vertex(li);

  c = tr.locate(Point(-2.,  2.,  -2.));
  Vertex_handle v4 = c->vertex(li);
  c = tr.locate(Point(-2.,  2.,   2.));
  Vertex_handle v5 = c->vertex(li);

  c = tr.locate(Point( 2.,  2., -2.));
  Vertex_handle v6 = c->vertex(li);
  c = tr.locate(Point( 2.,  2.,  2.));
  Vertex_handle v7 = c->vertex(li);

  // constrain cube edges
  add_edge(v0, v1, tr, constraints);
  add_edge(v1, v2, tr, constraints);
  add_edge(v2, v3, tr, constraints);
  add_edge(v3, v0, tr, constraints);

  add_edge(v4, v5, tr, constraints);
  add_edge(v5, v6, tr, constraints);
  add_edge(v6, v7, tr, constraints);
  add_edge(v7, v4, tr, constraints);

  add_edge(v0, v4, tr, constraints);
  add_edge(v1, v5, tr, constraints);
  add_edge(v2, v6, tr, constraints);
  add_edge(v3, v7, tr, constraints);
}

void set_subdomain(Remeshing_triangulation& tr, const int index)
{
  for (Remeshing_triangulation::Finite_cells_iterator cit = tr.finite_cells_begin();
       cit != tr.finite_cells_end(); ++cit)
  {
    cit->set_subdomain_index(index);
  }
}

int main(int argc, char* argv[])
{
  const char* filename     = (argc > 1) ? argv[1] : "data/sphere_in_cube.tr.cgal";
  double target_edge_length = (argc > 2) ? atof(argv[2]) : 0.02;
  int nb_iter = (argc > 3) ? atoi(argv[3]) : 1;

  std::ifstream input(filename, std::ios::in);
  if (!input)
  {
    std::cerr << "File " << filename << " could not be found" << std::endl;
    return EXIT_FAILURE;
  }

  Remeshing_triangulation t3;
  input >> t3;
  set_subdomain(t3, 1);
  boost::unordered_set<std::pair<Vertex_handle, Vertex_handle> > constraints;
  make_constraints_from_cube_edges(t3, constraints);

  CGAL_assertion(t3.is_valid());

  CGAL::tetrahedral_adaptive_remeshing(t3, target_edge_length,
    CGAL::parameters::edge_is_constrained_map(
      Constrained_edges_property_map(&constraints))
    .number_of_iterations(nb_iter));

  save_ascii_triangulation("tet_remeshing_with_features_after.mesh", t3);

  return EXIT_SUCCESS;
}

