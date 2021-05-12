#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#define CGAL_TETRAHEDRAL_REMESHING_GENERATE_INPUT_FILES

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>
#include <CGAL/Tetrahedral_remeshing/tetrahedral_remeshing_io.h>

#include <CGAL/Random.h>
#include <CGAL/property_map.h>

#include <boost/unordered_set.hpp>

#include <iostream>
#include <fstream>
#include <vector>
#include <utility>

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
    : m_set_ptr(nullptr)
  {}
  Constrained_edges_property_map(boost::unordered_set<key_type>* set_)
    : m_set_ptr(set_)
  {}

public:
  friend void put(Constrained_edges_property_map& map,
                  const key_type& k,
                  const bool b)
  {
    CGAL_assertion(map.m_set_ptr != nullptr);
    CGAL_assertion(k.first < k.second);
    if (b)  map.m_set_ptr->insert(k);
    else    map.m_set_ptr->erase(k);
  }

  friend value_type get(const Constrained_edges_property_map& map,
                              const key_type& k)
  {
    CGAL_assertion(map.m_set_ptr != nullptr);
    CGAL_assertion(k.first < k.second);
    return map.m_set_ptr->count(k) > 0;
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

void generate_input_cube(const std::size_t& n,
              Remeshing_triangulation& tr,
              boost::unordered_set<std::pair<Vertex_handle, Vertex_handle> >& constraints)
{
  CGAL::Random rng;

  // points in a sphere
  std::vector<Point> pts;
  while (pts.size() < n)
    pts.push_back(Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));
  tr.insert(pts.begin(), pts.end());

  // vertices of a larger cube
  Vertex_handle v0 = tr.insert(Point(-2., -2., -2.));
  Vertex_handle v1 = tr.insert(Point(-2., -2.,  2.));

  Vertex_handle v2 = tr.insert(Point( 2., -2., -2.));
  Vertex_handle v3 = tr.insert(Point( 2., -2.,  2.));

  Vertex_handle v4 = tr.insert(Point(-2.,  2.,  -2.));
  Vertex_handle v5 = tr.insert(Point(-2.,  2.,   2.));

  Vertex_handle v6 = tr.insert(Point( 2.,  2., -2.));
  Vertex_handle v7 = tr.insert(Point( 2.,  2.,  2.));

  CGAL_assertion(tr.is_valid(true));

  // writing file output
#ifdef CGAL_TETRAHEDRAL_REMESHING_GENERATE_INPUT_FILES
  std::ofstream outfile("data/sphere_in_cube.tr.cgal",
                        std::ios_base::out | std::ios_base::binary);
  CGAL::save_binary_triangulation(outfile, tr);
  outfile.close();
#endif

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

  CGAL_assertion(tr.is_valid(true));
}

void set_subdomain(Remeshing_triangulation& tr, const int index)
{
  for (Cell_handle c : tr.finite_cell_handles())
    c->set_subdomain_index(index);
}

int main(int argc, char* argv[])
{
  CGAL::Random rng;
  std::cout << "CGAL Random seed = " << CGAL::get_default_random().get_seed() << std::endl;

  Remeshing_triangulation tr;
  boost::unordered_set<std::pair<Vertex_handle, Vertex_handle> > constraints;
  generate_input_cube(1000, tr, constraints);

  const double target_edge_length = (argc > 1) ? atof(argv[1]) : 0.02;
  const int nb_iter = (argc > 2) ? atoi(argv[2]) : 1;

  set_subdomain(tr, 1);
  assert(tr.is_valid());

  CGAL::tetrahedral_isotropic_remeshing(tr, target_edge_length,
    CGAL::parameters::edge_is_constrained_map(
      Constrained_edges_property_map(&constraints))
    .number_of_iterations(nb_iter));

  return EXIT_SUCCESS;
}

