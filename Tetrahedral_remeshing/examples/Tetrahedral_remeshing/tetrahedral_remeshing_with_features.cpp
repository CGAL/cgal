#include <iostream>
#include <fstream>
#include <utility>

#define CGAL_TETRAHEDRAL_REMESHING_VERBOSE
#define CGAL_DUMP_REMESHING_STEPS

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Tetrahedral_remeshing/Remeshing_triangulation_3.h>
#include <CGAL/tetrahedral_remeshing.h>

//#include <CGAL/Tetrahedral_remeshing/internal/tetrahedral_remeshing_helpers.h>

#include <CGAL/Random.h>
#include <CGAL/property_map.h>

#include <boost/unordered_set.hpp>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Tetrahedral_remeshing::Remeshing_triangulation_3<K, int> Remeshing_triangulation;

typedef Remeshing_triangulation::Point         Point;
typedef Remeshing_triangulation::Vertex_handle Vertex_handle;
typedef Remeshing_triangulation::Cell_handle   Cell_handle;
typedef Remeshing_triangulation::Edge          Edge;

template <typename T3>
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

  friend const value_type get(const Constrained_edges_property_map& map,
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
    constraints.insert(std::make_pair(c->vertex(i), c->vertex(j)));
}

void generate_input(const std::size_t& n,
                    const char* filename,
                    boost::unordered_set<std::pair<Vertex_handle, Vertex_handle> >& constraints)
{
  Remeshing_triangulation tr;
  CGAL::Random rng;

  // points in a sphere
  while (tr.number_of_vertices() < n)
    tr.insert(Point(rng.get_double(-1., 1.), rng.get_double(-1., 1.), rng.get_double(-1., 1.)));
  // vertices of a larger cube
  Vertex_handle v0 = tr.insert(Point(-2., -2., -2.));
  Vertex_handle v1 = tr.insert(Point(-2., -2.,  2.));

  Vertex_handle v2 = tr.insert(Point( 2., -2., -2.));
  Vertex_handle v3 = tr.insert(Point( 2., -2.,  2.));
  
  Vertex_handle v4 = tr.insert(Point(-2.,  2.,  -2.));
  Vertex_handle v5 = tr.insert(Point(-2.,  2.,   2.));

  Vertex_handle v6 = tr.insert(Point( 2.,  2., -2.));
  Vertex_handle v7 = tr.insert(Point( 2.,  2.,  2.));

  // writing file output
  std::ofstream oFileT(filename, std::ios::out);
  oFileT << tr;
  oFileT.close();

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

int main(int argc, char* argv[])
{
  boost::unordered_set<std::pair<Vertex_handle, Vertex_handle> > constraints;
  generate_input(1000, "data/sphere_in_cube.tr.cgal", constraints);

  const char* filename     = (argc > 1) ? argv[1] : "data/sphere_in_cube.tr.cgal";
  float target_edge_length = (argc > 2) ? atof(argv[2]) : 0.1f;

  std::ifstream input(filename, std::ios::in);
  if (!input)
  {
    std::cerr << "File " << filename << " could not be found" << std::endl;
    return EXIT_FAILURE;
  }

  Remeshing_triangulation t3;
  input >> t3;
  CGAL_assertion(t3.is_valid());

  CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(t3,
    "tet_remeshing_with_features_before.mesh");

  CGAL::tetrahedral_adaptive_remeshing(t3, target_edge_length,
    CGAL::parameters::edge_is_constrained_map(
      Constrained_edges_property_map<Remeshing_triangulation>(&constraints)));

  std::ofstream oFileT("output.tr.cgal", std::ios::out);
  // writing file output;
  oFileT << t3;

  CGAL::Tetrahedral_remeshing::debug::dump_triangulation_cells(t3,
    "tet_remeshing_with_features_after.mesh");

  return EXIT_SUCCESS;
}

