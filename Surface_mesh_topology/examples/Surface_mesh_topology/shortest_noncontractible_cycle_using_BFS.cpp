#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <CGAL/Shortest_noncontractible_cycle.h>
#include <CGAL/IO/Color.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/squared_distance_3.h>
#define CGAL_USE_BASIC_VIEWER 1

using LCC_3 = CGAL::Linear_cell_complex_for_combinatorial_map<2, 3>;
using size_type = LCC_3::size_type;
using SNC = CGAL::Surface_mesh_topology::Shortest_noncontractible_cycle<LCC_3>;

LCC_3 lcc;

struct Draw_functor : public CGAL::DefaultDrawingFunctorLCC {
  Draw_functor(size_type edge_mark, size_type vertex_mark) : m_edge_mark(edge_mark), m_vertex_mark(vertex_mark) {}
  
  template<typename LCC>
  bool colored_vertex(const LCC& alcc, typename LCC::Dart_const_handle dh) const 
  {return alcc.is_marked(dh, m_vertex_mark);}
  
  template<typename LCC>
  CGAL::Color vertex_color(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  {return CGAL::Color(255, 0, 0);}

  template<typename LCC>
  bool colored_edge(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return alcc.is_marked(dh, m_edge_mark); }

  template<typename LCC>
  CGAL::Color edge_color(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  {return CGAL::Color(0, 0, 255);}

  template<typename LCC>
  bool colored_face(const LCC& alcc, typename LCC::Dart_const_handle dh) const {return true;}

  template<typename LCC>
  CGAL::Color face_color(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  {return CGAL::Color(211, 211, 211);}

  template<typename LCC>
  bool colored_volume(const LCC& alcc, typename LCC::Dart_const_handle dh) const { return false; }

private:
  size_type m_edge_mark, m_vertex_mark;
};

int main(int argc, char* argv[]) {
  std::cout << "Program shortest_noncontractible_cycle_using_BFS started.\n";
  std::ifstream inp;
  if (argc == 1) inp = std::ifstream("data/3torus-smooth.off");
  else inp = std::ifstream(argv[1]);
  CGAL::load_off(lcc, inp);
  std::cout << "File loaded. Running the main program...\n";
  
  SNC snc(lcc);
  SNC::Path cycle;
  SNC::Distance_type x;
  auto root = lcc.darts().begin();
  std::cout << "Finding the shortest noncontractible cycle...\n";
  snc.find_cycle(root, cycle, &x);
  if (cycle.size() == 0) {
    std::cout << "  Cannot find such cycle. Stop.\n";
    return 0;
  }
  size_type m = lcc.get_new_mark(), is_root = lcc.get_new_mark();
  for (auto e : cycle) {
    lcc.mark_cell<1>(e, m);
  }
  lcc.mark_cell<0>(root, is_root);
  std::cout << "  Number of edges in cycle: " << cycle.size() << std::endl;
  std::cout << "  Cycle length: " << x << std::endl;
  std::cout << "  Root: " << lcc.point_of_vertex_attribute(lcc.vertex_attribute(root)) << std::endl;
  Draw_functor df(m, is_root);
  CGAL::draw(lcc, "Hello", false, df);
  lcc.free_mark(m);
  lcc.free_mark(is_root);
}
