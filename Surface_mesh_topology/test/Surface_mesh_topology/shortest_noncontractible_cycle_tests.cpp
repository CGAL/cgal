#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_3.h>
#include <iostream>
#include <cstdlib>
#include <CGAL/Shortest_noncontractible_cycle.h>

/*
  Test the following cases:
  Input type: GMap, LCC_for_CMap, Surface_mesh
  Orientation: oriented, non-oriented
  Distance: weighted, unweighted
  Function: find_cycle, edge_width
*/

using GMap_2 = CGAL::Generalized_map<2>;
using LCC_for_CMap_2 = CGAL::Linear_cell_complex_for_combinatorial_map<2,3>;
using LCC_for_GMap_2 = CGAL::Linear_cell_complex_for_generalized_map<2,3>;
using Kernel = CGAL::Simple_cartesian<double>;
using Point = Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;

struct Weight_functor_for_SM {
  using Weight_t = double;
  Weight_functor_for_SM(const Surface_mesh& mesh) : m_mesh(mesh) {}
  double operator()(Surface_mesh::Halfedge_index he) const {
    Point A = m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 0));
    Point B = m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 1));
    return CGAL::sqrt(CGAL::squared_distance(A, B));
  }
private:
  Surface_mesh m_mesh;
};

template <class LCC_3>
struct Weight_functor_for_LCC {
  Weight_functor_for_LCC(const LCC_3& lcc) : m_lcc(lcc) { }
  using Weight_t = double;
  Weight_t operator()(typename LCC_3::Dart_const_handle dh) const {
    auto x = m_lcc.point_of_vertex_attribute(m_lcc.vertex_attribute(dh));
    auto y = m_lcc.point_of_vertex_attribute(m_lcc.vertex_attribute(m_lcc.next(dh)));
    return CGAL::sqrt(CGAL::squared_distance(x, y));
  }
private:
  const LCC_3& m_lcc;
};

bool get_data(Surface_mesh& sm) {
  std::ifstream in("./data/3torus.off");
  if (in.fail()) return false;
  in >> sm;
  return true;
}

template <class T>
bool get_data(T& lcc) {
  std::ifstream in("./data/3torus.off");
  if (in.fail()) return false;
  CGAL::load_off(lcc, in);
  return true;
}

bool find_cycle_in_unweighted_cmap() {
  LCC_for_CMap_2 lcc;
  if (!get_data(lcc)) {
    std::cout << "Fail find_cycle_in_unweighted_cmap: Cannot locate file data/3torus.off\n";
    return false;
  }
  CGAL::Shortest_noncontractible_cycle<LCC_for_CMap_2> snc(lcc);
  CGAL::Shortest_noncontractible_cycle<LCC_for_CMap_2>::Path cycle;
  CGAL::Shortest_noncontractible_cycle<LCC_for_CMap_2>::Distance_type cycle_length;
  LCC_for_CMap_2::Dart_handle root = lcc.darts().begin();
  snc.find_cycle(root, cycle, &cycle_length);
  if (cycle.size() != cycle_length) {
    std::cout << "Fail find_cycle_in_unweighted_cmap: cycle.size() != cycle_length\n";
    return false;
  }
  for (auto e : cycle)
    if (e == NULL) {
      std::cout << "Fail find_cycle_in_unweighted_cmap: NULL dart handle found in cycle\n";
      return false;
    }
  return true;
}

bool find_cycle_in_unweighted_gmap() { // Make a non-oriented case here
  return true;
}

bool edge_width_in_weighted_cmap_gmap_mesh() {
  LCC_for_CMap_2 lcc_cm;
  LCC_for_GMap_2 lcc_gm;
  Surface_mesh sm;
  if (!get_data(lcc_cm) || !get_data(lcc_gm) || !get_data(sm)) {
    std::cout << "Fail edge_width_in_weighted_cmap_gmap_mesh: Cannot locate file data/3torus.off\n";
    return false;
  }
  Weight_functor_for_LCC<LCC_for_CMap_2> wf_cm(lcc_cm);
  Weight_functor_for_LCC<LCC_for_GMap_2> wf_gm(lcc_gm);
  Weight_functor_for_SM wf_sm(sm);
  typedef CGAL::Shortest_noncontractible_cycle<LCC_for_CMap_2, Weight_functor_for_LCC<LCC_for_CMap_2> > SNC_1;
  typedef CGAL::Shortest_noncontractible_cycle<LCC_for_GMap_2, Weight_functor_for_LCC<LCC_for_GMap_2> > SNC_2;
  typedef CGAL::Shortest_noncontractible_cycle<Surface_mesh, Weight_functor_for_SM> SNC_3;
  SNC_1 snc1(lcc_cm, wf_cm);
  SNC_2 snc2(lcc_gm, wf_gm);
  SNC_3 snc3(sm, wf_sm);
  SNC_1::Path cycle1;
  SNC_2::Path cycle2;
  SNC_3::Path cycle3;
  double cycle_length1, cycle_length2, cycle_length3;
  snc1.edge_width(cycle1, &cycle_length1);
  snc2.edge_width(cycle2, &cycle_length2);
  snc3.edge_width(cycle3, &cycle_length3);
  for (auto e : cycle1)
    if (e == NULL) {
      std::cout << "Fail edge_width_in_weighted_cmap_gmap_mesh: NULL dart handle found in cycle\n";
      return false;
    }
  for (auto e : cycle2)
    if (e == NULL) {
      std::cout << "Fail edge_width_in_weighted_cmap_gmap_mesh: NULL dart handle found in cycle\n";
      return false;
    }
  if (cycle1.size() != cycle2.size() || cycle1.size() != cycle3.size()) {
    std::cout << "Fail edge_width_in_weighted_cmap_gmap_mesh: Inconsistency in number of edges of the edge-width "
              << "(" << cycle1.size() << ", " << cycle2.size() << ", " << cycle3.size() << ").\n";
    return false;
  }
  if (cycle_length1 - cycle_length2 > 1e-5 || cycle_length1 - cycle_length3 > 1e-5) {
    std::cout << "Fail edge_width_in_weighted_cmap_gmap_mesh: Inconsistency in the edge-width length"
              << std::fixed << std::setprecision(6)
              << "(" << cycle_length1 << ", " << cycle_length2 << ", " << cycle_length3 << ").\n";
    return false;
  }
  return true;
}

bool unsew_edge_width_repeatedly_in_unweighted_gmap() {
  LCC_for_GMap_2 lcc_gm;
  if (!get_data(lcc_gm)) {
    std::cout << "Fail unsew_edge_width_repeatedly_in_unweighted_gmap: Cannot locate file data/3torus.off\n";
    return false;
  }
  std::vector<unsigned int> cycle_lengths;
  unsigned int length;
  do {
    CGAL::Shortest_noncontractible_cycle<LCC_for_GMap_2> snc(lcc_gm);
    CGAL::Shortest_noncontractible_cycle<LCC_for_GMap_2>::Path cycle;
    snc.edge_width(cycle, &length);
    for (auto e : cycle) {
      if (e == NULL) {
        std::cout << "Fail unsew_edge_width_repeatedly_in_unweighted_gmap: NULL dart handle found in cycle\n";
        return false;
      }
      lcc_gm.unsew<2>(e);
    }
    cycle_lengths.push_back(length);
  } while (length != 0);
  for (int i = 1; i < cycle_lengths.size(); ++i)
    if (cycle_lengths[i] > cycle_lengths[i-1]) {
      std::cout << "Fail unsew_edge_width_repeatedly_in_unweighted_gmap: Edge width length decreases instead of increases\n";
      return false;
    }
  return true;
}

int main() {
  if (find_cycle_in_unweighted_cmap() && find_cycle_in_unweighted_gmap() &&
      edge_width_in_weighted_cmap_gmap_mesh() && unsew_edge_width_repeatedly_in_unweighted_gmap()) {
    std::cout << "All tests passed\n";
    return EXIT_SUCCESS;
  } else return EXIT_FAILURE;
}
