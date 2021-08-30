#include <CGAL/Generalized_map.h>
#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/squared_distance_3.h>
#include <iostream>
#include <cstdlib>
#include <CGAL/Curves_on_surface_topology.h>

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
using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Point = Kernel::Point_3;
using Surface_mesh = CGAL::Surface_mesh<Point>;
using Polyhedron = CGAL::Polyhedron_3<Kernel>;

struct Weight_functor_for_GM {
  using Weight_t = unsigned int;
  Weight_functor_for_GM(const GMap_2& gm, GMap_2::size_type amark) : m_gm(gm), m_mark(amark) {}
  unsigned int operator() (GMap_2::Dart_const_handle dh) const {
    if (m_gm.is_marked(dh, m_mark)) return 3;
    else return 4;
  }
private:
  const GMap_2&     m_gm;
  GMap_2::size_type m_mark;
};

struct Weight_functor_for_SM {
  using Weight_t = double;
  Weight_functor_for_SM(const Surface_mesh& mesh) : m_mesh(mesh) {}
  double operator()(Surface_mesh::Halfedge_index he) const {
    Point A = m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 0));
    Point B = m_mesh.point(m_mesh.vertex(m_mesh.edge(he), 1));
    return CGAL::sqrt(CGAL::squared_distance(A, B));
  }
private:
  const Surface_mesh& m_mesh;
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

bool get_data(Polyhedron& sm) {
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

bool find_cycle_in_unweighted_cmap_and_polyhedron() {
  LCC_for_CMap_2 lcc;
  if (!get_data(lcc)) {
    std::cerr << "Fail find_cycle_in_unweighted_cmap_and_polyhedron: Cannot locate file data/3torus.off\n";
    return false;
  }
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_for_CMap_2> cst1(lcc);
  LCC_for_CMap_2::Dart_handle root1 = lcc.darts().begin();
  Point R = lcc.point_of_vertex_attribute(lcc.vertex_attribute(root1));
  CGAL::Surface_mesh_topology::Path_on_surface<LCC_for_CMap_2> cycle1 = cst1.compute_shortest_non_contractible_cycle_with_base_point(root1);
  for (std::size_t i = 0; i < cycle1.length(); ++i) {
    auto e = cycle1[i];
    if (e == nullptr) {
      std::cerr << "Fail find_cycle_in_unweighted_cmap_and_polyhedron: NULL dart handle found in cycle1\n";
      return false;
    }
  }
  Polyhedron p;
  if (!get_data(p)) {
    std::cerr << "Fail find_cycle_in_unweighted_cmap_and_polyhedron: Cannot locate file data/3torus.off\n";
    return false;
  }
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Polyhedron> cst2(p);
  boost::graph_traits<Polyhedron>::halfedge_iterator root2 = p.halfedges_begin(), endit = p.halfedges_end();
  for (; root2 != endit; ++root2) {
    if ((*root2)->vertex()->point() == R) break;
  }
  if (root2 == endit) {
    std::cerr << "Fail find_cycle_in_unweighted_cmap_and_polyhedron: Cannot find CMap's root in the Polyhedron\n";
    return false;
  }
  CGAL::Surface_mesh_topology::Path_on_surface<Polyhedron> cycle2 = cst2.compute_shortest_non_contractible_cycle_with_base_point(*root2);
  for (std::size_t i = 0; i < cycle2.length(); ++i) {
    auto e = cycle2[i];
    if (e == nullptr) {
      std::cerr << "Fail find_cycle_in_unweighted_cmap_and_polyhedron: NULL dart handle found in cycle\n";
      return false;
    }
  }
  if (cycle1.length() != cycle2.length()) {
    std::cerr << "Fail find_cycle_in_unweighted_cmap_and_polyhedron: Inconsistency in cycle length\n";
    return false;
  }
  return true;
}

bool edge_width_in_unweighted_polyhedron() {
  Polyhedron p;
  if (!get_data(p)) {
    std::cerr << "Fail edge_width_in_unweighted_polyhedron: Cannot locate file data/3torus.off\n";
    return false;
  }
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Polyhedron> cst(p);
  CGAL::Surface_mesh_topology::Path_on_surface<Polyhedron> cycle = cst.compute_edge_width();
  for (std::size_t i = 0; i < cycle.length(); ++i) {
    auto e = cycle[i];
    if (e == nullptr) {
      std::cerr << "Fail edge_width_in_unweighted_polyhedron: NULL dart handle found in cycle\n";
      return false;
    }
  }
  return true;
}

bool find_cycle_in_nonorientable_gmap() { // Make a non-oriented case here
  // Sewing the Petersen graph embedded in a Klein bottle surface
  GMap_2 gm;
  std::vector<GMap_2::Dart_handle> faces;
  for (int i = 0; i < 6; ++i) faces.push_back(gm.make_combinatorial_polygon(5));
  gm.sew<2>(faces[0], faces[1]); // 1-2
  gm.sew<2>(gm.alpha<1>(faces[1]), gm.alpha<1>(faces[2])); // 1-6
  gm.sew<2>(gm.alpha<1>(faces[0]), faces[2]); // 1-5
  gm.sew<2>(gm.next(faces[0]), gm.next(faces[4])); // 2-3
  gm.sew<2>(gm.next(faces[1]), gm.alpha<0>(faces[4])); // 2-7
  gm.sew<2>(gm.alpha<1,0,1>(faces[1]), gm.alpha<0,1,0,1,0>(faces[3])); // 6-9
  gm.sew<2>(gm.alpha<1,0,1>(faces[2]), gm.alpha<1,0,1,0>(faces[3])); // 6-8
  gm.sew<2>(gm.alpha<1,0,1>(faces[0]), gm.alpha<1,0,1,0,1>(faces[5])); // 5-4
  gm.sew<2>(gm.next(faces[2]), gm.alpha<1,0,1,0>(faces[5])); // 5-10
  gm.sew<2>(gm.alpha<0,1,0,1>(faces[0]), faces[3]); // 3-4
  gm.sew<2>(gm.alpha<0,1,0,1>(faces[1]), faces[5]); // 7-9
  gm.sew<2>(gm.alpha<0,1,0,1,0>(faces[2]), gm.alpha<1,0,1,0>(faces[4])); // 8-10
  gm.sew<2>(gm.alpha<1>(faces[4]), gm.alpha<1>(faces[5])); // 7-10
  gm.sew<2>(gm.alpha<0,1,0,1>(faces[4]), gm.alpha<1>(faces[3])); // 3-8
  gm.sew<2>(gm.alpha<0,1>(faces[5]), gm.alpha<0,1,0>(faces[3])); // 9-4

  // gm.display_characteristics(std::cerr);
  // std::cerr << '\n';

  GMap_2::size_type chosen_cycle = gm.get_new_mark(), smallest_edge = gm.get_new_mark();
  gm.mark_cell<1>(gm.alpha<0,1>(faces[5]), smallest_edge); // 9-4

  Weight_functor_for_GM wf (gm, smallest_edge);
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<GMap_2> cst(gm);

  CGAL::Surface_mesh_topology::Path_on_surface<GMap_2> cycle = cst.compute_shortest_non_contractible_cycle_with_base_point(faces[0]);

  gm.mark_cell<1>(gm.alpha<1>(faces[1]), chosen_cycle); // 1-6
  gm.mark_cell<1>(gm.alpha<1,0,1>(faces[1]), chosen_cycle); // 6-9
  gm.mark_cell<1>(gm.alpha<0,1>(faces[5]), chosen_cycle); // 9-4
  gm.mark_cell<1>(gm.alpha<1,0,1,0>(faces[0]), chosen_cycle); // 4-5
  gm.mark_cell<1>(gm.alpha<0>(faces[2]), chosen_cycle); // 5-1

  unsigned int cycle_length = 0;
  for (std::size_t i = 0; i < cycle.length(); ++i) {
    cycle_length += wf(cycle[i]);
    auto dh = cycle[i];
    if (!gm.is_marked(dh, chosen_cycle)) {
      std::cerr << "Fail find_cycle_in_nonorientable_gmap: Cycle found is not the same as expected cycle.\n";
      return false;
    }
  }

  if (cycle_length != 19) {
    std::cerr << "Fail find_cycle_in_nonorientable_gmap: Cycle length (" << cycle_length << ") is not as expected (should be 19).\n";
    return false;
  }
  gm.free_mark(chosen_cycle);
  gm.free_mark(smallest_edge);
  return true;
}

bool edge_width_in_weighted_cmap_gmap_mesh() {
  LCC_for_CMap_2 lcc_cm;
  LCC_for_GMap_2 lcc_gm;
  Surface_mesh sm;
  if (!get_data(lcc_cm) || !get_data(lcc_gm) || !get_data(sm)) {
    std::cerr << "Fail edge_width_in_weighted_cmap_gmap_mesh: Cannot locate file data/3torus.off\n";
    return false;
  }
  Weight_functor_for_LCC<LCC_for_CMap_2> wf_cm(lcc_cm);
  Weight_functor_for_LCC<LCC_for_GMap_2> wf_gm(lcc_gm);
  Weight_functor_for_SM wf_sm(sm);
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_for_CMap_2> cst1(lcc_cm);
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_for_GMap_2> cst2(lcc_gm);
  CGAL::Surface_mesh_topology::Curves_on_surface_topology<Surface_mesh> cst3(sm);
  CGAL::Surface_mesh_topology::Path_on_surface<LCC_for_CMap_2> cycle1 = cst1.compute_shortest_non_contractible_cycle(wf_cm);
  CGAL::Surface_mesh_topology::Path_on_surface<LCC_for_GMap_2> cycle2 = cst2.compute_shortest_non_contractible_cycle(wf_gm);
  CGAL::Surface_mesh_topology::Path_on_surface<Surface_mesh> cycle3 = cst3.compute_shortest_non_contractible_cycle(wf_sm);

  if (cycle1.length()!=cycle2.length() || cycle1.length()!=cycle3.length())
  {
    std::cerr << "Fail edge_width_in_weighted_cmap_gmap_mesh: Inconsistency in number of edges of the edge-width "
              << "(" << cycle1.length() << ", " << cycle2.length() << ", " << cycle3.length() << ").\n";
    return false;
  }

  if (!cycle1.is_valid())
  {
    std::cerr << "Fail edge_width_in_weighted_cmap_gmap_mesh: cycle1 is ill-formed\n";
    return false;
  }
  if (!cycle2.is_valid())
  {
    std::cerr << "Fail edge_width_in_weighted_cmap_gmap_mesh: cycle2 is ill-formed\n";
    return false;
  }
  if (!cycle3.is_valid())
  {
    std::cerr << "Fail edge_width_in_weighted_cmap_gmap_mesh: cycle3 is ill-formed\n";
    return false;
  }

  double cycle_length1=0, cycle_length2=0, cycle_length3=0;
  std::vector<Point> v1, v2, v3;

  for (std::size_t i=0; i<cycle1.length(); ++i)
  {
    auto e=cycle1.get_ith_real_dart(i);
    cycle_length1+=wf_cm(e);
    v1.push_back(lcc_cm.point(e));
  }
  v1.push_back(lcc_cm.point(cycle1.get_ith_real_dart(0)));

  std::size_t i=0;
  while (lcc_gm.point(cycle2.get_ith_real_dart(i))!=v1[0])
  {
    if (i==cycle2.length())
    {
      std::cerr<<"Fail edge_width_in_weighted_cmap_gmap_mesh: cycle 2: can not find root "<<v1[0]<<std::endl;
      return false;
    }
    ++i;
  }
  for (std::size_t j=0; j<cycle2.length(); ++j, i=cycle2.next_index(i))
  {
    auto e=cycle2.get_ith_real_dart(i);
    cycle_length2 += wf_gm(e);
    v2.push_back(lcc_gm.point(e));
  }
  v2.push_back(lcc_gm.point(cycle2.get_ith_real_dart(i)));

  i=0;
  while (sm.point(sm.source(cycle3.get_ith_real_dart(i)))!=v1[0])
  {
    if (i==cycle3.length())
    {
      std::cerr<<"Fail edge_width_in_weighted_cmap_gmap_mesh: cycle 3: can not find root "<<v1[0]<<std::endl;
      return false;
    }
    ++i;
  }
  for (std::size_t j=0; j<cycle3.length(); ++j, i=cycle3.next_index(i))
  {
    auto e=cycle3.get_ith_real_dart(i);
    cycle_length3+=wf_sm(e);
    v3.push_back(sm.point(sm.source(e)));
  }
  v3.push_back(sm.point(sm.source(cycle3.get_ith_real_dart(i))));

  int v2orientation=0, v3orientation=0;
  bool same=true;
  for (std::size_t i=0; i<v1.size(); ++i)
  {
    std::size_t j=v1.size()-i-1;
    if (v2orientation==0)
    {
      if (v1[i]==v2[i] && v1[i]!=v2[j])      { v2orientation=1; }
      else if (v1[i]!=v2[i] && v1[i]==v2[j]) { v2orientation=2; }
      else if (v1[i]!=v2[i] && v1[i]!=v2[j]) { same=false; }
    }
    else if (v2orientation==1)
    { if (v1[i]!=v2[i]) { same=false; } }
    else if (v2orientation==2)
    { if (v1[i]!=v2[j]) { same=false; } }

    if (v3orientation==0)
    {
      if (v1[i]==v3[i] && v1[i]!=v3[j]) { v3orientation=1; }
      else if (v1[i]!=v3[i] && v1[i]==v3[j]) { v3orientation=2; }
      else if (v1[i]!=v3[i] && v1[i]!=v3[j]) { same=false; }
    }
    else if (v3orientation==1)
    { if (v1[i]!=v3[i]) { same=false; } }
    else if (v3orientation==2)
    { if (v1[i]!=v3[j]) { same=false; } }

    if (!same)
    {
      std::cerr << "Fail edge_width_in_weighted_cmap_gmap_mesh: Inconsistency in the vertex ordering.\n";
      for (std::size_t i=0; i<v1.size(); ++i) std::cerr<<v1[i]<<", ";
      std::cerr<<'\n';
      for (std::size_t i=0; i<v2.size(); ++i) std::cerr<<v2[i]<<", ";
      std::cerr<<'\n';
      for (std::size_t i=0; i<v3.size(); ++i) std::cerr<<v3[i]<<", ";
      std::cerr<<'\n';
      return false;
    }
  }
  if (v1[0]!=v1.back())
  {
    std::cerr << "Fail edge_width_in_weighted_cmap_gmap_mesh: The path is not a cycle.\n";
    return false;
  }
  if (cycle_length1-cycle_length2 > 1e-5 || cycle_length1-cycle_length3 > 1e-5)
  {
    std::cerr << "Fail edge_width_in_weighted_cmap_gmap_mesh: Inconsistency in the edge-width length"
              << std::fixed << std::setprecision(6)
              << "(" << cycle_length1 << ", " << cycle_length2 << ", " << cycle_length3 << ").\n";
    return false;
  }
  return true;
}

bool unsew_edge_width_repeatedly_in_unweighted_gmap() {
  LCC_for_GMap_2 lcc_gm;
  if (!get_data(lcc_gm)) {
    std::cerr << "Fail unsew_edge_width_repeatedly_in_unweighted_gmap: Cannot locate file data/3torus.off\n";
    return false;
  }
  std::vector<std::size_t> cycle_lengths;
  std::size_t length;
  do {
    CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_for_GMap_2> cst(lcc_gm);
    CGAL::Surface_mesh_topology::Path_on_surface<LCC_for_GMap_2> cycle = cst.compute_edge_width();
    length = cycle.length();
    LCC_for_GMap_2::size_type belong_to_cycle = lcc_gm.get_new_mark();
    for (std::size_t i = 0; i < cycle.length(); ++i) {
      auto e = cycle[i];
      if (e == nullptr) {
        std::cerr << "Fail unsew_edge_width_repeatedly_in_unweighted_gmap: NULL dart handle found in cycle\n";
        return false;
      }
      lcc_gm.mark(e, belong_to_cycle);
    }
    for (auto dh = lcc_gm.darts().begin(), dhend = lcc_gm.darts().end(); dh != dhend; ++dh)
      if (lcc_gm.is_marked(dh, belong_to_cycle))
        lcc_gm.unsew<2>(dh);
    lcc_gm.close<2>();
    cycle_lengths.push_back(length);
    lcc_gm.free_mark(belong_to_cycle);
  } while (length != 0);
  for (std::size_t i = 1; i < cycle_lengths.size(); ++i)
    if (cycle_lengths[i] > cycle_lengths[i-1]) {
      std::cerr << "Fail unsew_edge_width_repeatedly_in_unweighted_gmap: Edge width length decreases instead of increases\n";
      return false;
    }
  return true;
}

int main() {
  bool res=find_cycle_in_unweighted_cmap_and_polyhedron() &&
      edge_width_in_unweighted_polyhedron() &&
      // find_cycle_in_nonorientable_gmap() &&
      edge_width_in_weighted_cmap_gmap_mesh() &&
      unsew_edge_width_repeatedly_in_unweighted_gmap();

  if (res)
  {
    std::cout << "All tests passed\n";
    return EXIT_SUCCESS;
  }
  return EXIT_FAILURE;
}
