#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/IO/Color.h>
#include <CGAL/draw_linear_cell_complex.h>
#include <CGAL/squared_distance_3.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <unordered_set>
#define CGAL_USE_BASIC_VIEWER 1

using LCC_3 = CGAL::Linear_cell_complex_for_generalized_map<2, 3>;
using Dart_handle = LCC_3::Dart_handle;
using Dart_const_handle = LCC_3::Dart_const_handle;
using Dart_container = std::vector<Dart_handle>;
using Point = LCC_3::Point;
using Path_on_surface = CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;

struct Weight_functor {
  Weight_functor(const LCC_3& lcc) : m_lcc(lcc) { }
  using Weight_t = double;
  Weight_t operator()(Dart_const_handle dh) const {
    const Point& x = m_lcc.point(dh);
    const Point& y = m_lcc.point(m_lcc.template alpha<0>(dh));
    return CGAL::sqrt(CGAL::squared_distance(x, y));
  }
private:
  const LCC_3& m_lcc;
};

using CST = CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;

struct Draw_functor : public CGAL::DefaultDrawingFunctorLCC {
  Draw_functor(LCC_3::size_type am1, LCC_3::size_type am2) : is_root(am1),
                                                             belong_to_cycle(am2)
  {}

  template<typename LCC>
  bool colored_vertex(const LCC& alcc, typename LCC::Dart_const_handle dh) const 
  { return alcc.is_marked(dh, is_root); }
  
  template<typename LCC>
  CGAL::Color vertex_color(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0,255,0); }

  template<typename LCC>
  bool colored_edge(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return alcc.is_marked(dh, belong_to_cycle); }

  template<typename LCC>
  CGAL::Color edge_color(const LCC& /* alcc*/, typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0, 0, 255); }

  template<typename LCC>
  bool colored_face(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const {return true;}

  template<typename LCC>
  CGAL::Color face_color(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const
  {return CGAL::Color(211, 211, 211);}

  template<typename LCC>
  bool colored_volume(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const { return false; }
  
  LCC_3::size_type is_root;
  LCC_3::size_type belong_to_cycle;
};

int main(int argc, char* argv[]) {
  std::cout << "Program unsew_edgewidth_repeatedly started.\n";
  LCC_3 lccoriginal, lcccopy;
  std::ifstream inp;
  if (argc == 1) inp = std::ifstream("data/double-torus-example.off");
  else inp = std::ifstream(argv[1]);
  CGAL::load_off(lccoriginal, inp);

  boost::unordered_map<Dart_handle, Dart_handle> copy_to_origin;
  lcccopy.copy(lccoriginal, NULL, &copy_to_origin);

  LCC_3::size_type is_root=lccoriginal.get_new_mark();
  LCC_3::size_type belong_to_cycle=lccoriginal.get_new_mark();
  Draw_functor df(is_root, belong_to_cycle);

  std::cout << "File loaded. Running the main program...\n";
  for (int loop = 1; ; ++loop) {
    std::cout << "Finding #" << loop << " edge-width:\n";
    Weight_functor wf(lcccopy);
    CST cst(lcccopy);
    Path_on_surface cycle = cst.compute_edgewidth(wf);
    if (cycle.length() == 0) {
      std::cout << "  Cannot find edge-width. Stop.\n";
      break;
    }

    LCC_3::size_type is_root_copy = lcccopy.get_new_mark();
    LCC_3::size_type belong_to_cycle_copy = lcccopy.get_new_mark();

    lcccopy.mark_cell<0>(cycle[0], is_root_copy);
    double x = 0;
    for (int i = 0; i < cycle.length(); ++i) {
      x += wf(cycle[i]);
      if (!lcccopy.is_marked(cycle[i], belong_to_cycle_copy))
        lcccopy.mark_cell<1>(cycle[i], belong_to_cycle_copy);
    }

    for (auto dh = lcccopy.darts().begin(), dhend = lcccopy.darts().end(); dh != dhend; ++dh) {
      if (lcccopy.is_marked(dh, is_root_copy) && !lccoriginal.is_marked(copy_to_origin[dh], is_root))
        lccoriginal.mark(copy_to_origin[dh], is_root);
      if (lcccopy.is_marked(dh, belong_to_cycle_copy) && !lccoriginal.is_marked(copy_to_origin[dh], belong_to_cycle))
        lccoriginal.mark(copy_to_origin[dh], belong_to_cycle);
      if (lcccopy.is_marked(dh, belong_to_cycle_copy) && !lcccopy.is_free<2>(dh))
        lcccopy.unsew<2>(dh);
    }

    lcccopy.free_mark(belong_to_cycle_copy);
    lcccopy.free_mark(is_root_copy);
    
    std::cout << "  Number of edges in cycle: " << cycle.length() << std::endl;
    std::cout << "  Cycle length: " << x << std::endl;
    std::cout << "  Root: " << lcccopy.point_of_vertex_attribute(lcccopy.vertex_attribute(cycle[0])) << std::endl;

  }


  CGAL::draw(lccoriginal, "Hello", false, df);

  lccoriginal.free_mark(belong_to_cycle);
  lccoriginal.free_mark(is_root);  
}
