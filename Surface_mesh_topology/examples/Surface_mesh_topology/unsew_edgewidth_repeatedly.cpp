#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Shortest_noncontractible_cycle.h>
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

using SNC = CGAL::Surface_mesh_topology::Shortest_noncontractible_cycle<LCC_3, Weight_functor>;

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
  if (argc == 1) inp = std::ifstream("../../examples/Surface_mesh_topology/data/3torus.off");
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
    SNC snc(lcccopy, wf);
    SNC::Path cycle;
    SNC::Distance_type x;
    snc.edge_width(cycle, &x);
    if (cycle.size() == 0) {
      std::cout << "  Cannot find edge-width. Stop.\n";
      break;
    }

    for (LCC_3::Dart_of_cell_range<0>::iterator it=lcccopy.darts_of_cell<0>(cycle[0]).begin(),
           itend=lcccopy.darts_of_cell<0>(cycle[0]).end(); it!=itend; ++it)
    {
      if (copy_to_origin.count(it)>0 &&
          !lccoriginal.is_marked(copy_to_origin[it], is_root))
      {
        lccoriginal.mark_cell<0>(copy_to_origin[it], is_root);
        break;
      }
    }
    
    for (auto e : cycle)
    {
      for (auto it=lcccopy.darts_of_cell<1>(e).begin(),
             itend=lcccopy.darts_of_cell<1>(e).end(); it!=itend; ++it)
      {
        if (copy_to_origin.count(it)>0 &&
            !lccoriginal.is_marked(copy_to_origin[it], belong_to_cycle))
          {
            lccoriginal.mark_cell<1>(copy_to_origin[it], belong_to_cycle);
            break;
          }
      }
      if (!lcccopy.is_free<2>(e)) 
      { lcccopy.unsew<2>(e); }
    }
    
    std::cout << "  Number of edges in cycle: " << cycle.size() << std::endl;
    std::cout << "  Cycle length: " << x << std::endl;
    std::cout << "  Root: " << lcccopy.point_of_vertex_attribute(lcccopy.vertex_attribute(cycle[0])) << std::endl;
  }


  CGAL::draw(lccoriginal, "Hello", false, df);

  lccoriginal.free_mark(belong_to_cycle);
  lccoriginal.free_mark(is_root);  
}
