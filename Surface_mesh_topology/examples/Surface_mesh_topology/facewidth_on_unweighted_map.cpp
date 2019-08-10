#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_linear_cell_complex.h>
#define CGAL_USE_BASIC_VIEWER 1


using LCC_3 = CGAL::Linear_cell_complex_for_combinatorial_map<2, 3>;
using CST = CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;
using Path_on_surface = CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;

struct Draw_functor : public CGAL::DefaultDrawingFunctorLCC {
  Draw_functor(LCC_3::size_type amark1, LCC_3::size_type amark2) : m_belong_to_cycle(amark1), m_belong_to_facewidth(amark2)
  {}

  template<typename LCC>
  bool colored_vertex(const LCC& alcc, typename LCC::Dart_const_handle dh) const 
  { return false; }
  
  template<typename LCC>
  CGAL::Color vertex_color(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0,255,0); }

  template<typename LCC>
  bool colored_edge(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return alcc.is_marked(dh, m_belong_to_cycle); }

  template<typename LCC>
  CGAL::Color edge_color(const LCC& /* alcc*/, typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0, 0, 255); }

  template<typename LCC>
  bool colored_face(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const {return true;}

  template<typename LCC>
  CGAL::Color face_color(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  // {return CGAL::Color(211, 211, 211);}
  {return alcc.is_marked(dh, m_belong_to_facewidth) ? CGAL::Color(255, 0, 0) : CGAL::Color(211, 211, 211);}

  template<typename LCC>
  bool colored_volume(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const { return false; }
  
  LCC_3::size_type m_belong_to_cycle, m_belong_to_facewidth;
};

LCC_3 lcc;

int main(int argc, char* argv[])
{
  std::cout << "Program facewidth_on_unweighted_map started.\n";
  std::ifstream inp;
  if (argc == 1) inp = std::ifstream("data/double-torus-example.off");
  else inp = std::ifstream(argv[1]);
  if (inp.fail()) {
    std::cout << "Cannot read file. Exiting program\n";
    return EXIT_FAILURE;
  }
  CGAL::load_off(lcc, inp);
  std::cout << "File loaded. Running the main program...\n";
  
  CST                cst(lcc);
  Path_on_surface cycle = cst.compute_facewidth();
  std::cout << "Finding the facewidth...\n";

  if (cycle.length() == 0) {
    std::cout << "  Cannot find such cycle. Stop.\n";
    return 0;
  }
  LCC_3::size_type belong_to_cycle = lcc.get_new_mark();
  LCC_3::size_type belong_to_facewidth = lcc.get_new_mark();
  int nb_faces = 1;
  for (int i = 0; i < cycle.length(); ++i) {
    for (auto dh = lcc.one_dart_per_incident_cell<1,2>(cycle[i]).begin(),
              dhend = lcc.one_dart_per_incident_cell<1,2>(cycle[i]).end(); dh != dhend; ++dh)
    {
      if (!lcc.is_marked(dh, belong_to_cycle))
        lcc.mark_cell<1>(dh, belong_to_cycle);
    }
    for (auto dh = lcc.darts_of_cell<2>(cycle[i]).begin(),
              dhend = lcc.darts_of_cell<2>(cycle[i]).end(); dh != dhend; ++dh)
    {
      if (!lcc.is_marked(dh, belong_to_facewidth))
        lcc.mark(dh, belong_to_facewidth);
    }
    if (i > 0 && !lcc.belong_to_same_cell<2>(cycle[i], cycle[i-1])) ++nb_faces;
  }

  std::cout << "  Number of edges in cycle: " << nb_faces << std::endl;
  // for (int i = 0; i < cycle.length(); ++i)
  //   std::cout << lcc.point(cycle[i]) << " ---- ";
  // std::cout << lcc.point(lcc.next(cycle.back()));
  Draw_functor df(belong_to_cycle, belong_to_facewidth);
  CGAL::draw(lcc, "Hello", false, df);

  lcc.free_mark(belong_to_cycle);
  lcc.free_mark(belong_to_facewidth);
}
