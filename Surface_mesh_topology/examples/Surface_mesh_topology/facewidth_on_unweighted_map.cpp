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
using Dart_handle = LCC_3::Dart_handle;

struct Draw_functor : public CGAL::DefaultDrawingFunctorLCC {
  Draw_functor(LCC_3::size_type amark1, LCC_3::size_type amark2) : m_vertex_mark(amark1), m_face_mark(amark2)
  {}

  template<typename LCC>
  bool colored_vertex(const LCC& alcc, typename LCC::Dart_const_handle dh) const 
  { return alcc.is_marked(dh, m_vertex_mark); }
  
  template<typename LCC>
  CGAL::Color vertex_color(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0, 255, 0); }

  template<typename LCC>
  bool colored_edge(const LCC& alcc, typename LCC::Dart_const_handle dh) const { return false; }

  template<typename LCC>
  CGAL::Color edge_color(const LCC& /* alcc*/, typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0, 0, 255); }

  template<typename LCC>
  bool colored_face(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const {return true;}

  template<typename LCC>
  CGAL::Color face_color(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return alcc.is_marked(dh, m_face_mark) ? CGAL::Color(255, 0, 0) : CGAL::Color(211, 211, 211); }

  template<typename LCC>
  bool colored_volume(const LCC& /* alcc */, typename LCC::Dart_const_handle /* dh */) const { return false; }
  
  LCC_3::size_type m_vertex_mark, m_face_mark;
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
  
  CST cst(lcc);
  std::vector<Dart_handle> cycle = cst.compute_facewidth();
  std::cout << "Finding the facewidth...\n";

  if (cycle.size() == 0) {
    std::cout << "  Cannot find such cycle. Stop.\n";
    return 0;
  }
  LCC_3::size_type vertex_mark = lcc.get_new_mark();
  LCC_3::size_type face_mark = lcc.get_new_mark();
  for (int i = 0; i < cycle.size(); ++i)
  {
    if (i % 2 == 0)
    {
      // Color the edges of the face
      for (auto dh = lcc.darts_of_cell<0>(cycle[i]).begin(),
                dhend = lcc.darts_of_cell<0>(cycle[i]).end(); dh != dhend; ++dh)
      {
        if (!lcc.is_marked(dh, vertex_mark))
          lcc.mark(dh, vertex_mark);
      }
    } 
    else 
    {
      // Color the face
      for (auto dh = lcc.darts_of_cell<2>(cycle[i]).begin(),
                dhend = lcc.darts_of_cell<2>(cycle[i]).end(); dh != dhend; ++dh)
      {
        if (!lcc.is_marked(dh, face_mark))
          lcc.mark(dh, face_mark);
      }
    }
  }

  std::cout << "  Number of faces: " << cycle.size()/2 << std::endl;

  Draw_functor df(vertex_mark, face_mark);
  CGAL::draw(lcc, "Hello", false, df);

  lcc.free_mark(vertex_mark);
  lcc.free_mark(face_mark);
}
