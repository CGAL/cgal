#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/draw_linear_cell_complex.h>

using LCC_3            =CGAL::Linear_cell_complex_for_combinatorial_map<2, 3>;
using CST              =CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;
using Path_on_surface  =CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;
using Dart_const_handle=LCC_3::Dart_const_handle;

struct Draw_functor : public CGAL::DefaultDrawingFunctorLCC
{
  Draw_functor(LCC_3::size_type amark1, LCC_3::size_type amark2) :
    m_vertex_mark(amark1), m_face_mark(amark2)
  {}

  template<typename LCC>
  bool colored_vertex(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return alcc.is_marked(dh, m_vertex_mark); }

  template<typename LCC>
  CGAL::Color vertex_color(const LCC& /* alcc */,
                           typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0, 255, 0); }

  template<typename LCC>
  bool colored_edge(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return false; }

  template<typename LCC>
  CGAL::Color edge_color(const LCC& /* alcc*/,
                         typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0, 0, 255); }

  template<typename LCC>
  bool colored_face(const LCC& /* alcc */,
                    typename LCC::Dart_const_handle /* dh */) const
  {return true;}

  template<typename LCC>
  CGAL::Color face_color(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return alcc.is_marked(dh, m_face_mark)?CGAL::Color(255, 0, 0)
                                          :CGAL::Color(211, 211, 211); }

  template<typename LCC>
  bool colored_volume(const LCC& /* alcc */,
                      typename LCC::Dart_const_handle /* dh */) const
  { return false; }

  LCC_3::size_type m_vertex_mark, m_face_mark;
};

int main(int argc, char* argv[])
{
  std::cout<<"Program facewidth_on_unweighted_map started."<<std::endl;
  std::string filename(argc==1?"data/double-torus.off":argv[1]);
  bool draw=(argc<3?false:std::string(argv[2])=="-draw");

  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    return EXIT_FAILURE;
  }
  LCC_3 lcc;
  CGAL::load_off(lcc, inp);
  std::cout<<"File '"<<filename<<"' loaded. Finding the facewidth..."<<std::endl;

  CST cst(lcc, true);
  std::vector<Dart_const_handle> cycle=cst.compute_facewidth(true);

  if (cycle.size()==0)
  { std::cout<<"  Cannot find such cycle. Stop."<<std::endl; }
  else
  {
    std::cout<<"  Number of faces: "<<cycle.size()/2<<std::endl;

    if (draw)
    {
      LCC_3::size_type vertex_mark = lcc.get_new_mark();
      LCC_3::size_type face_mark = lcc.get_new_mark();
      for (int i=0; i<cycle.size(); ++i)
      {
        if (i%2==0)
        { // Color the vertex
          lcc.mark_cell<0>(cycle[i], vertex_mark);
        }
        else
        { // Color the face
          lcc.mark_cell<2>(cycle[i], face_mark);
        }
      }

      Draw_functor df(vertex_mark, face_mark);
      CGAL::draw(lcc, "Face width", false, df);

      lcc.free_mark(vertex_mark);
      lcc.free_mark(face_mark);
    }
  }

  return EXIT_SUCCESS;
}
