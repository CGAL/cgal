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
#include <unordered_map>

using LCC_3            =CGAL::Linear_cell_complex_for_generalized_map<2, 3>;
using Dart_descriptor      =LCC_3::Dart_descriptor;
using Dart_const_descriptor=LCC_3::Dart_const_descriptor;
using Dart_container   =std::vector<Dart_descriptor>;
using Point            =LCC_3::Point;
using Path_on_surface  =CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;
using CST              =CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;

struct Weight_functor
{
  Weight_functor(const LCC_3& lcc) : m_lcc(lcc) {}
  using Weight_t=double;
  Weight_t operator()(Dart_const_descriptor d) const
  {
    const Point& x=m_lcc.point(d);
    const Point& y=m_lcc.point(m_lcc.template alpha<0>(d));
    return CGAL::sqrt(CGAL::squared_distance(x, y));
  }
private:
  const LCC_3& m_lcc;
};

#ifdef CGAL_USE_BASIC_VIEWER

struct Draw_functor : public CGAL::DefaultDrawingFunctorLCC
{
  Draw_functor(LCC_3::size_type am1, LCC_3::size_type am2) : is_root(am1),
                                                             belong_to_cycle(am2)
  {}

  template<typename LCC>
  bool colored_vertex(const LCC& alcc, typename LCC::Dart_const_descriptor d) const
  { return alcc.is_marked(d, is_root); }

  template<typename LCC>
  CGAL::IO::Color vertex_color(const LCC& /* alcc */,
                           typename LCC::Dart_const_descriptor /* d */) const
  { return CGAL::IO::Color(0,255,0); }

  template<typename LCC>
  bool colored_edge(const LCC& alcc, typename LCC::Dart_const_descriptor d) const
  { return alcc.is_marked(d, belong_to_cycle); }

  template<typename LCC>
  CGAL::IO::Color edge_color(const LCC& /* alcc*/,
                         typename LCC::Dart_const_descriptor /* d */) const
  { return CGAL::IO::Color(0, 0, 255); }

  template<typename LCC>
  bool colored_face(const LCC& /* alcc */,
                    typename LCC::Dart_const_descriptor /* d */) const {return true;}

  template<typename LCC>
  CGAL::IO::Color face_color(const LCC& /* alcc */,
                         typename LCC::Dart_const_descriptor /* d */) const
  {return CGAL::IO::Color(211, 211, 211);}

  template<typename LCC>
  bool colored_volume(const LCC& /* alcc */,
                      typename LCC::Dart_const_descriptor /* d */) const { return false; }

  LCC_3::size_type is_root;
  LCC_3::size_type belong_to_cycle;
};

#endif // CGAL_USE_BASIC_VIEWER

int main(int argc, char* argv[])
{
  std::cout<<"Program unsew_edgewidth_repeatedly started."<<std::endl;
  std::string filename(argc==1?CGAL::data_file_path("meshes/double-torus-example.off"):argv[1]);

#ifdef CGAL_USE_BASIC_VIEWER
  bool draw=(argc<3?false:std::string(argv[2])=="-draw");
#endif // CGAL_USE_BASIC_VIEWER

  std::ifstream inp(filename);
  if (inp.fail())
  {
    std::cout<<"Cannot read file '"<<filename<<"'. Exiting program"<<std::endl;
    return EXIT_FAILURE;
  }

  LCC_3 lccoriginal, lcccopy;
  CGAL::load_off(lccoriginal, inp);
  std::cout<<"File '"<<filename<<"' loaded. Running the main program..."<<std::endl;

  std::unordered_map<Dart_descriptor, Dart_descriptor> origin_to_copy;
  lcccopy.copy(lccoriginal, &origin_to_copy, nullptr);

  LCC_3::size_type is_root=lccoriginal.get_new_mark();
  LCC_3::size_type belong_to_cycle=lccoriginal.get_new_mark();

  int loop=1;
  bool cycle_exist=true;
  Path_on_surface cycle(lcccopy);
  Weight_functor wf(lcccopy);
  do
  {
    std::cout<<"Finding #"<<loop++<<" edge-width:"<<std::endl;
    {
      CST cst(lcccopy);
      cycle=cst.compute_shortest_non_contractible_cycle(wf);
    }
    if (cycle.length()==0)
    { std::cout << "  Cannot find edge-width. Stop.\n"; cycle_exist=false; }
    else
    {
      LCC_3::size_type is_root_copy=lcccopy.get_new_mark();
      LCC_3::size_type belong_to_cycle_copy=lcccopy.get_new_mark();

      lcccopy.mark_cell<0>(cycle[0], is_root_copy);
      double cycle_length=0;
      for (std::size_t i=0; i<cycle.length(); ++i)
      {
        cycle_length+=wf(cycle[i]);
        if (!lcccopy.is_marked(cycle[i], belong_to_cycle_copy))
        { lcccopy.mark_cell<1>(cycle[i], belong_to_cycle_copy); }
      }

      for (auto d=lccoriginal.darts().begin(), dend=lccoriginal.darts().end();
           d!=dend; ++d)
      {
        if (lcccopy.is_marked(origin_to_copy[d], is_root_copy) &&
            !lccoriginal.is_marked(d, is_root))
        { lccoriginal.mark(d, is_root); }
        if (lcccopy.is_marked(origin_to_copy[d], belong_to_cycle_copy) &&
            !lccoriginal.is_marked(d, belong_to_cycle))
        { lccoriginal.mark(d, belong_to_cycle); }
        if (lcccopy.is_marked(origin_to_copy[d], belong_to_cycle_copy) &&
            !lcccopy.is_free<2>(origin_to_copy[d]))
        { lcccopy.unsew<2>(origin_to_copy[d]); }
      }
      lcccopy.close<2>();

      lcccopy.free_mark(belong_to_cycle_copy);
      lcccopy.free_mark(is_root_copy);

      std::cout<<"  Number of edges in cycle: "<<cycle.length()<<std::endl;
      std::cout<<"  Cycle length: "<<cycle_length<<std::endl;
      std::cout<<"  Root: "<<lcccopy.point(cycle[0])<<std::endl;
    }
  }
  while(cycle_exist);

#ifdef CGAL_USE_BASIC_VIEWER
  if (draw)
  {
    Draw_functor df(is_root, belong_to_cycle);
    CGAL::draw(lccoriginal, "Unsew edge width repeatdly", false, df);
  }
#endif // CGAL_USE_BASIC_VIEWER

  lccoriginal.free_mark(belong_to_cycle);
  lccoriginal.free_mark(is_root);

  return EXIT_SUCCESS;
}
