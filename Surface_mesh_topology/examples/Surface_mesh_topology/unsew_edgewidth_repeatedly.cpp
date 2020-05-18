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

using LCC_3            =CGAL::Linear_cell_complex_for_generalized_map<2, 3>;
using Dart_handle      =LCC_3::Dart_handle;
using Dart_const_handle=LCC_3::Dart_const_handle;
using Dart_container   =std::vector<Dart_handle>;
using Point            =LCC_3::Point;
using Path_on_surface  =CGAL::Surface_mesh_topology::Path_on_surface<LCC_3>;
using CST              =CGAL::Surface_mesh_topology::Curves_on_surface_topology<LCC_3>;

struct Weight_functor
{
  Weight_functor(const LCC_3& lcc) : m_lcc(lcc) {}
  using Weight_t=double;
  Weight_t operator()(Dart_const_handle dh) const
  {
    const Point& x=m_lcc.point(dh);
    const Point& y=m_lcc.point(m_lcc.template alpha<0>(dh));
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
  bool colored_vertex(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return alcc.is_marked(dh, is_root); }

  template<typename LCC>
  CGAL::Color vertex_color(const LCC& /* alcc */,
                           typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0,255,0); }

  template<typename LCC>
  bool colored_edge(const LCC& alcc, typename LCC::Dart_const_handle dh) const
  { return alcc.is_marked(dh, belong_to_cycle); }

  template<typename LCC>
  CGAL::Color edge_color(const LCC& /* alcc*/,
                         typename LCC::Dart_const_handle /* dh */) const
  { return CGAL::Color(0, 0, 255); }

  template<typename LCC>
  bool colored_face(const LCC& /* alcc */,
                    typename LCC::Dart_const_handle /* dh */) const {return true;}

  template<typename LCC>
  CGAL::Color face_color(const LCC& /* alcc */,
                         typename LCC::Dart_const_handle /* dh */) const
  {return CGAL::Color(211, 211, 211);}

  template<typename LCC>
  bool colored_volume(const LCC& /* alcc */,
                      typename LCC::Dart_const_handle /* dh */) const { return false; }

  LCC_3::size_type is_root;
  LCC_3::size_type belong_to_cycle;
};

#endif // CGAL_USE_BASIC_VIEWER

int main(int argc, char* argv[])
{
  std::cout<<"Program unsew_edgewidth_repeatedly started."<<std::endl;
  std::string filename(argc==1?"data/double-torus.off":argv[1]);

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

  boost::unordered_map<Dart_handle, Dart_handle> origin_to_copy;
  lcccopy.copy(lccoriginal, &origin_to_copy, NULL);

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

      for (auto dh=lccoriginal.darts().begin(), dhend=lccoriginal.darts().end();
           dh!=dhend; ++dh)
      {
        if (lcccopy.is_marked(origin_to_copy[dh], is_root_copy) &&
            !lccoriginal.is_marked(dh, is_root))
        { lccoriginal.mark(dh, is_root); }
        if (lcccopy.is_marked(origin_to_copy[dh], belong_to_cycle_copy) &&
            !lccoriginal.is_marked(dh, belong_to_cycle))
        { lccoriginal.mark(dh, belong_to_cycle); }
        if (lcccopy.is_marked(origin_to_copy[dh], belong_to_cycle_copy) &&
            !lcccopy.is_free<2>(origin_to_copy[dh]))
        { lcccopy.unsew<2>(origin_to_copy[dh]); }
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
