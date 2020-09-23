#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Curves_on_surface_topology.h>
#include <CGAL/Surface_mesh_topology/internal/Path_generators.h>
#include <vector>
#include <sstream>
#include <tuple>

#include "Creation_of_test_cases_for_paths.h"

/* If you want to use a viewer, you can use qglviewer. */
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/draw_face_graph_with_paths.h>
#endif

struct MyItems
{
  template <class CMap>
  struct Dart_wrapper
  {
#ifdef CGAL_PWRLE_TURN_V3
    typedef std::size_t Dart_info;
#endif // CGAL_PWRLE_TURN_V3
    typedef CGAL::Cell_attribute_with_point<CMap> Vertex_attrib;
    typedef std::tuple<Vertex_attrib> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_traits
<3, CGAL::Exact_predicates_inexact_constructions_kernel> MyTraits;

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2, 3,
                                            MyTraits, MyItems> LCC_3_cmap;

#define NB_TESTS 3 // 0 ... 2
static int nbtests=0;
static int starting_seed;

using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
bool unit_test_canonize(Curves_on_surface_topology<LCC_3_cmap>& cst,
                        std::vector<Path_on_surface<LCC_3_cmap> >& paths,
                        const char* msg,
                        bool
#ifdef CGAL_USE_BASIC_VIEWER
                        draw
#endif
                        , int testtorun)
{
  CGAL_USE(msg);
  bool res=true;

  if (testtorun==-1 || nbtests==testtorun)
  {
#ifdef CGAL_TRACE_PATH_TESTS
    std::cout<<"[Test "<<nbtests<<"] "<<msg<<": "<<std::flush;
#endif

    for (unsigned int i=0; i<paths.size(); ++i)
    {
      if (i>0 && !cst.are_freely_homotopic(paths[0], paths[i]))
      {
        std::cout<<"[Test "<<nbtests<<"] ERROR: ";
        std::cout<<"paths["<<i<<"]"
                 <<" is not homotopic to paths[0] while it should be."
                 <<std::endl;
        res=false;
      }
    }

#ifdef CGAL_USE_BASIC_VIEWER
    if (draw)
    {
      std::string title="Test "+std::to_string(nbtests);
      CGAL::draw(paths[0].get_map(), paths, title.c_str());
    }
#endif // CGAL_USE_BASIC_VIEWER

#ifdef CGAL_TRACE_PATH_TESTS
    if (res) { std::cout<<std::endl; }
#endif
  }

  ++nbtests;
  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_double_torus_quad(bool draw, int testtorun)
{
  bool res=true;
  std::vector<Path_on_surface<LCC_3_cmap> > paths;

  // Test 0 (double torus, three G1 cycles)
  {
    LCC_3_cmap lcc;
    if (!CGAL::load_off(lcc, "./data/double-torus.off"))
    {
      std::cout<<"PROBLEM reading file ./data/double-torus.off"<<std::endl;
      return false;
    }
    Curves_on_surface_topology<LCC_3_cmap> cst(lcc);

    for (unsigned int i=0; i<3; ++i)
    {
      Path_on_surface<LCC_3_cmap> p(lcc);
      CGAL::generate_g1_double_torus(p, i);
      paths.push_back(p);
    }

    if (!unit_test_canonize(cst, paths,
                            "canonize paths on double torus gen1",
                            draw, testtorun))
    { res=false; }
  }

  // Test 1 (double torus, one random path, deformed randomly two times)
  {
    paths.clear();
    LCC_3_cmap lcc;
    if (!CGAL::load_off(lcc, "./data/double-torus.off")) // "./data/double-torus-smooth.off"))
    {
      std::cout<<"PROBLEM reading file ./data/double-torus.off"<<std::endl; // ./data/double-torus-smooth.off"<<std::endl;
      return false;
    }
    Curves_on_surface_topology<LCC_3_cmap> cst(lcc);

    if (testtorun==-1 || nbtests==testtorun)
    {
      CGAL::Random random(static_cast<unsigned int>(starting_seed+nbtests));
      Path_on_surface<LCC_3_cmap> p(lcc);

      internal::generate_random_closed_path
        (p, static_cast<std::size_t>(random.get_int(5, 20)), random); // random path, length between 5 and 20
      paths.push_back(p);

      p.update_path_randomly(random);
      paths.push_back(p);

      p.update_path_randomly(random);
      paths.push_back(p);
    }

    if (!unit_test_canonize(cst, paths,
                            "random canonize paths on double torus",
                            draw, testtorun))
    { res=false; }
  }

  // Test 2 (3torus-smooth, one random path, deformed randomly two times)
  {
    paths.clear();
    LCC_3_cmap lcc;
    if (!CGAL::load_off(lcc, "./data/3torus-smooth.off")) // 3torus.off
    {
      std::cout<<"PROBLEM reading file ./data/3torus-smooth.off"<<std::endl; // 3torus.off
      return false;
    }
    Curves_on_surface_topology<LCC_3_cmap> cst(lcc);

    if (testtorun==-1 || nbtests==testtorun)
    {
      CGAL::Random random(static_cast<unsigned int>(starting_seed+nbtests));
      Path_on_surface<LCC_3_cmap> p(lcc);

      internal::generate_random_closed_path
        (p, static_cast<std::size_t>(random.get_int(5, 200)), random); // random path, length between 5 and 200
      paths.push_back(p);

      p.update_path_randomly(random);
      paths.push_back(p);

      p.update_path_randomly(random);
      paths.push_back(p);

      if (!unit_test_canonize(cst, paths,
                              "random canonize paths on double torus",
                              draw, testtorun))
      { res=false; }
    }
  }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
[[noreturn]] void usage(int /*argc*/, char** argv)
{
  std::cout<<"usage: "<<argv[0]<<" [-draw] [-test N] [-seed S]"<<std::endl
           <<"   Test path homotopy method."
           <<std::endl
           <<"   -draw: draw mesh and paths."
           <<std::endl
           <<"   -test N: only run test number N (0<=N<"<<NB_TESTS<<")."
           <<std::endl
           <<"   -seed S: uses S as seed of random generator. Otherwise use a different seed at each run (based on time)."
           <<std::endl
           <<std::endl;
  exit(EXIT_FAILURE);
}
///////////////////////////////////////////////////////////////////////////////
[[noreturn]] void error_command_line(int argc, char** argv, const char* msg)
{
  std::cout<<"ERROR: "<<msg<<std::endl;
  usage(argc, argv);
}
///////////////////////////////////////////////////////////////////////////////
int main(int argc, char** argv)
{
  bool draw=false;
  std::string arg;
  int testN=-1;

  CGAL::Random r; // Used when user do not provide its own seed.
  starting_seed=r.get_int(0,INT_MAX);

  for (int i=1; i<argc; ++i)
  {
    arg=argv[i];
    if (arg=="-draw")
    { draw=true; }
    else if (arg=="-test")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "Error: no number after -test option."); }
      testN=std::stoi(std::string(argv[++i]));
      if (testN<0 || testN>=NB_TESTS)
      { error_command_line(argc, argv, "Error: invalid value for -test option."); }
    }
    else if (arg=="-seed")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "Error: no number after -seed option."); }
      starting_seed=std::stoi(std::string(argv[++i]));
    }
    else if (arg=="-h" || arg=="--help" || arg=="-?")
    { usage(argc, argv); }
    else if (arg[0]=='-')
    { std::cout<<"Unknown option "<<arg<<", ignored."<<std::endl; }
  }

  std::cout<<"Start tests (seed="<<starting_seed<<"): "<<std::flush;

#ifdef CGAL_TRACE_PATH_TESTS
  std::cout<<std::endl;
#endif

  if (!test_double_torus_quad(draw, testN))
  {
    std::cout<<"TEST DOUBLE TORUS FAILED."<<std::endl;
    return EXIT_FAILURE;
  }

  if (testN==-1)
  { std::cout<<"all the "<<nbtests<<" tests OK."<<std::endl; }
  else
  { std::cout<<"test "<<testN<<" OK."<<std::endl; }

  return EXIT_SUCCESS;
}
