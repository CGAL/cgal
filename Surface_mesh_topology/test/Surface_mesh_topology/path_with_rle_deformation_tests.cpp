#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Surface_mesh_topology/internal/Path_generators.h>
#include <CGAL/Surface_mesh_topology/internal/Path_on_surface_with_rle.h>
#include <vector>
#include <sstream>
#include <tuple>

#include "Creation_of_test_cases_for_paths.h"

// If you want to use a viewer, you can use qglviewer.
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

#define NB_TESTS 21 // 0 ... 20
static unsigned int nbtests=0;

static const unsigned int ALL_TESTS=(std::numeric_limits<unsigned int>::max)();

enum Transformation // enum for the type of transformations
{
  REDUCTION,
  PUSH,
  FULL_SIMPLIFICATION
};

using namespace CGAL::Surface_mesh_topology;

///////////////////////////////////////////////////////////////////////////////
void transform_path(Path_on_surface<LCC_3_cmap>& path, Transformation t,
                    bool use_only_positive,
                    bool use_only_negative,
                    bool
#ifdef CGAL_USE_BASIC_VIEWER
                    draw
#endif
                    =false,
                    std::size_t repeat=0) // If 0, repeat as long as there is one modifcation;
                                           // otherwise repeat the given number of times
{
#ifdef CGAL_USE_BASIC_VIEWER
  std::vector<Path_on_surface<LCC_3_cmap> > v;
  if (draw)
  {
    v.push_back(path);
    // CGAL::draw(path.get_map(), v);
  }
#endif // CGAL_USE_BASIC_VIEWER

  internal::Light_MQ<LCC_3_cmap> lmq(path.get_map());
  Path_on_surface<LCC_3_cmap> prepp(path.get_map());
  Path_on_surface<LCC_3_cmap> prevp=path;
  internal::Path_on_surface_with_rle<internal::Light_MQ<LCC_3_cmap> >curp(lmq);
  std::size_t nb=0;
  bool modified=false;
  do
  {
    curp=internal::Path_on_surface_with_rle<internal::Light_MQ<LCC_3_cmap> >
        (lmq, prevp, use_only_positive, use_only_negative);

    modified=false;
    /* curp->display_negative_turns();
    std::cout<<"  "; curp->display_positive_turns();
    std::cout<<" -> "<<std::flush; */

    if (t==REDUCTION || t==FULL_SIMPLIFICATION)
    {
      modified=curp.remove_brackets(false);
      if (!modified)
      { modified=curp.remove_spurs(false); }
    }
    if (t==PUSH || t==FULL_SIMPLIFICATION)
    {
      if (!modified)
      { modified=curp.right_push(false); }
    }

    if (modified)
    {
      prevp=Path_on_surface<LCC_3_cmap>(curp);
#ifdef CGAL_USE_BASIC_VIEWER
      if (draw) { v.push_back(prevp); }
#endif // CGAL_USE_BASIC_VIEWER
    }

    // if (draw /* && nbtest==1*/)
    // CGAL::draw(path.get_map(), v);

    ++nb;
  }
  while((repeat==0 && modified) || (nb<repeat));

#ifdef CGAL_USE_BASIC_VIEWER
  if (draw)
  {
    std::string title="Test "+std::to_string(nbtests);
    CGAL::draw(path.get_map(), v, title.c_str());
  }
#endif // CGAL_USE_BASIC_VIEWER

  path.swap(prevp);
  CGAL_assertion(path.is_valid(true));
}
///////////////////////////////////////////////////////////////////////////////
bool unit_test(Path_on_surface<LCC_3_cmap>& path, Transformation t,
               std::size_t repeat,
               const char* msg, const char* expected_result,
               bool draw, unsigned int testtorun,
               bool use_only_positive,
               bool use_only_negative)
{
  CGAL_USE(msg);
  bool res=true;

  if (testtorun==ALL_TESTS || nbtests==testtorun)
  {
#ifdef CGAL_TRACE_PATH_TESTS
    std::cout<<"[Test "<<nbtests<<"] "<<msg<<": "<<std::flush;
    path.display_pos_and_neg_turns();
#else
    std::cout<<"."<<std::flush;
#endif

    transform_path(path, t, use_only_positive, use_only_negative, draw, repeat);

    if (!path.same_turns(expected_result))
    {
      std::cout<<"[Test "<<nbtests<<"] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of ("<<expected_result<<")"<<std::endl;
      res=false;
    }

#ifdef CGAL_TRACE_PATH_TESTS
    std::cout<<" -> "<<std::flush; path.display_pos_and_neg_turns();
    std::cout<<std::endl;
#endif
  }

  ++nbtests;
  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_all_cases_spurs_and_bracket(bool draw, unsigned int testtorun)
{
  bool res=true;
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    return false;
  }

  Path_on_surface<LCC_3_cmap> path(lcc);

  generate_one_positive_spur(path); // Test 0
  if (!unit_test(path, REDUCTION, 1, "Positive spur (2^6 1 0 2^4)",
                 "2 2 2 2 2 2 3 2 2 2", draw, testtorun, true, false))
  { res=false; }

  generate_one_negative_spur(path); // Test 1
  if (!unit_test(path, REDUCTION, 1, "Negative spur (-2^6 -1 0 -2^4)",
                 "-2 -2 -2 -2 -2 -2 -3 -2 -2 -2", draw, testtorun, false, true))
  { res=false; }

  generate_cyclic_spur(path); // Test  2
  if (!unit_test(path, REDUCTION, 1, "Cyclic spur (0 0)",
                 "", draw, testtorun, false, false))
  { res=false; }

  generate_one_positive_bracket(path); // Test 3
  if (!unit_test(path, REDUCTION, 1, "Positive bracket (2^3 3 1 2^6 1 3 2^2)",
                 "2 2 2 2 -2 -2 -2 -2 -2 -2 2 2 2", draw, testtorun, true, false))
  { res=false; }

  generate_one_negative_bracket(path); // Test 4
  if (!unit_test(path, REDUCTION, 1, "Negative bracket (-2^3 -1 -2^6 -1 -2^2)",
                 "-2 -2 3 2 2 2 2 2 2 3 -2 -2", draw, testtorun, false, true))
  { res=false; }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/spiral-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/spiral-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  generate_positive_bracket_special1(path); // Test 5
  if (!unit_test(path, REDUCTION, 1, "Positive special case 1 (4 1 2^7 1)",
                 "-2 -2 -2 -2 -2 -2 -2 2", draw, testtorun, true, false))
  { res=false; }

  generate_negative_bracket_special1(path); // Test 6
  if (!unit_test(path, REDUCTION, 1, "Negative special case 1 (-4 -1 -2^7 -1)",
                 "2 2 2 2 2 2 2 -2", draw, testtorun, false, true))
  { res=false; }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/loop-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/spiral-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  generate_positive_bracket_special2(path); // Test 7
  if (!unit_test(path, REDUCTION, 1, "Positive special case 2 (1 2^11)",
                 "-2 -2 -2 -2 -2 -2 -2 -2 -2 -3", draw, testtorun, true, false))
  { res=false; }

  generate_negative_bracket_special2(path); // Test 8
  if (!unit_test(path, REDUCTION, 1, "Negative special case 2 (-1 -2^11)",
                 "2 2 2 2 2 2 2 2 2 3", draw, testtorun, false, true))
  { res=false; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_all_cases_l_shape(bool draw, unsigned int testtorun)
{
  bool res=true;
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    return false;
  }
  Path_on_surface<LCC_3_cmap> path(lcc);

  generate_one_l_shape(path); // Test 9
  if (!unit_test(path, PUSH, 1, "L-shape (-2^2 -3 -2^8 -1 -2^5 -3 -2^3)",
                 "-2 -2 2 1 2 2 2 2 2 2 2 3 2 2 2 2 1 2 -2 -2 -2",
                 draw, testtorun, false, true))
  { res=false; }

  generate_l_shape_case2(path); // Test 10
  if (!unit_test(path, PUSH, 1, "L-shape (-2^2 -3 -1 -2^5 -3 -2^3)",
                 "-2 -2 2 2 2 2 2 2 1 2 2 2 2", draw, testtorun, false, true))
  { res=false; }

  generate_l_shape_case3(path); // Test 11
  if (!unit_test(path, PUSH, 1, "L-shape (-2^2 -3 -2^5 -1 -3 -2^3)",
                 "-2 -2 2 1 2 2 2 2 2 2 -2 -2 -2", draw, testtorun, false, true))
  { res=false; }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/case4-right-shift-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/case4-right-shift-squared.off"
             <<std::endl;
    exit(EXIT_FAILURE);
  }
  lcc.reverse_orientation();

  generate_l_shape_case4(path); // Test 12
  if (!unit_test(path, PUSH, 1, "L-shape (-4 -2^7 -1 -2^3)",
                 "4 1 2 2 2 2 2 2 3 2 2 1", draw, testtorun, false, true))
  { res=false; }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/cases5-6-right-shift-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cases5-6-right-shift-squared.off"
             <<std::endl;
    exit(EXIT_FAILURE);
  }
  lcc.reverse_orientation();

  generate_l_shape_case5(path); // Test 13
  if (!unit_test(path, PUSH, 1, "L-shape (-4 -1 -2^12)",
                 "4 2 2 2 2 2 2 2 2 2 2 2 2 1", draw, testtorun, false, true))
  { res=false; }

  generate_l_shape_case6(path); // Test 14
  if (!unit_test(path, PUSH, 1, "L-shape (-4 -2^12 -1)",
                 "4 1 2 2 2 2 2 2 2 2 2 2 2 2", draw, testtorun, false, true))
  { res=false; }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/case7-right-shift-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/case7-right-shift-squared.off"
             <<std::endl;
    exit(EXIT_FAILURE);
  }

  generate_l_shape_case7(path); // Test 15
  if (!unit_test(path, PUSH, 1, "L-shape (-3 -2^7 -1 -2^3), false",
                 "1 2 2 2 2 2 2 2 3 2 2 2", draw, testtorun, false, true))
  { res=false; }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  generate_l_shape_case8(path); // Test 16
  if (!unit_test(path, PUSH, 1, "L-shape (-2^20)",
                 "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2",
                 draw, testtorun, false, true))
  { res=false; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_some_random_paths_on_cube(bool draw, unsigned int testtorun)
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    return false;
  }

  Path_on_surface<LCC_3_cmap> path(lcc);
  bool res=true;

  generate_path(path, // Test 17
                {-8,8,6,0,-10,10,-7,0,2,2,-7,10,-2,3,-9,2,2,8,-6,5,0,
                    8,3,2,2,10,-10,7,9,4,-10,0,10,6,-5,8,-5,6,3,6}, // pair of numbers: (turn, length of flat)
                {3,0,-6,4,-8,6,10,1,-6,7,-5,0,5,1,2,5,-10,10,10,6,-7,6
                    ,1,10,8,5,7,4,2,7,10,7,10,0,-8,6,7,0,1,0,-2,4,-4,9,4,9,0,18});
  if (!unit_test(path, FULL_SIMPLIFICATION, 0, "first random path on cube",
                 "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2", draw, testtorun, true, false))
  { res=false; }

  generate_path(path, {1,3,0,6,9,4,1,10,-3,3}, // Test 18
                {-3,4,-10,4,-3,0,-6,10,-10,4});
  if (!unit_test(path, FULL_SIMPLIFICATION, 0, "second random path on cube",
                   "1 2 2 2 2", draw, testtorun, true, false))
  { res=false; }

  generate_path(path, {10,6,-4,5,-5,8,-2,2,4,9,5,3,8,4,-7,7}, // Test 19
                {6,0,-7,7,8,9,3,0,-9,1,-8,2,-3,1,4,8});
  if (!unit_test(path, FULL_SIMPLIFICATION, 0, "third random path on cube",
                 "2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 1 2 2 2 2 2 2",
                 draw, testtorun, true, false))
  { res=false; }

  generate_path(path, // Test 20
                {-2,10,6,5,8,1,-2,10,3,8,-7,7,5,10,1,5,9,0,-5,2,-2,4,7,8,-9,9,
                    -9,0,10,3,1,10,-4,8,5,4,-3,5,7,5,-4,7,-4,10,-8,8,2,0,2,1},
                {6,0,0,9,-3,7,-10,10,6,7,1,9,10,7,2,6,1,0,-6,8,-8,8,9,9,-8,2
                    ,-2,2,4,8,-2,6,3,10,-2,8,1,6,-7,5,7,6,5,5,6,7,-2,6,1,4});
  if (!unit_test(path, FULL_SIMPLIFICATION, 0, "fourth random path on cube",
                   "1 2 2 2 2 2 2", draw, testtorun, true, false))
  { res=false; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
[[noreturn]] void usage(int /*argc*/, char** argv)
{
  std::cout<<"usage: "<<argv[0]<<" [-draw] [-test N]"<<std::endl
           <<"   Test several path transformations "
           <<"(bracket flattening, spurs removal and right shift of l-shape)."
           <<std::endl
           <<"   -draw: draw mesh and paths before and after the transformation"
           <<std::endl
           <<"   -test N: only run test number N (0<=N<"<<NB_TESTS<<")."
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
  unsigned int testN=ALL_TESTS;

  for (int i=1; i<argc; ++i)
  {
    arg=argv[i];
    if (arg=="-draw")
    { draw=true; }
    else if (arg=="-test")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "Error: no number after -test option."); }
      testN=static_cast<unsigned int>(std::stoi(std::string(argv[++i])));
      if (testN>=NB_TESTS)
      { error_command_line(argc, argv, "Error: invalid value for -test option."); }
    }
    else if (arg=="-h" || arg=="--help" || arg=="-?")
    { usage(argc, argv); }
    else if (arg[0]=='-')
    { std::cout<<"Unknown option "<<arg<<", ignored."<<std::endl; }
  }

#ifdef CGAL_TRACE_PATH_TESTS
  std::cout<<std::endl;
#endif

  if (!test_all_cases_spurs_and_bracket(draw, testN))
  {
    std::cout<<"TEST SPURS AND BRACKET FAILED."<<std::endl;
    return EXIT_FAILURE;
  }

  if (!test_all_cases_l_shape(draw, testN))
  {
    std::cout<<"TEST L_SHAPE FAILED."<<std::endl;
    return EXIT_FAILURE;
  }

  if (!test_some_random_paths_on_cube(draw, testN))
  {
    std::cout<<"TEST RANDOM PATHS ON CUBE FAILED."<<std::endl;
    return EXIT_FAILURE;
  }

  if (testN==ALL_TESTS)
  { std::cout<<"all the "<<nbtests<<" tests OK."<<std::endl; }
  else
  { std::cout<<"test "<<testN<<" OK."<<std::endl; }

  return EXIT_SUCCESS;
}
