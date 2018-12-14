#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Path_generators.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Path_on_surface_with_rle.h>
#include <vector>
#include <sstream>

#include "Creation_of_test_cases_for_paths.h"

/* If you want to use a viewer, you can use qglviewer. */
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/draw_lcc_with_paths.h>
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
    typedef CGAL::cpp11::tuple<Vertex_attrib> Attributes;
  };
};

typedef CGAL::Linear_cell_complex_traits
<3, CGAL::Exact_predicates_inexact_constructions_kernel> MyTraits;

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2, 3,
                                            MyTraits, MyItems> LCC_3_cmap;

#define NB_TESTS 23 // 0 ... 22
static int nbtests=0;
static int starting_seed;

enum Transformation // enum for the type of transformations
{
  REDUCTION,
  PUSH,
  FULL_SIMPLIFICATION
};

///////////////////////////////////////////////////////////////////////////////
void transform_path(CGAL::Path_on_surface<LCC_3_cmap>& path, Transformation t,
                    bool start_with_positive_turn,
                    bool draw=false,
                    std::size_t repeat=0) // If 0, repeat as long as there is one modifcation;
                                           // otherwise repeat the given number of times
{
  std::vector<CGAL::Path_on_surface<LCC_3_cmap> > v;
#ifdef CGAL_USE_BASIC_VIEWER
  if (draw)
  {
    v.push_back(path);
    // CGAL::draw(path.get_map(), v);
  }
#endif // CGAL_USE_BASIC_VIEWER

  CGAL::Path_on_surface<LCC_3_cmap> prevp=path;
  CGAL::Path_on_surface_with_rle<LCC_3_cmap> curp(path.get_map());
  std::size_t nb=0;
  bool modified=false;
  do
  {
    curp=CGAL::Path_on_surface_with_rle<LCC_3_cmap>(prevp,
                                                    start_with_positive_turn);
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
      prevp=CGAL::Path_on_surface<LCC_3_cmap>(curp);
#ifdef CGAL_USE_BASIC_VIEWER
      if (draw) { v.push_back(prevp); }
#endif // CGAL_USE_BASIC_VIEWER

      /* curp->display_negative_turns();
      std::cout<<"  "; curp->display_positive_turns();
      std::cout<<std::endl; */
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
}
///////////////////////////////////////////////////////////////////////////////
bool unit_test(CGAL::Path_on_surface<LCC_3_cmap>& path, Transformation t,
               std::size_t repeat,
               const char* msg, const char* expected_result,
               bool draw, int testtorun,
               bool start_with_positive_turn)
{
  bool res=true;

  if (testtorun==-1 || nbtests==testtorun)
  {
#ifdef CGAL_TRACE_PATH_TESTS
    std::cout<<"[Test "<<nbtests<<"] "<<msg<<": "<<std::flush;
    path.display_pos_and_neg_turns();
#else
    std::cout<<"."<<std::flush;
#endif

    transform_path(path, t, start_with_positive_turn, draw, repeat);

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
bool test_all_cases_spurs_and_bracket(bool draw, int testtorun)
{
  bool res=true;
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    return false;
  }

  CGAL::Path_on_surface<LCC_3_cmap> path(lcc);

  generate_one_positive_spur(path); // Test 0
  if (!unit_test(path, REDUCTION, 1, "Positive spur (2^6 1 0 2^4)",
                 "2 2 2 2 2 2 3 2 2 2", draw, testtorun, true))
  { res=false; }

  generate_one_negative_spur(path); // Test 1
  if (!unit_test(path, REDUCTION, 1, "Negative spur (-2^6 -1 0 -2^4)",
                 "-2 -2 -2 -2 -2 -2 -3 -2 -2 -2", draw, testtorun, false))
  { res=false; }

  generate_cyclic_spur(path); // Test  2
  if (!unit_test(path, REDUCTION, 1, "Cyclic spur (0 0)",
                 "", draw, testtorun, true))
  { res=false; }

  generate_one_positive_bracket(path); // Test 3
  if (!unit_test(path, REDUCTION, 1, "Positive bracket (2^3 3 1 2^6 1 3 2^2)",
                 "2 2 2 2 -2 -2 -2 -2 -2 -2 2 2 2", draw, testtorun, true))
  { res=false; }

  generate_one_negative_bracket(path); // Test 4
  if (!unit_test(path, REDUCTION, 1, "Negative bracket (-2^3 -1 -2^6 -1 -2^2)",
                 "-2 -2 3 2 2 2 2 2 2 3 -2 -2", draw, testtorun, false))
  { res=false; }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/spiral-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/spiral-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  generate_positive_bracket_special1(path); // Test 5
  if (!unit_test(path, REDUCTION, 1, "Positive special case 1 (4 1 2^7 1)",
                 "-2 -2 -2 -2 -2 -2 -2 2", draw, testtorun, true))
  { res=false; }

  generate_negative_bracket_special1(path); // Test 6
  if (!unit_test(path, REDUCTION, 1, "Negative special case 1 (-4 -1 -2^7 -1)",
                 "2 2 2 2 2 2 2 -2", draw, testtorun, false))
  { res=false; }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/loop-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/spiral-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  generate_positive_bracket_special2(path); // Test 7
  if (!unit_test(path, REDUCTION, 1, "Positive special case 2 (1 2^11)",
                 "-2 -2 -2 -2 -2 -2 -2 -2 -2 -3", draw, testtorun, true))
  { res=false; }

  generate_negative_bracket_special2(path); // Test 8
  if (!unit_test(path, REDUCTION, 1, "Negative special case 2 (-1 -2^11)",
                 "2 2 2 2 2 2 2 2 2 3", draw, testtorun, false))
  { res=false; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_all_cases_l_shape(bool draw, int testtorun)
{
  bool res=true;
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    return false;
  }
  CGAL::Path_on_surface<LCC_3_cmap> path(lcc);

  generate_one_l_shape(path); // Test 9
  if (!unit_test(path, PUSH, 1, "L-shape (-2^2 -3 -2^8 -1 -2^5 -3 -2^3)",
                 "-2 -2 2 1 2 2 2 2 2 2 2 3 2 2 2 2 1 2 -2 -2 -2",
                 draw, testtorun, false))
  { res=false; }

  generate_l_shape_case2(path); // Test 10
  if (!unit_test(path, PUSH, 1, "L-shape (-2^2 -3 -1 -2^5 -3 -2^3)",
                 "-2 -2 2 2 2 2 2 2 1 2 2 2 2", draw, testtorun, false))
  { res=false; }

  generate_l_shape_case3(path); // Test 11
  if (!unit_test(path, PUSH, 1, "L-shape (-2^2 -3 -2^5 -1 -3 -2^3)",
                 "-2 -2 2 1 2 2 2 2 2 2 -2 -2 -2", draw, testtorun, false))
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
                 "4 1 2 2 2 2 2 2 3 2 2 1", draw, testtorun, false))
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
                 "4 2 2 2 2 2 2 2 2 2 2 2 2 1", draw, testtorun, false))
  { res=false; }

  generate_l_shape_case6(path); // Test 14
  if (!unit_test(path, PUSH, 1, "L-shape (-4 -2^12 -1)",
                 "4 1 2 2 2 2 2 2 2 2 2 2 2 2", draw, testtorun, false))
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
                 "1 2 2 2 2 2 2 2 3 2 2 2", draw, testtorun, false))
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
                 draw, testtorun, false))
  { res=false; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_some_random_paths_on_cube(bool draw, int testtorun)
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    return false;
  }

  CGAL::Path_on_surface<LCC_3_cmap> path(lcc);
  bool res=true;

  CGAL::Random random(starting_seed+nbtests); // fix seed
  generate_random_bracket(path, 2, 6, 3, random); // Test 17
  if (!unit_test(path, FULL_SIMPLIFICATION, 0, "(2^1 1 2^6 1 2^3 ... )",
                 "2 2", draw, testtorun, true))
  { res=false; }

  random=CGAL::Random(starting_seed+nbtests);
  generate_random_bracket(path, 3, 8, 4, random); // Test 18
  if (!unit_test(path, FULL_SIMPLIFICATION, 0, "(2^2 1 2^8 1 2^4 ... )",
                   "", draw, testtorun, true))
  { res=false; }

  random=CGAL::Random(starting_seed+nbtests);
  generate_random_bracket(path, 5, 12, 8, random); // Test 19
  if (!unit_test(path, FULL_SIMPLIFICATION, 0, "(2^4 1 2^12 1 2^8 ...)",
                 "1", draw, testtorun, true))
  { res=false; }

  random=CGAL::Random(starting_seed+nbtests);
  generate_random_bracket(path, 5, 12, 8, random); // Test 20
  if (!unit_test(path, FULL_SIMPLIFICATION, 0, "(2^4 1 2^12 1 2^8 ...)",
                   "2", draw, testtorun, true))
  { res=false; }

  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_right_push(bool draw, int testtorun)
{
  bool res=true;

  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    return false;
  }

  CGAL::Path_on_surface<LCC_3_cmap> path(lcc);

  if (testtorun==-1 || nbtests==testtorun)
  {
    generate_one_l_shape(path); // Test 21
    transform_path(path, PUSH, false, draw, 0); // TODO COMPARE RESULT WITH GOLD STANDARD
  }
  ++nbtests;

  if (testtorun==-1 || nbtests==testtorun)
  {
    generate_l_shape_case2(path); // Test 22
    transform_path(path, PUSH, false, draw, 0); // TODO COMPARE RESULT WITH GOLD STANDARD
  }
  ++nbtests;

  return res;
}
///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
void usage(int /*argc*/, char** argv)
{
  std::cout<<"usage: "<<argv[0]<<" [-draw] [-test N] [-seed S]"<<std::endl
           <<"   Test several path transformations "
           <<"(bracket flattening, spurs removal and right shift of l-shape)."
           <<std::endl
           <<"   -draw: draw mesh and paths before and after the transformation"
           <<std::endl
           <<"   -test N: only run test number N (0<=N<"<<NB_TESTS<<")."
           <<std::endl
           <<"   -seed S: uses S as seed of random generator. Otherwise use a different seed at each run (based on time)."
           <<std::endl
           <<std::endl;
  exit(EXIT_FAILURE);
}
///////////////////////////////////////////////////////////////////////////////
void error_command_line(int argc, char** argv, const char* msg)
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

  std::cout<<"Start tests (seed="<<starting_seed<<"):"<<std::flush;

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

  if (!test_right_push(draw, testN))
  {
    std::cout<<"TEST RIGHT PUSH FAILED."<<std::endl;
    return EXIT_FAILURE;
  }

  if (testN==-1)
  { std::cout<<"all the "<<nbtests<<" tests OK."<<std::endl; }
  else
  { std::cout<<"test "<<testN<<" OK."<<std::endl; }

  return EXIT_SUCCESS;
}
