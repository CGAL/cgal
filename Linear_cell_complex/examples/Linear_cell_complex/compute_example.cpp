#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Combinatorial_map_functionalities.h>
#include <vector>
#include <sstream>

/* If you want to use a viewer, you can use qglviewer. */
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/LCC_with_paths.h>
#endif

#include <CGAL/Path_generators.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Creation_of_test_cases_for_paths.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
typedef CGAL::Linear_cell_complex_for_generalized_map<2,3> LCC_3_gmap;

#define NB_TESTS 25
int nbtests=0;

///////////////////////////////////////////////////////////////////////////////
void simplify_path(CGAL::Path_on_surface<LCC_3_cmap>& path,
                   bool draw=false,
                   unsigned int repeat=0) // If 0, repeat as long as there is one modifcation;
                                          // otherwise repeat the given number of times
{
  std::vector<const CGAL::Path_on_surface<LCC_3_cmap>*> v;
  if (draw)
  {
    v.push_back(&path);
    // display(path.get_map(), v);
  }

  CGAL::Path_on_surface<LCC_3_cmap>* prevp=&path;
  CGAL::Path_on_surface<LCC_3_cmap>* curp=NULL;
  unsigned int nb=0;
  do
  {
    curp=new CGAL::Path_on_surface<LCC_3_cmap>(*prevp);
    /* curp->display_negative_turns();
    std::cout<<"  "; curp->display_positive_turns();
    std::cout<<" -> "<<std::flush; */

    if (curp->bracket_flattening_one_step())
    {
      if (draw) { v.push_back(curp); }
      prevp=curp;

      /* curp->display_negative_turns();
      std::cout<<"  "; curp->display_positive_turns();
      std::cout<<std::endl; */
    }
    else
    {
      if (curp->remove_spurs())
      {
        if (draw) { v.push_back(curp); }
        prevp=curp;

       /* curp->display_negative_turns();
        std::cout<<"  "; curp->display_positive_turns();
        std::cout<<std::endl; */
      }
      else
      {
        delete curp;
        curp=NULL;
        // std::cout<<"unchanged."<<std::endl;
      }
    }
     // if (draw /* && nbtest==1*/)
     // display(path.get_map(), v);

    ++nb;
  }
  while((repeat==0 && curp!=NULL) || (nb<repeat));

  if (draw)
  {
    std::string title="Test "+std::to_string(nbtests-1);
    display(path.get_map(), v, title.c_str());
  }

  path.swap(*prevp);
}
///////////////////////////////////////////////////////////////////////////////
void push_l_shape(CGAL::Path_on_surface<LCC_3_cmap>& path,
                  bool draw=false,
                  unsigned int repeat=0) // If 0, repeat as long as there is one modifcation;
                                         // otherwise repeat the given number of times
{
  std::vector<const CGAL::Path_on_surface<LCC_3_cmap>*> v;
  if (draw)
  {
    v.push_back(&path);
    // display(path.get_map(), v);
  }

  CGAL::Path_on_surface<LCC_3_cmap>* prevp=&path;
  CGAL::Path_on_surface<LCC_3_cmap>* curp=NULL;
  unsigned int nb=0;
  do
  {
    curp=new CGAL::Path_on_surface<LCC_3_cmap>(*prevp);
    /* curp->display_negative_turns();
    std::cout<<"  "; curp->display_positive_turns();
    std::cout<<" -> "<<std::flush; */

    if (curp->right_push_one_step())
    {
      if (draw) { v.push_back(curp); }
      prevp=curp;

      /* curp->display_negative_turns();
      std::cout<<"  "; curp->display_positive_turns();
      std::cout<<std::endl; */
    }
    else
    {
      delete curp;
      curp=NULL;
      // std::cout<<"unchanged."<<std::endl;
    }
    // if (draw /* && nbtest==1*/)
    // display(path.get_map(), v);

    ++nb;
  }
  while((repeat==0 && curp!=NULL) || (nb<repeat));

  if (draw)
  {
    std::string title="Test "+std::to_string(nbtests-1);
    display(path.get_map(), v, title.c_str());
  }

  path.swap(*prevp);
}
///////////////////////////////////////////////////////////////////////////////
bool test_simplify_random_path(const LCC_3_cmap& lcc,
                               std::size_t nb1, std::size_t nb2, std::size_t nb3,
                               CGAL::Random& random,
                               bool draw=false)
{
  CGAL::Path_on_surface<LCC_3_cmap> path(lcc);
  CGAL::initialize_path_random_starting_dart(path, random);
  CGAL::extend_straight_positive(path, nb1-1);
  CGAL::create_braket_positive(path, nb2);
  CGAL::extend_straight_positive(path, nb3);
  CGAL::generate_random_path(path, random.get_int(0, 15), random);

  simplify_path(path, draw);

  return path.is_valid();
}
///////////////////////////////////////////////////////////////////////////////
bool test_all_cases_spurs_and_bracket(bool draw, int testtorun)
{
  bool res=true;
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  /* std::cout<<"Initial map 1: ";
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid() << std::endl;
  */

  CGAL::Path_on_surface<LCC_3_cmap> path(lcc);

  if (nbtests++==testtorun || testtorun==-1) // Test 0
  {
    std::cout<<"."<<std::flush;
    generate_one_positive_spur(path);
    // std::cout<<"Case 1: Positive spur (2^6 1 0 2^4): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns("2 2 2 2 2 2 3 2 2 2"))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 1] ERROR: ";
      std::cout<<"we obtained "; path.display_positive_turns();
      std::cout<<" instead of +(2 2 2 2 2 2 3 2 2 2)"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 1
  {
    std::cout<<"."<<std::flush;
    generate_one_negative_spur(path);
    // std::cout<<"Case 2: Negative spur (-2^6 -1 0 -2^4): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns("-2 -2 -2 -2 -2 -2 -3 -2 -2 -2"))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 2] ERROR: ";
      std::cout<<"we obtained "; path.display_negative_turns();
      std::cout<<" instead of -(2 2 2 2 2 2 3 2 2 2)"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 2
  {
    std::cout<<"."<<std::flush;
    generate_cyclic_spur(path);
    // std::cout<<"Case 3: Cyclic spur (0 0): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns(""))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 3] ERROR: ";
      std::cout<<"we obtained "; path.display_positive_turns();
      std::cout<<" instead of +()"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 3
  {
    std::cout<<"."<<std::flush;
    generate_one_positive_bracket(path);
    // std::cout<<"Case 4: Positive bracket (2^3 3 1 2^6 1 3 2^2): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns("2 2 2 2 -2 -2 -2 -2 -2 -2 2 2 2"))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 4] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (2 2 2 2 -2 -2 -2 -2 -2 -2 2 2 2)"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 4
  {
    std::cout<<"."<<std::flush;
    generate_one_negative_bracket(path);
    // std::cout<<"Case 5: Negative bracket (-2^3 -1 -2^6 -1 -2^2): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns("-2 -2 3 2 2 2 2 2 2 3 -2 -2"))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 5] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (-2 -2 3 2 2 2 2 2 2 3 -2 -2)"<<std::endl;
      res=false;
    }
  }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/spiral-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/spiral-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 5
  {
    std::cout<<"."<<std::flush;
    generate_positive_bracket_special1(path);
    // std::cout<<"Case 6: Positive special case 1 (3 1 2^8 1): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns("-2 -2 -2 -2 -2 -2 -2 -2 1"))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 6] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (-2 -2 -2 -2 -2 -2 -2 -2 1)"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 6
  {
    std::cout<<"."<<std::flush;
    generate_negative_bracket_special1(path);
    // std::cout<<"Case 7: Negative special case 1 (-3 -1 -2^8 -1): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns("2 2 2 2 2 2 2 2 6"))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 7] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (2 2 2 2 2 2 2 2 6)"<<std::endl;
      res=false;
    }
  }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/loop-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/spiral-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 7
  {
    std::cout<<"."<<std::flush;
    generate_positive_bracket_special2(path);
    // std::cout<<"Case 8: Positive special case 1 (1 2^11): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns("-2 -2 -2 -2 -2 -2 -2 -2 -2 -3"))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 8] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (-2 -2 -2 -2 -2 -2 -2 -2 -2 -3)"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 8
  {
    std::cout<<"."<<std::flush;
    generate_negative_bracket_special2(path);
    // std::cout<<"Case 9: Negative special case 1 (-1 -2^11): "<<std::flush;
    simplify_path(path, draw, 1);
    if (!path.same_turns("2 2 2 2 2 2 2 2 2 3"))
    {
      std::cout<<"[test_all_cases_spurs_and_bracket case 9] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (2 2 2 2 2 2 2 2 2 3)"<<std::endl;
      res=false;
    }
  }

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
    exit(EXIT_FAILURE);
  }
  CGAL::Path_on_surface<LCC_3_cmap> path(lcc);

  if (nbtests++==testtorun || testtorun==-1) // Test 9
  {
    std::cout<<"."<<std::flush;
    generate_one_l_shape(path);
    // std::cout<<"Case 1: L-shape (-2^2 -3 -2^8 -1 -2^5 -3 -2^3): "<<std::flush;
    push_l_shape(path, draw, 1);
    if (!path.same_turns("-2 -2 2 1 2 2 2 2 2 2 2 3 2 2 2 2 1 2 -2 -2 -2"))
    {
      std::cout<<"[test_all_cases_l_shape case 1] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (-2 -2 2 1 2 2 2 2 2 2 2 3 2 2 2 2 1 2 -2 -2 -2)"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 10
  {
    std::cout<<"."<<std::flush;
    generate_l_shape_case2(path);
    // std::cout<<"Case 2: L-shape (-2^2 -3 -1 -2^5 -3 -2^3): "<<std::flush;
    push_l_shape(path, draw, 1);
    if (!path.same_turns("-2 -2 2 2 2 2 2 2 1 2 2 2 2"))
    {
      std::cout<<"[test_all_cases_l_shape case 2] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (-2 -2 2 2 2 2 2 2 1 2 2 2 2)"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 11
  {
    std::cout<<"."<<std::flush;
    generate_l_shape_case3(path);
    // std::cout<<"Case 3: L-shape (-2^2 -3 -2^5 -1 -3 -2^3): "<<std::flush;
    push_l_shape(path, draw, 1);
    if (!path.same_turns("-2 -2 2 1 2 2 2 2 2 2 -2 -2 -2"))
    {
      std::cout<<"[test_all_cases_l_shape case 3] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (-2 -2 2 1 2 2 2 2 2 2 -2 -2 -2)"<<std::endl;
      res=false;
    }
  }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/case4-right-shift-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/case4-right-shift-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  lcc.reverse_orientation();

  if (nbtests++==testtorun || testtorun==-1) // Test 12
  {
    std::cout<<"."<<std::flush;
    generate_l_shape_case4(path);
    // std::cout<<"Case 4: L-shape (-4 -2^7 -1 -2^3): "<<std::flush;
    push_l_shape(path, draw, 1);
    if (!path.same_turns("4 1 2 2 2 2 2 2 3 2 2 1"))
    {
      std::cout<<"[test_all_cases_l_shape case 4] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (4 1 2 2 2 2 2 2 3 2 2 1)"<<std::endl;
      res=false;
    }
  }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/cases5-6-right-shift-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cases5-6-right-shift-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  lcc.reverse_orientation();

  if (nbtests++==testtorun || testtorun==-1) // Test 13
  {
    std::cout<<"."<<std::flush;
    generate_l_shape_case5(path);
    // std::cout<<"Case 5: L-shape (-4 -1 -2^12): "<<std::flush;
    push_l_shape(path, draw, 1);
    if (!path.same_turns("4 2 2 2 2 2 2 2 2 2 2 2 2 1"))
    {
      std::cout<<"[test_all_cases_l_shape case 5] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (4 2 2 2 2 2 2 2 2 2 2 2 2 1)"<<std::endl;
      res=false;
    }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 14
  {
    std::cout<<"."<<std::flush;
    generate_l_shape_case6(path);
    // std::cout<<"Case 6: L-shape (-4 -2^12 -1): "<<std::flush;
    push_l_shape(path, draw, 1);
    if (!path.same_turns("4 1 2 2 2 2 2 2 2 2 2 2 2 2"))
    {
      std::cout<<"[test_all_cases_l_shape case 6] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (4 1 2 2 2 2 2 2 2 2 2 2 2 2)"<<std::endl;
      res=false;
    }
  }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/case7-right-shift-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/case7-right-shift-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 15
  {
    std::cout<<"."<<std::flush;
    generate_l_shape_case7(path);
    // std::cout<<"Case 7: L-shape (-3 -2^7 -1 -2^3): "<<std::flush;
    push_l_shape(path, draw, 1);
    if (!path.same_turns("1 2 2 2 2 2 2 2 3 2 2 2"))
    {
      std::cout<<"[test_all_cases_l_shape case 7] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (1 2 2 2 2 2 2 2 3 2 2 2)"<<std::endl;
      res=false;
    }
  }

  lcc.clear();
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 16
  {
    std::cout<<"."<<std::flush;
    generate_l_shape_case8(path);
    // std::cout<<"Case 8: L-shape (-2^20): "<<std::flush;
    push_l_shape(path, draw, 1);
    if (!path.same_turns("2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2"))
    {
      std::cout<<"[test_all_cases_l_shape case 8] ERROR: ";
      std::cout<<"we obtained "; path.display_pos_and_neg_turns();
      std::cout<<" instead of (2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2)"<<std::endl;
      res=false;
    }
  }

  // path.display_pos_and_neg_turns();std::cout<<std::endl;
  return res;
}
///////////////////////////////////////////////////////////////////////////////
bool test_some_random_paths_on_cube(bool draw, int testtorun)
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  CGAL::Random random(1); // fix seed

  if (nbtests++==testtorun || testtorun==-1) // Test 17
  {
    std::cout<<"."<<std::flush;
    // path 1: 2 straight darts to begin; 1 + turn; 6 straight; 1 + turn; 3 straight
    if (!test_simplify_random_path(lcc, 2, 6, 3, random, draw))
    { return false; }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 18
  {
    std::cout<<"."<<std::flush;
    // path 2: 3 straight darts to begin; 1 + turn; 8 straight; 1 + turn; 4 straight
    if (!test_simplify_random_path(lcc, 3, 8, 4, random, draw))
    { return false; }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 19
  {
    std::cout<<"."<<std::flush;
    // path 3: 5 straight darts to begin; 1 + turn; 12 straight; 1 + turn; 8 straight
    if (!test_simplify_random_path(lcc, 5, 12, 8, random, draw))
    { return false; }
  }

  if (nbtests++==testtorun || testtorun==-1) // Test 20
  {
    std::cout<<"."<<std::flush;
    // path 4: 5 straight darts to begin; 1 + turn; 12 straight; 1 + turn; 8 straight
    if (!test_simplify_random_path(lcc, 5, 12, 8, random, draw))
    { return false; }
  }

  return true;
}
///////////////////////////////////////////////////////////////////////////////
bool test_file(const std::string& filename, bool draw, int testtorun)
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, filename.c_str()))
  {
    std::cout<<"PROBLEM reading file "<<filename<<std::endl;
    exit(EXIT_FAILURE);
  }

 /* std::cout<<"Initial map: ";
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid() << std::endl; */

  CGAL::Random random(1); // fix seed
  CGAL::Path_on_surface<LCC_3_cmap> p1(lcc);
  generate_random_path(p1, 10, random); // 10 %
  CGAL::Path_on_surface<LCC_3_cmap> p2(lcc);
  generate_random_path(p2, 15, random);
  CGAL::Path_on_surface<LCC_3_cmap> p3(lcc);
  generate_random_path(p3, 10, random); // 10 %
  CGAL::Path_on_surface<LCC_3_cmap> p4(lcc);
  generate_random_path(p4, 15, random);

  CGAL::Combinatorial_map_tools<LCC_3_cmap> cmt(lcc);

  CGAL::Path_on_surface<LCC_3_cmap>
      pp1=cmt.transform_original_path_into_quad_surface(p1),
      pp2=cmt.transform_original_path_into_quad_surface(p2),
      pp3=cmt.transform_original_path_into_quad_surface(p3),
      pp4=cmt.transform_original_path_into_quad_surface(p4);

  if (nbtests++==testtorun || testtorun==-1) // Test 21
  { std::cout<<"."<<std::flush; simplify_path(pp1, draw); }
  if (nbtests++==testtorun || testtorun==-1) // Test 22
  { std::cout<<"."<<std::flush; simplify_path(pp2, draw); }
  if (nbtests++==testtorun || testtorun==-1) // Test 23
  { std::cout<<"."<<std::flush; simplify_path(pp3, draw); }
  if (nbtests++==testtorun || testtorun==-1) // Test 24
  { std::cout<<"."<<std::flush; simplify_path(pp4, draw); }

  return pp1.is_valid() && pp2.is_valid() && pp3.is_valid() && pp4.is_valid();
}
///////////////////////////////////////////////////////////////////////////////
void usage(int /*argc*/, char** argv)
{
  std::cout<<"usage: "<<argv[0]<<" [-draw] [filename.off] [-case N]"<<std::endl
           <<"   Test several path transformations "
           <<"(bracket flattening, spurs removal and right shift of l-shape)."
           <<std::endl
           <<"   -draw: draw mesh and paths before and after the transformation"
           <<std::endl
           <<"   filename.off: give the off file to use for the last test "
           <<"(random path on a given surface). Use data/elephant.off by default."
           <<std::endl
           <<"   -test N: only run test number N (0<=N<"<<NB_TESTS<<")."
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
  std::string filename="./data/elephant.off";
  std::string arg;
  int testN=-1;

  for (int i=1; i<argc; ++i)
  {
    arg=argv[i];
    if (arg=="-draw")
    { draw=true; }
    if (arg=="-test")
    {
      if (i==argc-1)
      { error_command_line(argc, argv, "Error: no number after -test option."); }
      testN=std::stoi(std::string(argv[++i]));
      if (testN<0 || testN>=NB_TESTS)
      { error_command_line(argc, argv, "Error: invalid value for -test option."); }
    }
    else if (arg=="-h" || arg=="--help" || arg=="-?")
    { usage(argc, argv); }
    else if (arg[0]=='-')
    { std::cout<<"Unknown option "<<arg<<", ignored."<<std::endl; }
    else { filename=argv[i]; }
  }

  std::cout<<"Start tests: "<<std::flush;

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

  if (!test_file(filename, draw, testN))
  {
    std::cout<<"TEST FILE FAILED."<<std::endl;
    return EXIT_FAILURE;
  }

  std::cout<<": all the "<<nbtests<<" tests OK."<<std::endl;


  return EXIT_SUCCESS;
}
