#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Combinatorial_map_functionalities.h>
#include <vector>

/* If you want to use a viewer, you can use qglviewer. */
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/LCC_with_paths.h>
#endif

#include <CGAL/Path_generators.h>
#include <CGAL/Path_on_surface.h>
#include <CGAL/Creation_of_test_cases_for_paths.h>

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
typedef CGAL::Linear_cell_complex_for_generalized_map<2,3> LCC_3_gmap;

void simplify_path(CGAL::Path_on_surface<LCC_3_cmap>& path,
                   bool draw=false)
{
  std::vector<const CGAL::Path_on_surface<LCC_3_cmap>*> v;
  if (draw)
  {
    v.push_back(&path);
    display(path.get_map(), v);
  }

  CGAL::Path_on_surface<LCC_3_cmap>* prevp=&path;
  CGAL::Path_on_surface<LCC_3_cmap>* curp=NULL;
  do
  {
    curp=new CGAL::Path_on_surface<LCC_3_cmap>(*prevp);
    if (curp->bracket_flattening_one_step())
    {
      if (draw) { v.push_back(curp); }
      prevp=curp;
    }
    else
    {
      delete curp;
      curp=NULL;
    }
    // if (nbtest==1)
     // display(lcc, v);
  }
  while(curp!=NULL);

  curp=new CGAL::Path_on_surface<LCC_3_cmap>(*prevp);
  if (curp->remove_spurs())
  {
    if (draw) { v.push_back(curp); }
    prevp=curp;
  }
  else
  {
    delete curp;
    curp=NULL;
  }

  path.swap(*prevp);
  if (draw)
  { display(path.get_map(), v); }
}

void test_file(int argc, char** argv)
{
  if (argc!=2)
  {
    std::cout<<"Usage: "<<argv[0]<<" filename.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, argv[1]))
  {
    std::cout<<"PROBLEM reading file "<<argv[1]<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout<<"Initial map: ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  CGAL::Random random;
  CGAL::Path_on_surface<LCC_3_cmap> p1(lcc);
  generate_random_path(p1, 10, random); // 10 %
  CGAL::Path_on_surface<LCC_3_cmap> p2(lcc);
  generate_random_path(p2, 15, random);
  CGAL::Path_on_surface<LCC_3_cmap> p3(lcc);
  generate_random_path(p3, 10, random); // 10 %
  CGAL::Path_on_surface<LCC_3_cmap> p4(lcc);
  generate_random_path(p4, 15, random);

  std::vector<const CGAL::Path_on_surface<LCC_3_cmap>*> v;
  v.push_back(&p1);
  v.push_back(&p2);
  v.push_back(&p3);
  v.push_back(&p4);  

#ifdef CGAL_USE_BASIC_VIEWER
  display(lcc, v);
#endif // CGAL_USE_BASIC_VIEWER
    
  CGAL::Combinatorial_map_tools<LCC_3_cmap> cmt(lcc);

  CGAL::Path_on_surface<LCC_3_cmap>
      pp1=cmt.transform_original_path_into_quad_surface(p1);

  std::cout<<"Original path has "<<pp1.length()<<" darts."<<std::endl;
  simplify_path(pp1, false);
  std::cout<<"After bracket flattening, the path has "<<pp1.length()<<" darts."<<std::endl;
}

void test_simplify_random_path(const LCC_3_cmap& lcc,
                               std::size_t nb1, std::size_t nb2, std::size_t nb3,
                               CGAL::Random& random,
                               bool draw=false)
{
  static std::size_t nbtest=0;
  std::cout<<"[BEGIN] TEST "<<nbtest<<"."<<std::endl;

  CGAL::Path_on_surface<LCC_3_cmap> path(lcc);
  CGAL::initialize_path_random_starting_dart(path, random);
  CGAL::extend_straight_positive(path, nb1-1);
  CGAL::create_braket_positive(path, nb2);
  CGAL::extend_straight_positive(path, nb3);
  CGAL::generate_random_path(path, random.get_int(0, 15), random);

  simplify_path(path, draw);

  std::cout<<"[END] TEST "<<nbtest++<<"."<<std::endl;
}

void test_some_random_paths_on_cube()
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout<<"Initial map: ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  CGAL::Random random(1); // fix seed

  // path 1: 2 straight darts to begin; 1 + turn; 6 straight; 1 + turn; 3 straight
  test_simplify_random_path(lcc, 2, 6, 3, random, true);

  // path 2: 3 straight darts to begin; 1 + turn; 8 straight; 1 + turn; 4 straight
  test_simplify_random_path(lcc, 3, 8, 4, random, true);

  // path 3: 5 straight darts to begin; 1 + turn; 12 straight; 1 + turn; 8 straight
  test_simplify_random_path(lcc, 5, 12, 8, random, true);

  // path 4: 5 straight darts to begin; 1 + turn; 12 straight; 1 + turn; 8 straight
  test_simplify_random_path(lcc, 5, 12, 8, random, true);
}

void test_all_cases_spurs_and_bracket()
{
  LCC_3_cmap lcc1;
  if (!CGAL::load_off(lcc1, "./data/cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./data/cube-mesh-5-5.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  LCC_3_cmap lcc2;
  if (!CGAL::load_off(lcc2, "./data/spiral-squared.off"))
  {
    std::cout<<"PROBLEM reading file ./data/spiral-squared.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  // lcc2.reverse_orientation();
  /* std::cout<<"Initial map 1: ";
  lcc.display_characteristics(std::cout) << ", valid="
                                         << lcc.is_valid() << std::endl;
  */

  CGAL::Path_on_surface<LCC_3_cmap> path1(lcc1);

  generate_one_positive_spur(path1);
  std::cout<<"Posivite spur (2^6 1 0 2^4): "<<std::flush;
  path1.display_positive_turns();
  simplify_path(path1, false);
  std::cout<<" -> "; path1.display_positive_turns(); std::cout<<std::endl;

  generate_one_negative_spur(path1);
  std::cout<<"Negative spur (-2^6 -1 0 -2^4): "<<std::flush;
  path1.display_negative_turns();
  simplify_path(path1, false);
  std::cout<<" -> "; path1.display_negative_turns(); std::cout<<std::endl;

  generate_cyclic_spur(path1);
  std::cout<<"Cyclic spur (0 0): "<<std::flush;
  path1.display_positive_turns();
  simplify_path(path1, false);
  std::cout<<" -> "; path1.display_positive_turns(); std::cout<<std::endl;

  generate_one_positive_bracket(path1);
  std::cout<<"Positive bracket (2^3 3 1 2^6 1 3 2^2): "<<std::flush;
  path1.display_positive_turns();
  simplify_path(path1, false);
  std::cout<<" -> "; path1.display_positive_turns(); std::cout<<std::endl;

  generate_one_negative_bracket(path1);
  std::cout<<"Negative bracket (-2^3 -3 -1 -2^6 -1 -3 -2^2): "<<std::flush;
  path1.display_negative_turns();
  simplify_path(path1, false);
  std::cout<<" -> "; path1.display_negative_turns(); std::cout<<std::endl;

  CGAL::Path_on_surface<LCC_3_cmap> path2(lcc2);

  generate_positive_bracket_special1(path2);
  std::cout<<"Positive special case 1 (3 1 2^8 1): "<<std::flush;
  path2.display_positive_turns();
  simplify_path(path2, false);
  std::cout<<" -> +"; path2.display_positive_turns();
  std::cout<<" -"; path2.display_negative_turns(); std::cout<<std::endl;

  generate_negative_bracket_special1(path2);
  std::cout<<"Negative special case 1 (-3 -1 -2^8 -1): "<<std::flush;
  path2.display_negative_turns();
  simplify_path(path2, true);
  std::cout<<" -> +"; path2.display_positive_turns();
  std::cout<<" -"; path2.display_negative_turns(); std::cout<<std::endl;
}

int main(int argc, char** argv)
{
  // test_file(argc, argv);
  // test_some_random_paths_on_cube();
   test_all_cases_spurs_and_bracket();

  return EXIT_SUCCESS;
}
