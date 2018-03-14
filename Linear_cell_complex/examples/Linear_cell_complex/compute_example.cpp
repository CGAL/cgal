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

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
typedef CGAL::Linear_cell_complex_for_generalized_map<2,3> LCC_3_gmap;

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
  pp1.bracket_flattening();
  std::cout<<"After bracket flattening, the path has "<<pp1.length()<<" darts."<<std::endl;
}

void test_simplify_path(const LCC_3_cmap& lcc,
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

  std::vector<const CGAL::Path_on_surface<LCC_3_cmap>*> v;
  v.push_back(&path);
  // display(lcc, v);

  CGAL::Path_on_surface<LCC_3_cmap>* prevp=&path;
  CGAL::Path_on_surface<LCC_3_cmap>* curp=NULL;
  do
  {
    curp=new CGAL::Path_on_surface<LCC_3_cmap>(*prevp);
    if (curp->bracket_flattening_one_step())
    { v.push_back(curp); }
    else
    { curp=NULL; }
    prevp=curp;
    // if (nbtest==1)
     // display(lcc, v);
  }
  while(curp!=NULL);
  
  if (draw)
  { display(lcc, v); }

  std::cout<<"[END] TEST "<<nbtest++<<"."<<std::endl;
}

void test_square()
{
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, "./cube-mesh-5-5.off"))
  {
    std::cout<<"PROBLEM reading file ./cube-mesh-5-5.off"<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout<<"Initial map: ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  CGAL::Random random(1); // fix seed

  // path 1: 2 straight darts to begin; 1 + turn; 6 straight; 1 + turn; 3 straight
  test_simplify_path(lcc, 2, 6, 3, random, false);

  // path 2: 3 straight darts to begin; 1 + turn; 8 straight; 1 + turn; 4 straight
  test_simplify_path(lcc, 3, 8, 4, random, false);

  // path 3: 5 straight darts to begin; 1 + turn; 12 straight; 1 + turn; 8 straight
  test_simplify_path(lcc, 5, 12, 8, random, false);

  // path 4: 5 straight darts to begin; 1 + turn; 12 straight; 1 + turn; 8 straight
  test_simplify_path(lcc, 5, 12, 8, random, false);
}

int main(int argc, char** argv)
{
  // test_file(argc, argv);
  test_square();
   
  return EXIT_SUCCESS;
}
