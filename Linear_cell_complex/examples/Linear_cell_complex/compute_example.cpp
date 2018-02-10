#include <CGAL/Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/Linear_cell_complex_for_generalized_map.h>
#include <CGAL/Linear_cell_complex_constructors.h>
#include <CGAL/Combinatorial_map_functionalities.h>
#include <vector>

/* If you want to use a viewer, you can use qglviewer. */
#ifdef CGAL_USE_BASIC_VIEWER
#include <CGAL/LCC_with_paths.h>
#endif

typedef CGAL::Linear_cell_complex_for_combinatorial_map<2,3> LCC_3_cmap;
typedef CGAL::Linear_cell_complex_for_generalized_map<2,3> LCC_3_gmap;

int main(int argc, char** argv)
{
  if (argc!=2)
  {
    std::cout<<"Usage: "<<argv[0]<<" filename.off"<<std::endl;
    exit(EXIT_FAILURE);
  }
  
  LCC_3_cmap lcc;
  if (!CGAL::load_off(lcc, argv[1]))
  {
    std::cout<<"PROBLEM reading gile "<<argv[1]<<std::endl;
    exit(EXIT_FAILURE);
  }

  std::cout<<"Initial map: ";
  lcc.display_characteristics(std::cout) << ", valid=" 
                                         << lcc.is_valid() << std::endl;

  CGAL::Random random;
  CGAL::Path_on_surface<LCC_3_cmap> p1(lcc);
  p1.generate_random_path(10, random); // 10 %
  CGAL::Path_on_surface<LCC_3_cmap> p2(lcc);
  p2.generate_random_path(15, random);
  CGAL::Path_on_surface<LCC_3_cmap> p3(lcc);
  p3.generate_random_path(10, random); // 10 %
  CGAL::Path_on_surface<LCC_3_cmap> p4(lcc);
  p4.generate_random_path(15, random);
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

  pp1.bracket_flattening();

  return EXIT_SUCCESS;
}
