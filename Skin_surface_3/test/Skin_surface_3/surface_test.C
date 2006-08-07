#define  CGAL_SKIN_SURFACE_USE_EXACT_IMPLICIT_SURFACE

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/mesh_skin_surface_3.h>

#include <list>
#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel   K;
typedef CGAL::Mixed_complex_traits_3<K>                     Traits;
typedef CGAL::Skin_surface_3<Traits>                        Skin_surface_3;
typedef Skin_surface_3::RT                                  RT;
typedef Skin_surface_3::Weighted_point                      Weighted_point;
typedef Weighted_point::Point                               Bare_point;
typedef CGAL::Polyhedron_3<Traits>                          Polyhedron;
typedef Skin_surface_3::Quadratic_surface                   Quadratic_surface;


typedef Skin_surface_3::Triangulated_mixed_complex  Triangulated_mixed_complex;
typedef Triangulated_mixed_complex::Vertex_handle   Tmc_Vertex_handle;
typedef Triangulated_mixed_complex::Finite_vertices_iterator 
                                                    Tmc_Finite_vertices_iterator;
typedef Triangulated_mixed_complex::Finite_cells_iterator 
                                                    Tmc_Finite_cells_iterator;

#include <fstream>

class Test_file {
public:
  Test_file(double shrink) : s(shrink) {
  }
  void operator()(char * filename) {
    std::cout << s << " " << filename << std::endl;
    
    std::ifstream in(filename);
    
    std::list<Weighted_point> l;
    double x,y,z,w;
    while (in >> x >> y >> z >> w) 
      l.push_front(Weighted_point(Bare_point(x,y,z),w));
    
    Skin_surface_3 skin_surface(l.begin(), l.end(), s);
    
    Skin_surface_3::Triangulated_mixed_complex tmc;
    triangulate_mixed_complex_3(skin_surface.get_regular_triangulation(), 
				skin_surface.get_shrink_factor(), 
				tmc, false);

    for (Tmc_Finite_vertices_iterator vit = tmc.finite_vertices_begin();
	 vit != tmc.finite_vertices_end(); vit++) {
      if (tmc.is_infinite(vit->cell())) {
	std::cerr << "ERROR: is_infinite (main)" << std::endl;
      }
      Quadratic_surface::RT val = vit->cell()->surf->value(vit->point());
      std::list<Triangulated_mixed_complex::Cell_handle> cells;
      tmc.incident_cells(vit, std::back_inserter(cells));
      for (std::list<Triangulated_mixed_complex::Cell_handle>::iterator cell = 
	     cells.begin();
	   cell != cells.end(); cell++) {
	if (!tmc.is_infinite(*cell)) {
	  Quadratic_surface::RT val2 = (*cell)->surf->value(vit->point());
 	  CGAL_assertion(val == val2);
	}
      }
    }

    for (Tmc_Finite_cells_iterator cit = tmc.finite_cells_begin();
	 cit != tmc.finite_cells_end(); cit++) {
      Bare_point baryc = (cit->vertex(0)->point() + 
			  (cit->vertex(1)->point()-cit->vertex(0)->point())/4
			  +
			  (cit->vertex(2)->point()-cit->vertex(0)->point())/4
			  +
			  (cit->vertex(3)->point()-cit->vertex(0)->point())/4);
      Quadratic_surface::RT val1 = cit->surf->value(baryc);
      Quadratic_surface::RT val2 = skin_surface.value(baryc);
      CGAL_assertion(val1==val2);
    }

  }
private:
  double s;
};

int main(int argc, char *argv[]) {
  std::vector<char *> filenames;
   filenames.push_back("data/degenerate.cin");
   filenames.push_back("data/test1.cin");
   filenames.push_back("data/test2.cin");
   filenames.push_back("data/test3.cin");
   filenames.push_back("data/test4.cin");
   filenames.push_back("data/test5.cin");
   filenames.push_back("data/test6.cin");
   filenames.push_back("data/test7.cin");
   filenames.push_back("data/test8.cin");
   filenames.push_back("data/test9.cin");
   filenames.push_back("data/test10.cin");
   filenames.push_back("data/test11.cin");
  
  std::for_each(filenames.begin(), filenames.end(), Test_file(.85));
  std::for_each(filenames.begin(), filenames.end(), Test_file(.5));

  return 0;
}
