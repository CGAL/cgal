#include <CGAL/Gmpz.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>
#include <CGAL/convex_hull_3.h>
#include <list>
#include <fstream>

#include <CGAL/Nef_S2/Normalizing.h>
#include <CGAL/IO/Polyhedron_iostream.h>

typedef CGAL::Gmpz NT;
typedef CGAL::Homogeneous<NT> Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef CGAL::Nef_polyhedron_3<Kernel> Nef_polyhedron_3;
typedef Nef_polyhedron_3::Vertex_const_iterator Vertex_const_iterator;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Line_3 Line_3;
typedef Kernel::Aff_transformation_3 Aff_transformation_3;

typedef std::list<Plane_3> Plane_list;
typedef std::list<Plane_list> Hsp;
typedef Plane_list::const_iterator Plane_list_iterator;
typedef Hsp::const_iterator Hsp_iterator;
Hsp hsp;


void loadHsp() {

  int hulls;
  std::cin >> hulls;
  
  for(int h=0; h<hulls; ++h) {
    std::cerr << "read hull " << h << std::endl;
    Plane_3 pl;
    Plane_list plist;

    int planes;
    std::cin >> planes;
    for(int p=0; p<planes; ++p) {
      std::cerr << "read plane " << p << std::endl;
      std::cin >> pl;
      plist.push_back(pl);
    }
    hsp.push_back(plist);
  }
}

Nef_polyhedron_3 create_from_halfspaces(const Plane_list& plane_list, bool invert) {
  
  Nef_polyhedron_3 cube;
  std::ifstream in("Nef3/centered_cube.nef3");
  in >> cube;

  int mag = 1000000;
  Aff_transformation_3 scale(mag, 0, 0, 
			     0, mag, 0, 
			     0, 0, mag, 1);
  cube.transform(scale);

  Plane_list_iterator pli;
  for(pli = plane_list.begin(); pli != plane_list.end(); ++pli) {
    if(invert) {      
      Plane_3 inverted_plane(-pli->a(), -pli->b(),
			     -pli->c(), -pli->d());
      cube = cube.intersection(inverted_plane);
    } else
      cube = cube.intersection(*pli);
  }  
  return cube;
}

bool test_convex_hull(const Nef_polyhedron_3& N, const Plane_list& plane_list, bool inverse) {

  Plane_list_iterator pli;
  for(pli = plane_list.begin(); pli != plane_list.end(); ++pli) {
    //    std::cerr << "plane " << *pli << std::endl;
    Vertex_const_iterator vi;
    for(vi = N.vertices_begin(); vi != N.vertices_end(); ++vi) {
      //      std::cerr << CGAL::to_double(vi->point().x()) << ", " 
      //		<< CGAL::to_double(vi->point().y()) << ", " 
      //		<< CGAL::to_double(vi->point().z()) << ": " 
      //		<< pli->oriented_side(vi->point()) << std::endl; 
      if((!inverse  && pli->oriented_side(vi->point()) == CGAL::ON_POSITIVE_SIDE) ||
	 (inverse && pli->oriented_side(vi->point()) == CGAL::ON_NEGATIVE_SIDE))
	return false;
    }
  }

  return true;
}

int main(int argc, char* argv[]) {

  Nef_polyhedron_3 result;

  loadHsp();
  Hsp_iterator hi;
  int i = 0;
  for(hi = hsp.begin(); hi != hsp.end(); ++hi) {
    std::cerr << "polyhedron " << ++i << std::endl;
    if(hi == hsp.begin()) {
      std::cerr << "create 1 " << std::endl;
      result = create_from_halfspaces(*hi, false);
      std::cerr << "test 1 " << std::endl;
      if(!test_convex_hull(result, *hi, false))
	std::cerr << "convex hull incorrect" << std::endl;
    } else {
      std::cerr << "create 2 " << std::endl;
      Nef_polyhedron_3 tmp = create_from_halfspaces(*hi, false);
      std::cerr << "size of obstacle " << tmp.number_of_vertices() << std::endl;
      if(!test_convex_hull(tmp, *hi, false))
	std::cerr << "convex hull incorrect" << std::endl;
      result = result - tmp;
    }
  }

  std::cout << result;
}
