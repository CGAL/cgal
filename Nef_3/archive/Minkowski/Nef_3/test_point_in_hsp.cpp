#include <CGAL/Cartesian.h>
#include <list>
#include <fstream>

typedef CGAL::Cartesian<double> Kernel;
typedef Kernel::Plane_3 Plane_3;
typedef Kernel::Point_3 Point_3;

typedef std::list<Plane_3> Plane_list;
typedef std::list<Plane_list> Hsp;
typedef Plane_list::const_iterator Plane_list_iterator;
typedef Hsp::const_iterator Hsp_iterator;
Hsp hsp;


void loadHsp() {

  int hulls;
  std::cin >> hulls;

  for(int h=0; h<hulls; ++h) {
    //    std::cerr << "read hull " << h << std::endl;
    Plane_3 pl;

    for(int j=0; j<2; ++j) {
      Plane_list plist;
      int planes;
      std::cin >> planes;
      for(int p=0; p<planes; ++p) {
        //      std::cerr << "read plane " << p << std::endl;
        std::cin >> pl;
        plist.push_back(pl);
      }
      hsp.push_back(plist);
    }
  }
}

bool test_point_in_polyhedron(const Point_3 p,
                              const Plane_list& bad,
                              const Plane_list& good,
                              bool turn) {

  bool result = true;
  Plane_list_iterator pli;
  for(pli = bad.begin(); pli != bad.end(); ++pli) {
    CGAL::Oriented_side os = pli->oriented_side(p);
    if(os == CGAL::ON_ORIENTED_BOUNDARY)
      std::cerr << "attention: point on boundary " << *pli << std::endl;
    if(turn) {
      if(os != CGAL::ON_NEGATIVE_SIDE) continue;
    } else {
      if(os != CGAL::ON_POSITIVE_SIDE) continue;
    }
    std::cerr << "fail " << *pli << std::endl;
    result = false;
  }

  for(pli = good.begin(); pli != good.end(); ++pli) {
    CGAL::Oriented_side os = pli->oriented_side(p);
    if(turn) {
      if(os == CGAL::ON_POSITIVE_SIDE) continue;
    } else {
      if(os == CGAL::ON_NEGATIVE_SIDE) continue;
    }
    if(os == CGAL::ON_ORIENTED_BOUNDARY)
      std::cerr << "on csp boundary " << std::endl;
    std::cerr << "fail " << *pli << std::endl;
    result = false;
  }

  return result;
}

int main(int argc, char* argv[]) {

  CGAL_assertion(argc > 3 && argc < 7);
  double x = std::atof(argv[1]);
  double y = std::atof(argv[2]);
  double z = std::atof(argv[3]);
  Point_3 p(x,y,z);

  std::cerr << " test point " << p << std::endl;

  bool turn_hull = (argc>4 && std::atoi(argv[4]) == 1) ? true : false;
  bool turn_obs  = (argc>5 && std::atoi(argv[5]) == 1) ? true : false;

  std::cerr << " turn hull " << turn_hull << std::endl;
  std::cerr << " turn obs  " << turn_obs << std::endl;

  loadHsp();

  CGAL_assertion((hsp.size()%2) == 0);
  Plane_list emptyList;

  int problem = 0;
  Hsp_iterator hi;
  int i = 0;
  for(hi = hsp.begin(); hi != hsp.end(); ++hi) {
    std::cerr << "polyhedron " << ++i << std::endl;
    const Plane_list& bad = *hi;
    if(hi == hsp.begin()) {
      ++hi;
      CGAL_assertion(hi->size() == 0);
      if(!test_point_in_polyhedron(p, bad, emptyList, turn_hull)) {
        ++problem;
        std::cerr << "point not in hull " << std::endl;
      }
    } else {
      ++hi;
      const Plane_list& good = *hi;
      if(test_point_in_polyhedron(p, bad, good, turn_obs)) {
        ++problem;
        std::cerr << "point inside obstacle " << std::endl;
      }
    }
  }

  CGAL_assertion(problem < 2);

  if(problem)
    std::cerr << "point NOT in hsp " << std::endl;
  else
    std::cerr << "point in hsp " << std::endl;
}
