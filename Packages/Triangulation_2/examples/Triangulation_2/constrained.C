// file : examples/Triangulation_2/constrained.C
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <fstream>
#include <list>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>

typedef double Coord_type;
typedef CGAL::Cartesian<Coord_type>  Gt;
typedef Gt::Point_2    Point;
typedef CGAL::Constrained_Delaunay_triangulation_2<Gt>  CDT;
typedef CDT::Constraint     Constraint;

int
main( )
{
  std::ifstream is("data/constrained.cin");
  std::list<Constraint> lc;
  int n;
  Point p,q;
  is >> n;
  std::cerr << "Reading " << n << " constraints" << std::endl;
  for(; n > 0; n--) {
    is >> p >> q;
    lc.push_back(std::make_pair(p,q));
  }

  CDT cdt(lc.begin(),lc.end());
  assert(cdt.is_valid());
  return 0;
}
  
