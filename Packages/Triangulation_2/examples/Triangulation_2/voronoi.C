#include <CGAL/basic.h>
#include <fstream>


// Define shorter names to please linker (g++/egcs)
#define Cartesian Cart

#include <CGAL/Cartesian.h>

#include <CGAL/squared_distance_2.h>   // to avoid a g++ problem
#include <CGAL/Point_2.h>
#include <CGAL/predicates_on_points_2.h>
#include <CGAL/Triangle_2.h>
#include <CGAL/Segment_2.h>

#include <CGAL/Triangulation_euclidean_traits_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

typedef double coord_type;
typedef CGAL::Cartesian<coord_type>  Rpst;

typedef CGAL::Triangulation_euclidean_traits_2<Rpst> Gt;
typedef Gt::Point_2       Point;
typedef Gt::Segment_2     Segment;
typedef Gt::Ray_2         Ray;
typedef Gt::Line_2        Line;
typedef Gt::Triangle_2    Triangle;

typedef CGAL::Triangulation_vertex_base_2<Gt> Vb;
typedef CGAL::Triangulation_face_base_2<Gt>  Fb;
typedef CGAL::Triangulation_default_data_structure_2<Gt,Vb,Fb> Tds;
typedef CGAL::Delaunay_triangulation_2<Gt,Tds>  Triangulation;

typedef Triangulation::Face_iterator  Face_iterator;
typedef Triangulation::Edge_iterator  Edge_iterator;

int
main( )
{
  //    Options opt;
  //  parse(argc, argv, opt);

  std::ifstream in("data/voronoi.cin");
  std::istream_iterator<Point> begin(in);
  std::istream_iterator<Point> end;
  Triangulation T;
  T.insert(begin, end);
    
  std::cout << "Vertices :" << std::endl << "=============" << std::endl;
  Face_iterator fit =T.faces_begin()  ; 
  for (  ; fit != T.faces_end(); ++fit) 
    std::cout << T.dual(fit) << std::endl;

  std::cout << "Segment :" << std::endl << "=============" << std::endl;
  {
    Edge_iterator eit =T.edges_begin();
    for ( ; eit !=T.edges_end(); ++eit) {
 	CGAL::Object o = T.dual(eit);
	Triangulation::Segment s;
	if (CGAL::assign(s,o)) std::cout << s << std::endl;
    }
  }

  std::cout << "Rays :" << std::endl << "=============" << std::endl;
  {
     Edge_iterator eit =T.edges_begin();
    for ( ; eit !=T.edges_end(); ++eit) {
      CGAL::Object o = T.dual(eit);
      Ray r;
      if (CGAL::assign(r,o)) std::cout << r << std::endl;
    }
  }

    return 0;
}
