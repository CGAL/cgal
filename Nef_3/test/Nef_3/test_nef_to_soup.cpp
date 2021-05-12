#include <CGAL/boost/graph/convert_nef_polyhedron_to_polygon_mesh.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Nef_polyhedron_3.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Nef_polyhedron_iostream_3.h>

#include <fstream>

typedef CGAL::Exact_predicates_exact_constructions_kernel EPEC;
typedef CGAL::Exact_predicates_inexact_constructions_kernel EPIC;

int main()
{
  typedef CGAL::Nef_polyhedron_3< EPEC > Nef_polyhedron;
  typedef CGAL::Polyhedron_3< EPEC >     Polyhedron;
  typedef CGAL::Polyhedron_3< EPIC >     PolygonMesh;
  typedef typename EPIC::Point_3         Point;
  typedef std::list<std::list<std::size_t> > PolygonRange;
  typedef std::list<Point> PointRange;

  typename EPEC::RT n, d;
  std::istringstream str_n("6369051672525773");
  str_n >> n;
  std::istringstream str_d("4503599627370496");
  str_d >> d;

  EPEC::Point_3 p(n, 0, 0, d);
  EPEC::Point_3 q(0, n, 0, d);
  EPEC::Point_3 r(0, 0, n, d);
  EPEC::Point_3 s(0, 0, 0, 1);

  std::cout << "    build...\n";
  Polyhedron P;
  P.make_tetrahedron( p, q, r, s);
  Nef_polyhedron nef( P );
  PointRange points;
  PolygonRange polygons;
    std::cout << "    convert...\n";
  CGAL::convert_nef_polyhedron_to_polygon_soup(nef, points, polygons);
  CGAL_assertion(points.size() == 4);
  CGAL_assertion(polygons.size() == 4);
  PolygonMesh pm;
  CGAL::convert_nef_polyhedron_to_polygon_mesh<
      Nef_polyhedron,
      PolygonMesh>(nef, pm);

  return 0;
}
