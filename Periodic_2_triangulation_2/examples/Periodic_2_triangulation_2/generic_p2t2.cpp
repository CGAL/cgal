#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/internal/Generic_P2T2/Periodic_2_Delaunay_triangulation_2_generic.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K                                                   GT;
typedef GT::Vector_2                                        Vector;

typedef CGAL::Periodic_2_triangulation_vertex_base_2_generic<GT>    Vb;
typedef CGAL::Periodic_2_triangulation_face_base_2_generic<GT>      Fb;

typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
typedef CGAL::Periodic_2_Delaunay_triangulation_2_generic<GT, Tds>  PDT;

typedef PDT::Vertex_handle                                  Vertex_handle;
typedef PDT::Point                                          Point;

std::vector<Point> generate_lattice()
{
  std::vector<Point> pts;
  return pts;
}

int main()
{
  std::pair<Vector, Vector> basis;
  basis = std::make_pair(Vector(-0.5, 1), Vector(1.5, 0));

  std::vector<Point> pts { Point(0, 0), Point(1, 0) };
//  std::vector<Point> pts = generate_latice();

  PDT T(pts.begin(), pts.end(), basis);

  return EXIT_SUCCESS;
}
