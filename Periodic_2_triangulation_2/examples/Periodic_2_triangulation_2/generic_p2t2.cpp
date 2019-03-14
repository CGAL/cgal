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

  for(int i=0; i<10; ++i)
    for(int j=0; j<10; ++j)
      pts.emplace_back(i * 0.14, j * 0.32);

  return pts;
}

int main()
{
  CGAL::cpp11::array<Vector, 2> basis;
  basis = CGAL::make_array(Vector(-0.5, 1), Vector(1.5, 0));

#if 0
  std::vector<Point> pts { Point(0, 0), Point(-0.2, -0.6) };
#else
  std::vector<Point> pts = generate_lattice();
#endif

  PDT T(pts.begin(), pts.end(), basis);

//  PDT::Vertex_handle vh = T.insert(Point(12345.6, 5432.1));
//  std::set<PDT::Face_handle> incident_faces = T.incident_faces(vh);

  draw(T.dt2, "Post-Insertion");

//  T.convert_to_p2t2();
//  T.p2t2.insert(Point(0.3, 0.12));

  return EXIT_SUCCESS;
}
