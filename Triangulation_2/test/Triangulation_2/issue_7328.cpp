#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel EPECK;

typedef CGAL::Exact_intersections_tag                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<EPECK, CGAL::Default, Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>       CDTplus;
typedef CDT::Vertex_handle                                Vertex_handle;

int main() {
  CDTplus cdt;
  std::vector<Vertex_handle> vertices;
  vertices.push_back(cdt.insert(EPECK::Point_2(0.867732088229836496, -1.098635766973843797)));
  vertices.push_back(cdt.insert(EPECK::Point_2(0.868834588233415861, -1.100000000000000533)));
  vertices.push_back(cdt.insert(EPECK::Point_2(0.729063637498132522, -0.927047486193771419)));
  vertices.push_back(cdt.insert(EPECK::Point_2(0.760857518227448626, -0.918203415668045420)));

  cdt.insert_constraint(vertices[0], vertices[2]);
  cdt.insert_constraint(vertices[0], vertices[3]);

  CDTplus cdtC = cdt;

  EPECK::Point_2 p(0.868834588233415861, -1.100000000000000533), q(0.729063637498132522, -0.927047486193771419);
  cdtC.insert_constraint(p,q);
  return 0;
}
