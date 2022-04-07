#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb>           TDS;
typedef CGAL::Exact_predicates_tag                               Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, Itag> Triangulation;
typedef CGAL::Constrained_triangulation_plus_2<Triangulation> Delaunay;

typedef K::Point_2 Point_2;

int main()
{
  Delaunay dt;

  dt.insert(Point_2(0,0));
  dt.insert(Point_2(10,0));
  dt.insert(Point_2(0,10));
  dt.insert(Point_2(10,10));

  std::vector<Point_2> vec_constraint;

  vec_constraint.push_back(Point_2(1, 2));
  vec_constraint.push_back(Point_2(2, 2));
  vec_constraint.push_back(Point_2(3, 2));

  dt.insert_constraint(vec_constraint.begin(), vec_constraint.end());

  dt.insert_constraint(vec_constraint.begin(), vec_constraint.end());

  return 0;
}
