#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constraint_hierarchy_2.h>
#include <CGAL/Polyline_constrained_triangulation_2.h>
#include <CGAL/Triangulation_conformer_2.h>
#include <list>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K>              Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>    Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>       TDS;
typedef CGAL::Exact_predicates_tag                        Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS, Itag> CDT;
typedef CGAL::Polyline_constrained_triangulation_2<CDT>   PCT;
typedef PCT::Point                                    Point;
typedef PCT::Constraint_id                            Cid;
typedef PCT::Vertices_in_constraint                   Vertices_in_constraint;
typedef PCT::Vertex_handle                            Vertex_handle;
typedef PCT::Face_handle                              Face_handle;
typedef PCT::Edge                                     Edge;

void 
print(const PCT& cdt, Cid cid)
{
  std::cout << "Polyline constraint:" << std::endl;
  for(Vertices_in_constraint it = cdt.vertices_in_constraint_begin(cid);
      it != cdt.vertices_in_constraint_end(cid);
      it++){
    Vertex_handle vh = it->vertex;
    std::cout << vh->point() << std::endl;
  }
}

int
main( )
{
  PCT cdt;

  std::vector<Point> points;
  points.push_back(Point(1,1));
  points.push_back(Point(5,2));
  points.push_back(Point(15,3));
  Cid id1 = cdt.insert_polyline(points.begin(), points.end());

  print(cdt, id1);

  points.clear();
  points.push_back(Point(15,3));
  points.push_back(Point(15,7));
  points.push_back(Point(7,7));

  Cid id2 = cdt.insert_polyline(points.begin(), points.end());
  print(cdt, id2);

  cdt.concatenate(id1, id2);
  print(cdt, id1);

  return 0;
}
