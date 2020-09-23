#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel                 K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K>                     CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>                       CDTP;

typedef CDTP::Point                                                       Point;
typedef CDTP::Constraint_id                                               Cid;
typedef CDTP::Vertex_handle                                               Vertex_handle;

int main()
{
  CDTP cdtp;

  // First polyline
  cdtp.insert_constraint(Point(0,0), Point(1,1));
  cdtp.insert_constraint(Point(1,1), Point(2,2));
  cdtp.insert_constraint(Point(2,2), Point(3,3));

  // Second polyline that cuts the first one in 2
  cdtp.insert_constraint(Point(1,1), Point(1,2));
  cdtp.insert_constraint(Point(1,2), Point(1,3));

  // Third polyline
  cdtp.insert_constraint(Point(10,10), Point(20,20));
  cdtp.insert_constraint(Point(20,20), Point(30,30));

  // Fourth polyline
  cdtp.insert_constraint(Point(100,100), Point(200,200));

  // Segment soup of 8 segments as input
  std::cout << "Input CDT+ has " << cdtp.number_of_constraints() << " constraints/subconstraints" << std::endl;

  std::cout << "Splitting subconstraints graph into constraints" << std::endl;
  cdtp.split_subconstraint_graph_into_constraints();

  // 5 polylines as output
  std::cout << "Output CDT+ has " << cdtp.number_of_constraints() << " constraints:" << std::endl;

  for (CDTP::Constraint_id cid : cdtp.constraints())
  {
    std::cout << " *";
    for (CDTP::Vertex_handle vh : cdtp.vertices_in_constraint(cid))
      std::cout << " (" << vh->point() << ")";
    std::cout << std::endl;
  }

  return EXIT_SUCCESS;
}
