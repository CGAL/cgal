#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>

#include <vector>

typedef CGAL::Exact_predicates_exact_constructions_kernel                 K;
typedef CGAL::Polygon_2<K>                                                Polygon_2;
typedef CGAL::Exact_intersections_tag                                     Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,CGAL::Default, Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>                       CDTP;

typedef CDTP::Point                                                       Point;
typedef CDTP::Constraint_id                                               Cid;
typedef CDTP::Vertex_handle                                               Vertex_handle;

void
print(const CDTP& cdtp, Cid cid)
{
  std::cout << "Polyline constraint:" << std::endl;
  for(Vertex_handle vh : cdtp.vertices_in_constraint(cid)){
    std::cout << vh->point() << std::endl;
  }
}


void
contexts(const CDTP& cdtp)
{
  for(auto sc : cdtp.subconstraints()){
    Vertex_handle vp = sc.first.first, vq = sc.first.second;

    if(cdtp.number_of_enclosing_constraints(vp, vq) == 2){
      std::cout << "subconstraint " << vp->point() << " " << vq->point()
                << " is on constraints starting at:\n";
      for(const CDTP::Context& c : cdtp.contexts(vp,vq)){
        std::cout << (*(c.vertices_begin()))->point() << std::endl;
      }
    }
  }
}

int
main( )
{
  CDTP cdtp;

  cdtp.insert_constraint(Point(0,0), Point(1,1));

  std::vector<Point> points;
  points.push_back(Point(1,1));
  points.push_back(Point(5,2));
  points.push_back(Point(6,0));
  points.push_back(Point(3,0));
  Cid id1 = cdtp.insert_constraint(points.begin(), points.end());

  print(cdtp, id1);

  Polygon_2 poly;
  poly.push_back(Point(2,3));
  poly.push_back(Point(4,0));
  poly.push_back(Point(5,0));
  poly.push_back(Point(6,2));

  Cid id2 = cdtp.insert_constraint(poly.vertices_begin(), poly.vertices_end(), true);

  print(cdtp, id1);
  print(cdtp, id2);

  contexts(cdtp);

  return 0;
}
