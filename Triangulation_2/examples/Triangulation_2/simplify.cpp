#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polyline_constrained_triangulation_2.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_vertex_base_2<K>                Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K>      Fb;
typedef CGAL::Triangulation_data_structure_2<Vb,Fb>         TDS;
typedef CGAL::Exact_predicates_tag                          Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS, Itag> CDT;
typedef CGAL::Polyline_constrained_triangulation_2<CDT>     PCT;
typedef PCT::Point                                          Point;
typedef PCT::Vertices_in_constraint_iterator                Vertices_in_constraint_iterator;
typedef PCT::Points_in_constraint_iterator                  Points_in_constraint_iterator;
CGAL::Polyline_simplification_2::Stop_below_count_ratio_threshold Stop;
typedef CGAL::Polyline_simplification_2::Squared_distance_cost<double> Cost;

void print(const PCT& pct)
{
  for(Vertices_in_constraint_iterator vit = pct.vertices_in_constraint_begin();
      vit != pct.vertices_in_constraint_end();
      ++vit){
    std::cout << (*vit)->point() << std::endl ;
  }
  for(Points_in_constraint_iterator pit = pct.vertices_in_constraint_begin();
      pit != pct.vertices_in_constraint_end();
      ++pit){
    std::cout << *pit << std::endl ;
  }
}

int main( )
{
  Point points[] = { Point(0,1)
                   , Point(1,2)
                   , Point(2,1)
                   , Point(3,3)
                   , Point(4,1)
                   , Point(5,4)
                   , Point(6,1)
                   , Point(3,0)
                   };
  PCT pct;
  pct.insert_polyline(points,points+8);
  
  Stop stop(0.5);
  Cost cost;
  
  pct.simplify(stop, cost);
  print(pct);

  pct.remove_points();
  print(pct);  

  return 0;
}


