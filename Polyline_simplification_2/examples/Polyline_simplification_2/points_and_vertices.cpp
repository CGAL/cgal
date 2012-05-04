#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Polyline_constrained_triangulation_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>

namespace PS = CGAL::Polyline_simplification_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;
typedef CGAL::Exact_predicates_tag                          Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,CGAL::Default, Itag> CDT;
typedef CGAL::Polyline_constrained_triangulation_2<CDT>     PCT;
typedef PCT::Point                           Point;
typedef PCT::Constraint_id                   Constraint_id;
typedef PCT::Constraint_iterator             Constraint_iterator;
typedef PCT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
typedef PCT::Points_in_constraint_iterator   Points_in_constraint_iterator;
typedef PS::Stop_below_count_ratio_threshold Stop;
typedef PS::Squared_distance_cost Cost;

void print(const PCT& pct, Constraint_id cid)
{
  std::cout << "simplified polyline" <<std::endl;
  for(Vertices_in_constraint_iterator vit = 
        pct.vertices_in_constraint_begin(cid);
      vit != pct.vertices_in_constraint_end(cid);
      ++vit){
    std::cout << vit->point() << std::endl ;
  }
  
  std::cout << "original points" <<std::endl;
  for(Points_in_constraint_iterator pit = 
        pct.vertices_in_constraint_begin(cid);
      pit != pct.vertices_in_constraint_end(cid);
      ++pit){
    std::cout << *pit << std::endl ;
  }
  
}


int main( )
{
  const bool keep_points = true;
  PCT pct;
  Polygon_2 P;
  Constraint_id cid;
  int largest = 0;
  while(std::cin >> P){
    Constraint_id cid2 = pct.insert_constraint(P);
    if(P.size() > largest){
      cid = cid2;
    }
  }

  PS::mark_vertices_unremovable(pct,cid);
  PS::simplify(pct, cid, Cost(), Stop(0.5), keep_points);
  print(pct, cid);
  PS::simplify(pct, cid, Cost(), Stop(0.5), keep_points);
  pct.remove_points_from_constraint(cid);
  print(pct, cid);  

  return 0;
}


