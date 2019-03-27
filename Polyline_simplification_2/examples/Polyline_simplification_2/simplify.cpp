#include <iostream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Polyline_simplification_2/Squared_distance_cost.h>

namespace PS = CGAL::Polyline_simplification_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>                                  Polygon_2;

typedef PS::Vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS, CGAL::Exact_predicates_tag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>     CT;
typedef CT::Point                           Point;
typedef CT::Constraint_iterator             Constraint_iterator;
typedef CT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
typedef CT::Points_in_constraint_iterator   Points_in_constraint_iterator;
typedef PS::Stop_below_count_ratio_threshold Stop;
typedef PS::Squared_distance_cost Cost;

int main()
{
  CT ct;
  Polygon_2 P;
  while(std::cin >> P){
    ct.insert_constraint(P);
  }
  PS::simplify(ct, Cost(), Stop(0.5));

  for(Constraint_iterator cit = ct.constraints_begin();
      cit != ct.constraints_end();
      ++cit) {
    std::cout << "simplified polyline" << std::endl;
    for(Points_in_constraint_iterator vit = 
          ct.points_in_constraint_begin(*cit);
        vit != ct.points_in_constraint_end(*cit);
        ++vit)
      std::cout << *vit << std::endl;
  }
  return 0;
}


