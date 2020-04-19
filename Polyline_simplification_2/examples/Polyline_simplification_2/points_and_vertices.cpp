#include <iostream>
#include <fstream>
#include <boost/config.hpp>
#include <boost/version.hpp>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
#include <CGAL/IO/WKT.h>
#endif

namespace PS = CGAL::Polyline_simplification_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Polygon_2<K>    Polygon_2;
typedef CGAL::Polygon_with_holes_2<K>    Polygon_with_holes_2;
typedef PS::Vertex_base_2<K>  Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Exact_predicates_tag                          Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS, Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>     CT;
typedef CT::Point                           Point;
typedef CT::Constraint_id                   Constraint_id;
typedef CT::Constraint_iterator             Constraint_iterator;
typedef CT::Vertices_in_constraint_iterator Vertices_in_constraint_iterator;
typedef CT::Points_in_constraint_iterator   Points_in_constraint_iterator;
typedef PS::Stop_below_count_ratio_threshold Stop;
typedef PS::Squared_distance_cost Cost;

void print(const CT& ct, Constraint_id cid)
{
  std::cout << "simplified polyline" <<std::endl;
  for(Vertices_in_constraint_iterator vit =
        ct.vertices_in_constraint_begin(cid);
      vit != ct.vertices_in_constraint_end(cid);
      ++vit){
    std::cout << (*vit)->point() << std::endl ;
  }

  std::cout << "original points" <<std::endl;
  for(Points_in_constraint_iterator pit =
        ct.points_in_constraint_begin(cid);
      pit != ct.points_in_constraint_end(cid);
      ++pit){
    std::cout << *pit << std::endl ;
  }

}


int main(int argc, char* argv[])
{
  std::ifstream ifs( (argc==1)?"data/polygon.wkt":argv[1]);
#if BOOST_VERSION >= 105600 && (! defined(BOOST_GCC) || BOOST_GCC >= 40500)
  const bool remove_points = false;
  CT ct;
  Polygon_with_holes_2 P;
  Constraint_id cid;
  std::size_t largest = 0;
  while(CGAL::read_polygon_WKT(ifs, P)){
    const Polygon_2& poly = P.outer_boundary();
    Constraint_id cid2 = ct.insert_constraint(poly);
    if(poly.size() > largest){
      cid = cid2;
    }
  }

  PS::simplify(ct, cid, Cost(), Stop(0.5), remove_points);
  print(ct, cid);
  PS::simplify(ct, cid, Cost(), Stop(0.5), remove_points);
  ct.remove_points_without_corresponding_vertex(cid);
  print(ct, cid);
#else
  ifs.close();
#endif
  return 0;
}


