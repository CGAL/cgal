#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <CGAL/Timer.h>

namespace PS = CGAL::Polyline_simplification_2;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef PS::Vertex_base_2<K>  Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Exact_predicates_tag                          Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K,TDS, Itag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>     CT;
typedef CGAL::Polygon_2<K>                   Polygon_2;
typedef PS::Stop_above_cost_threshold Stop;
typedef PS::Squared_distance_cost            Cost;



void test(char* fname)
{
  CGAL::Timer timer;
  std::cerr << "simplify " << fname << std::endl;
  std::ifstream in(fname);
  int n;
  CT ct;
  Polygon_2 P;
  in >> n; // number of polygons
  while(in >> P){
    ct.insert_constraint(P);
  }
  std::cerr << ct.number_of_vertices() << " vertices\n";
  timer.start();
  PS::simplify(ct, Cost(), Stop(0.5));
  std::cerr << ct.number_of_vertices() << " vertices\n";
  std::cerr << timer.time() << " sec.\n";
}

int main(int argc, char* argv[])
{

  for(int i= 1;i < argc; i++){
    test(argv[i]);
  }

  return 0;
}




