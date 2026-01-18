#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polyline_simplification_2/simplify.h>
#include <vector>

namespace PS = CGAL::Polyline_simplification_2;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef PS::Vertex_base_2<K> Vb;
typedef CGAL::Constrained_triangulation_face_base_2<K> Fb;
typedef CGAL::Triangulation_data_structure_2<Vb, Fb> TDS;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, TDS , CGAL::Exact_predicates_tag> CDT;
typedef CGAL::Constrained_triangulation_plus_2<CDT>     CT;
typedef CT::Point                            Point;
typedef PS::Stop_above_cost_threshold        Stop;
typedef PS::Squared_distance_cost            Cost;

int main()
{
  double tolerance = 100;
  CT ct;
  std::vector<CT::Point> pts;

  pts.push_back(CT::Point(0, 0));
  pts.push_back(CT::Point(2, 0));
  pts.push_back(CT::Point(1, 0));
  pts.push_back(CT::Point(3, 0));
  pts.push_back(CT::Point(4, 1));
  tolerance = 100;
  ct.insert_constraint(pts.begin(), pts.end(), false);

  PS::simplify(ct, Cost(), Stop(tolerance * tolerance), false);

  return 0;
}
