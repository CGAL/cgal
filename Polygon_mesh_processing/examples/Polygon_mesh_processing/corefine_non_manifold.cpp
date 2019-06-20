#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/stitch_borders.h>
#include <CGAL/Polygon_mesh_processing/triangulate_faces.h>

#include <CGAL/boost/graph/helpers.h>

#include <fstream>
#include <iostream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef K::Point_3 P;
typedef CGAL::Surface_mesh<P> Surface_mesh;

namespace PMP = CGAL::Polygon_mesh_processing;

int main()
{
// {
//   Surface_mesh tm;
//   CGAL::make_triangle (P(0,0,0), P(5,0,0), P(0,5,0), tm);
//   CGAL::make_triangle (P(1,1,-1), P(2,1,-1), P(1.5,1,1), tm);
// std::ofstream("/tmp/tm.off") << tm;
//   PMP::experimental::autorefine(tm);
// std::ofstream("/tmp/tm_aurefined.off") << tm;
// }

  Surface_mesh tm1;
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,4,0), P(0,4,0), tm1);
  CGAL::make_quad(P(0,-4,0), P(4,-4,0), P(4,0,0), P(0,0,0), tm1);
  PMP::stitch_borders(tm1);
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,0,4), P(0,0,4), tm1);
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,0,-4), P(0,0,-4), tm1);
  PMP::triangulate_faces(tm1);

  Surface_mesh tm2;
//  CGAL::make_quad(P(2,2,-2), P(2,-2,-2), P(2,-2,2), P(2,2,2), tm2);  //TODO: test me + also test splitting the diagonal
// TODO: test more case including non-manifold vertices
  CGAL::make_quad(P(2,3,-2), P(2,-2,-2), P(2,-2,2), P(2,3,2), tm2);
  PMP::triangulate_faces(tm2);

  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_1(tm1, get(boost::vertex_point, tm1));

  // std::vector< std::vector<P> > polylines;
  // PMP::surface_intersection(tm1, tm2, std::back_inserter(polylines), CGAL::parameters::non_manifold_feature_map(nm_map_1));
  //
  // //dump polylines
  // std::ofstream output("intersection_polylines.cgal");
  // for(const std::vector<P>& polyline : polylines)
  // {
  //   output << polyline.size() << " ";
  //   std::copy(polyline.begin(), polyline.end(),std::ostream_iterator<P>(output," "));
  //   output << "\n";
  // }

  PMP::corefine(tm1, tm2, CGAL::parameters::non_manifold_feature_map(nm_map_1));

  std::ofstream("tm1_corefined.off") << tm1;
  std::ofstream("tm2_corefined.off") << tm2;
}
