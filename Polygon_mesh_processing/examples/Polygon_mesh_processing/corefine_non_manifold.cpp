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
// polyline intersection with a non-manifold edge
{
  std::cout << "running polyline test\n";
  Surface_mesh tm1;
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,4,0), P(0,4,0), tm1);
  CGAL::make_quad(P(0,-4,0), P(4,-4,0), P(4,0,0), P(0,0,0), tm1);
  PMP::stitch_borders(tm1);
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,0,4), P(0,0,4), tm1);
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,0,-4), P(0,0,-4), tm1);
  PMP::triangulate_faces(tm1);

  Surface_mesh tm2;
  CGAL::make_quad(P(2,3,-2), P(2,-2,-2), P(2,-2,2), P(2,3,2), tm2);
  PMP::triangulate_faces(tm2);

  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_1(tm1, get(boost::vertex_point, tm1));

  std::vector< std::vector<P> > polylines;
  PMP::surface_intersection(tm1, tm2, std::back_inserter(polylines), CGAL::parameters::non_manifold_feature_map(nm_map_1));

  //dump polylines
  std::ofstream output("intersection_polylines.cgal");
  for(const std::vector<P>& polyline : polylines)
  {
    output << polyline.size() << " ";
    std::copy(polyline.begin(), polyline.end(),std::ostream_iterator<P>(output," "));
    output << "\n";
  }
}

// simple case with only one non-manifold edge
{
  std::cout << "running corefinement test 1\n";
  Surface_mesh tm1;
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,4,0), P(0,4,0), tm1);
  CGAL::make_quad(P(0,-4,0), P(4,-4,0), P(4,0,0), P(0,0,0), tm1);
  PMP::stitch_borders(tm1);
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,0,4), P(0,0,4), tm1);
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,0,-4), P(0,0,-4), tm1);
  PMP::triangulate_faces(tm1);

  Surface_mesh tm2;
  CGAL::make_quad(P(2,3,-2), P(2,-2,-2), P(2,-2,2), P(2,3,2), tm2);
  PMP::triangulate_faces(tm2);

  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_1(tm1, get(boost::vertex_point, tm1));

  PMP::corefine(tm1, tm2, CGAL::parameters::non_manifold_feature_map(nm_map_1));

  std::ofstream("t1_tm1_corefined.off") << tm1;
  std::ofstream("t1_tm2_corefined.off") << tm2;
}

// edge-edge intersection on a non-manifold edge
{
  std::cout << "running corefinement test 2\n";
  Surface_mesh tm1;
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,4,0), P(0,4,0), tm1);
  CGAL::make_quad(P(0,-4,0), P(4,-4,0), P(4,0,0), P(0,0,0), tm1);
  PMP::stitch_borders(tm1);
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,0,4), P(0,0,4), tm1);
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,0,-4), P(0,0,-4), tm1);
  PMP::triangulate_faces(tm1);

  Surface_mesh tm2;
  CGAL::make_quad(P(2,2,-2), P(2,-2,-2), P(2,-2,2), P(2,2,2), tm2);  //TODO: test me + also test splitting the diagonal
  PMP::triangulate_faces(tm2);

  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_1(tm1, get(boost::vertex_point, tm1));

  PMP::corefine(tm1, tm2, CGAL::parameters::non_manifold_feature_map(nm_map_1));

  std::ofstream("t2_tm1_corefined.off") << tm1;
  std::ofstream("t2_tm2_corefined.off") << tm2;
}

// coplanar edge and non-manifold edge
{
  std::cout << "running corefinement test 3\n";
  Surface_mesh tm1;
  CGAL::make_quad(P(0,0,0), P(8,0,0), P(8,8,0), P(0,8,0), tm1);
  CGAL::make_quad(P(8,0,0), P(12,0,0), P(12,8,0), P(8,8,0), tm1);
  PMP::stitch_borders(tm1);
  CGAL::make_quad(P(8,0,0), P(8,8,0), P(8,8,4), P(8,0,4), tm1);
  CGAL::make_quad(P(8,0,0), P(8,8,0), P(8,8,-4), P(8,0,-4), tm1);
  PMP::triangulate_faces(tm1);

  Surface_mesh tm2;
  CGAL::make_quad(P(6,1,0), P(14,3,0), P(13,7,4), P(5,5,0), tm2);
  CGAL::make_quad(P(6,1,0), P(5,5,0), P(5,5,2), P(6,1,3), tm2);
  PMP::stitch_borders(tm2);
  PMP::triangulate_faces(tm2);

  std::ofstream("t3_tm1.off") << tm1;
  std::ofstream("t3_tm2.off") << tm2;

  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_1(tm1, get(boost::vertex_point, tm1));

  PMP::corefine(tm1, tm2, CGAL::parameters::non_manifold_feature_map(nm_map_1));

  std::ofstream("t3_tm1_corefined.off") << tm1;
  std::ofstream("t3_tm2_corefined.off") << tm2;
}

//TODO: add more tests nm-edge on vertex, nm-edge vs nm-edge, ...


// coplanar face and non-manifold edge
{
  std::cout << "running corefinement test 4\n";
  Surface_mesh tm1;
  CGAL::make_quad(P(0,0,0), P(8,0,0), P(8,8,0), P(0,8,0), tm1);
  CGAL::make_quad(P(8,0,0), P(12,0,0), P(12,8,0), P(8,8,0), tm1);
  PMP::stitch_borders(tm1);
  CGAL::make_quad(P(8,0,0), P(8,8,0), P(8,8,4), P(8,0,4), tm1);
  CGAL::make_quad(P(8,0,0), P(8,8,0), P(8,8,-4), P(8,0,-4), tm1);
  PMP::triangulate_faces(tm1);

  Surface_mesh tm2;
  CGAL::make_quad(P(6,1,0), P(14,3,0), P(13,7,0), P(5,5,0), tm2);
  CGAL::make_quad(P(6,1,0), P(5,5,0), P(5,5,2), P(6,1,3), tm2);
  PMP::stitch_borders(tm2);
  PMP::triangulate_faces(tm2);

  std::ofstream("t4_tm1.off") << tm1;
  std::ofstream("t4_tm2.off") << tm2;

  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_1(tm1, get(boost::vertex_point, tm1));

  PMP::corefine(tm1, tm2, CGAL::parameters::non_manifold_feature_map(nm_map_1));

  std::ofstream("t4_tm1_corefined.off") << tm1;
  std::ofstream("t4_tm2_corefined.off") << tm2;
}
//
// // coplanar face and non-manifold edge and regular intersection with incident face
{
  std::cout << "running corefinement test 5\n";
  Surface_mesh tm1;
  CGAL::make_quad(P(0,0,0), P(4,0,0), P(4,8,0), P(0,8,0), tm1);
  CGAL::make_quad(P(4,0,0), P(12,0,0), P(12,8,0), P(4,8,0), tm1);
  PMP::stitch_borders(tm1);
  CGAL::make_quad(P(4,0,0), P(4,8,0), P(4,8,4), P(4,0,4), tm1);
  CGAL::make_quad(P(4,0,0), P(4,8,0), P(4,8,-4), P(4,0,-4), tm1);
  PMP::triangulate_faces(tm1);

  Surface_mesh tm2;
  CGAL::make_quad(P(2,1,0), P(14,3,0), P(13,7,0), P(5,5,0), tm2);
  CGAL::make_quad(P(2,1,0), P(5,5,0), P(5,5,2), P(2,1,3), tm2);
  PMP::stitch_borders(tm2);
  PMP::triangulate_faces(tm2);

  std::ofstream("t5_tm1.off") << tm1;
  std::ofstream("t5_tm2.off") << tm2;

  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_1(tm1, get(boost::vertex_point, tm1));

  PMP::corefine(tm1, tm2, CGAL::parameters::non_manifold_feature_map(nm_map_1));

  std::ofstream("t5_tm1_corefined.off") << tm1;
  std::ofstream("t5_tm2_corefined.off") << tm2;
}

// coplanar face and 2 non-manifold edges
{
  std::cout << "running corefinement test 6\n";
  Surface_mesh tm1;
  CGAL::make_quad(P(0,0,0), P(8,0,0), P(8,8,0), P(0,8,0), tm1);
  CGAL::make_quad(P(8,0,0), P(12,0,0), P(12,8,0), P(8,8,0), tm1);
  PMP::stitch_borders(tm1);
  CGAL::make_quad(P(8,0,0), P(8,8,0), P(8,8,4), P(8,0,4), tm1);
  CGAL::make_quad(P(8,0,0), P(8,8,0), P(8,8,-4), P(8,0,-4), tm1);
  PMP::triangulate_faces(tm1);

  Surface_mesh tm2;
  CGAL::make_quad(P(6,1,0), P(14,3,0), P(13,7,0), P(5,5,0), tm2);
  CGAL::make_quad(P(6,1,0), P(5,5,0), P(5,5,2), P(6,1,3), tm2);
  PMP::stitch_borders(tm2);
  CGAL::make_quad(P(6,1,0), P(5,5,0), P(5,5,-2), P(6,1,-3), tm2);
  CGAL::make_quad(P(6,1,0), P(5,5,0), P(3,5,0), P(3,1,0), tm2);
  PMP::triangulate_faces(tm2);

  std::ofstream("t6_tm1.off") << tm1;
  std::ofstream("t6_tm2.off") << tm2;

  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_1(tm1, get(boost::vertex_point, tm1));
  PMP::Non_manifold_feature_map<Surface_mesh> nm_map_2(tm2, get(boost::vertex_point, tm2));

#if 0
  std::vector< std::vector<P> > polylines;
  PMP::surface_intersection(tm1, tm2, std::back_inserter(polylines),
                                      CGAL::parameters::non_manifold_feature_map(nm_map_1),
                                      CGAL::parameters::non_manifold_feature_map(nm_map_2));

  //dump polylines
  std::ofstream output("intersection_polylines.cgal");
  for(const std::vector<P>& polyline : polylines)
  {
    output << polyline.size() << " ";
    std::copy(polyline.begin(), polyline.end(),std::ostream_iterator<P>(output," "));
    output << "\n";
  }
#else
  PMP::corefine(tm1, tm2, CGAL::parameters::non_manifold_feature_map(nm_map_1),
                          CGAL::parameters::non_manifold_feature_map(nm_map_2));

  std::ofstream("t6_tm1_corefined.off") << tm1;
  std::ofstream("t6_tm2_corefined.off") << tm2;
#endif
}


}
