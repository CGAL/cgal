#include <CGAL/Surface_mesh.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Vector_graphics_on_surfaces/locally_shortest_path.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/boost/graph/IO/polygon_mesh_io.h>

#include <iostream>

namespace PMP = CGAL::Polygon_mesh_processing;

using K = CGAL::Exact_predicates_inexact_constructions_kernel;
using Mesh = CGAL::Surface_mesh<K::Point_3>;
using Face_location = PMP::Face_location<Mesh, double>;
using Edge_location = PMP::Edge_location<Mesh, double>;

int main()
{
  Mesh mesh;
  Mesh::Halfedge_index h = CGAL::make_quad(K::Point_3(0,0,0), K::Point_3(1,0,0),
                                           K::Point_3(1,1,0),K::Point_3(0,1,0),
                                           mesh);
  CGAL::Euler::split_face(h, next(next(h,mesh), mesh), mesh);
  PMP::isotropic_remeshing(faces(mesh), 0.15, mesh);

  std::ofstream("square.off") << std::setprecision(17) << mesh;

  // test starting on edges
  //int runid=0;
  for (Mesh::Halfedge_index h : halfedges(mesh))
  {
    if (is_border(h,mesh)) continue;

    Mesh::Face_index f = face(h, mesh);
    Face_location src(f, CGAL::make_array(0.5, 0.5, 0));

    // K::Point_3 src_pt = PMP::construct_point(src, mesh);
    // std::cout << "src = " << src_pt << "\n";

    double target_distance = 0.3;

    //std::ofstream out("straightest_geodesic_path_"+std::to_string(runid)+".polylines.txt");
    for (double n=0; n<8; ++n)
    {
      double theta = 2*CGAL_PI/8*n;
      K::Vector_2 dir(std::cos(theta), std::sin(theta));
      std::vector<Face_location> path = PMP::straightest_geodesic<K>(src, dir, target_distance, mesh);

      //TODO: check the output is of correct length

/*
      std::vector<K::Point_3> poly;
      poly.reserve(path.size());
      PMP::convert_path_to_polyline(path, mesh, std::back_inserter(poly));


      out << path.size() << " ";

      for (auto p : poly)
        out << " " << p;
      out << "\n";
*/
    }
//    ++runid;
  }

  // test starting on vertices
  //int runid=0;
  for (Mesh::Halfedge_index h : halfedges(mesh))
  {
    if (is_border(h,mesh)) continue;

    Mesh::Face_index f = face(h, mesh);
    Face_location src(f, CGAL::make_array(0., 1., 0.));

    // K::Point_3 src_pt = PMP::construct_point(src, mesh);
    // std::cout << "src = " << src_pt << "\n";

    double target_distance = 0.3;

    //std::ofstream out("straightest_geodesic_path_"+std::to_string(runid)+".polylines.txt");
    for (double n=0; n<8; ++n)
    {
      double theta = 2*CGAL_PI/16*n;
      K::Vector_2 dir(std::cos(theta), std::sin(theta));
      std::vector<Face_location> path = PMP::straightest_geodesic<K>(src, dir, target_distance, mesh);

      //TODO: check the output is of correct length

/*
      std::vector<K::Point_3> poly;
      poly.reserve(path.size());
      PMP::convert_path_to_polyline(path, mesh, std::back_inserter(poly));


      out << path.size() << " ";

      for (auto p : poly)
        out << " " << p;
      out << "\n";
*/
    }
//    ++runid;
  }

  //TODO: we still need to check that the output is correct (which is currently not the case at vertex)
  return 1;
}
