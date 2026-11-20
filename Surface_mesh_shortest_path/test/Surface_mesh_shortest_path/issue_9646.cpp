#include <CGAL/Surface_mesh_shortest_path.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>

#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

#include <iostream>
#include <iomanip>
#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel   Kernel;
typedef CGAL::Surface_mesh<Kernel::Point_3>                   Mesh;

typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits>              SMSP;

typedef CGAL::AABB_face_graph_triangle_primitive<Mesh>        Primitive;
typedef CGAL::AABB_traits_3<Kernel, Primitive>                AABB_traits;
typedef CGAL::AABB_tree<AABB_traits>                          Tree;

typedef boost::graph_traits<Mesh>::vertex_descriptor          vertex_descriptor;

int main(int, char**)
{
  std::cout.precision(17);
  std::cerr.precision(17);

  for (int n : {4, 5, 6, 7, 8, 9})
  {
    // planar unit square, n x n vertices, two triangles per cell
    Mesh mesh;
    std::vector<Mesh::Vertex_index> vi;
    for (int i = 0; i < n; ++i)
      for (int j = 0; j < n; ++j)
        vi.push_back(mesh.add_vertex(Kernel::Point_3(double(i)/(n-1), double(j)/(n-1), 0.0)));
    for (int i = 0; i < n - 1; ++i)
      for (int j = 0; j < n - 1; ++j)
      {
        int a = i*n+j, b = i*n+j+1, c = (i+1)*n+j, d = (i+1)*n+j+1;
        mesh.add_face(vi[a], vi[b], vi[d]);
        mesh.add_face(vi[a], vi[d], vi[c]);
      }

    const vertex_descriptor vd(static_cast<std::size_t>(n - 1)); // corner (0, 1, 0)

    // ---
    SMSP by_vertex(mesh);
    by_vertex.add_source_point(vd);
    const double d_vertex = by_vertex.shortest_distance_to_source_points(vd).first;
    const auto loc_v = by_vertex.face_location(vd);

    // ---
//  SMSP by_point(mesh);
//  Tree tree;
//  by_point.build_aabb_tree(tree);
//  const auto loc_p = by_point.locate(mesh.point(vd), tree);
//  by_point.add_source_point(loc_p);
//  const double d_point = by_point.shortest_distance_to_source_points(vd).first;

    // ---
    SMSP by_snapped_point(mesh);
    const auto snapped_loc_p = SMSP::locate<AABB_traits>(mesh.point(vd), mesh, CGAL::parameters::snapping_tolerance(1e-10));
    by_snapped_point.add_source_point(snapped_loc_p);
    const double d_snapped_point = by_snapped_point.shortest_distance_to_source_points(vd).first;

    // ---
    std::cout << "\nn=" << n
              << "  d(vertex source) = " << d_vertex
//            << "  d(point source) = "  << d_point
              << "  d(snapped point source) = "  << d_snapped_point
              << "\n    loc_v::bary = "
              << loc_v.second[0] << ", " << loc_v.second[1] << ", " << loc_v.second[2]
//            << "\n    loc_p::bary = "
//            << loc_p.second[0] << ", " << loc_p.second[1] << ", " << loc_p.second[2]
              << "\n    snapped_loc_p::bary = "
              << snapped_loc_p.second[0] << ", " << snapped_loc_p.second[1] << ", " << snapped_loc_p.second[2]
              << "\n";

    assert(snapped_loc_p.second[2] == 1);
  }

  std::cout << "Done" << std::endl;
  return EXIT_SUCCESS;
}
