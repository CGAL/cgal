#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_shortest_path.h>

#include <CGAL/AABB_face_graph_triangle_primitive.h>
#include <CGAL/AABB_traits_3.h>
#include <CGAL/AABB_tree.h>

#include <iostream>
#include <fstream>
#include <limits>
#include <stdlib.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;

typedef CGAL::Surface_mesh<Kernel::Point_3> Triangle_mesh;

typedef CGAL::Surface_mesh_shortest_path_traits<Kernel, Triangle_mesh> Traits;
typedef CGAL::Surface_mesh_shortest_path<Traits> Surface_mesh_shortest_path;

typedef boost::property_map<Triangle_mesh, CGAL::vertex_point_t>::const_type VPM;
typedef CGAL::AABB_face_graph_triangle_primitive<Triangle_mesh, VPM> AABB_face_graph_primitive;
typedef CGAL::AABB_traits_3<Kernel, AABB_face_graph_primitive> AABB_face_graph_traits;

void test_all_pairs()
{
  CGAL::Surface_mesh<Kernel::Point_3> mesh;
  std::ifstream input("data/test_mesh_6.off");
  input >> mesh;
  input.close();

  std::cout << "Input mesh: " << num_vertices(mesh) << " nv" << std::endl;

  for (Triangle_mesh::Vertex_index v1 : vertices(mesh))
  {
    Surface_mesh_shortest_path shortest_paths(mesh);
    shortest_paths.add_source_point(v1);

    for (Triangle_mesh::Vertex_index v2 : vertices(mesh))
    {
      const double sq_dist = CGAL::square(shortest_paths.shortest_distance_to_source_points(v2).first);
      const double lower_bound = CGAL::squared_distance(mesh.point(v1), mesh.point(v2));
      std::cout << "sq_dist(" << v1 << ", " << v2 << ") = " << sq_dist << std::endl;
      std::cout << "lower bound: " << lower_bound << std::endl;

      if(v1 == v2)
        assert(sq_dist == 0.);
      else
        assert(sq_dist > (1. - 1e-7) * lower_bound); // numerical errors

      CGAL_USE(sq_dist);
    }
  }
}

void test_flat_donut()
{
  CGAL::Surface_mesh<Kernel::Point_3> mesh;
  std::ifstream input("data/flat_donut.off");
  input >> mesh;
  input.close();

  std::cout << "Input mesh: " << num_vertices(mesh) << " nv" << std::endl;

  {
    Triangle_mesh::Face_index f0(0), f16(16);
    auto h0 = halfedge(f0, mesh), h16 = halfedge(f16, mesh);

    for(std::size_t i=0; i<3; ++i)
    {
      Triangle_mesh::Vertex_index v0(0);
      Surface_mesh_shortest_path shortest_paths(mesh);
      shortest_paths.add_source_point(v0);

      double dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(3)).first;
      assert(dist == 3);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(12)).first;
      assert(dist == 3);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(1)).first;
      assert(dist == 1);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(4)).first;
      assert(dist == 1);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(5)).first;
      assert(CGAL::abs(dist - sqrt(2.)) < 1e-7);

      // change the canonical halfedges to test barycentric coordinates
      h0 = next(h0, mesh);
      set_halfedge(f0, h0, mesh);
      h16 = next(h16, mesh);
      set_halfedge(f16, h16, mesh);
    }
  }

  {
    Triangle_mesh::Vertex_index v0(3);

    Surface_mesh_shortest_path shortest_paths(mesh);
    shortest_paths.add_source_point(v0);

    double dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(0)).first;
    assert(dist == 3);

    dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(15)).first;
    assert(dist == 3);

    dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(2)).first;
    assert(dist == 1);

    dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(7)).first;
    assert(dist == 1);

    dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(6)).first;
    assert(CGAL::abs(dist - sqrt(2.)) < 1e-7);
  }

  {
    Surface_mesh_shortest_path shortest_paths(mesh);
    shortest_paths.add_source_point(Triangle_mesh::Vertex_index(0));
    shortest_paths.add_source_point(Triangle_mesh::Vertex_index(3));

    double dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(0)).first;
    assert(dist == 0);

    dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(15)).first;
    assert(dist == 3);

    dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(2)).first;
    assert(dist == 1);

    dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(7)).first;
    assert(dist == 1);

    dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(6)).first;
    assert(CGAL::abs(dist - sqrt(2.)) < 1e-7);
  }

  {
    Triangle_mesh::Face_index f16(16), f1(1), f29(29), f6(6);
    auto h16 = halfedge(f16, mesh), h1 = halfedge(f1, mesh), h29 = halfedge(f29, mesh), h6 = halfedge(f6, mesh);

    for(std::size_t i=0; i<3; ++i)
    {
      Surface_mesh_shortest_path shortest_paths(mesh);
      shortest_paths.add_source_point(Triangle_mesh::Vertex_index(1));

      double dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(0)).first;
      assert(dist == 1);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(1)).first;
      assert(dist == 0);

      auto loc = shortest_paths.locate<AABB_face_graph_traits>(Kernel::Point_3(0.5, 0, 0));
      dist = shortest_paths.shortest_distance_to_source_points(loc.first, loc.second).first;
      assert(dist == 0.5);

      loc = shortest_paths.locate<AABB_face_graph_traits>(Kernel::Point_3(2.5, 0, 0));
      dist = shortest_paths.shortest_distance_to_source_points(loc.first, loc.second).first;
      assert(dist == 1.5);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(5)).first;
      assert(dist == 1);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(9)).first;
      assert(dist == 2);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(13)).first;
      assert(dist == 3);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(2)).first;
      assert(dist == 1);

      h16 = next(h16, mesh);
      set_halfedge(f16, h16, mesh);
      h1 = next(h1, mesh);
      set_halfedge(f1, h1, mesh);
      h29 = next(h29, mesh);
      set_halfedge(f29, h29, mesh);
      h6 = next(h6, mesh);
      set_halfedge(f6, h6, mesh);
    }
  }

  {
    Triangle_mesh::Face_index f6(6);
    auto h6 = halfedge(f6, mesh);

    for(std::size_t i=0; i<3; ++i)
    {
      Surface_mesh_shortest_path shortest_paths(mesh);
      shortest_paths.add_source_point(shortest_paths.locate<AABB_face_graph_traits>(Kernel::Point_3(1.5, 0, 0)));

      double dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(0)).first;
      assert(dist == 1.5);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(3)).first;
      assert(dist == 1.5);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(1)).first;
      assert(dist == 0.5);

      dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(26)).first;
      assert((dist - sqrt(2.)) < 1e-7);

      auto loc = shortest_paths.locate<AABB_face_graph_traits>(Kernel::Point_3(1.5, 1, 0));
      dist = shortest_paths.shortest_distance_to_source_points(loc.first, loc.second).first;
      assert(dist == 1);

      loc = shortest_paths.locate<AABB_face_graph_traits>(Kernel::Point_3(0.5, 1, 0));
      dist = shortest_paths.shortest_distance_to_source_points(loc.first, loc.second).first;
      assert(CGAL::abs(dist - sqrt(2.)) < 1e-7);

      h6 = next(h6, mesh);
      set_halfedge(f6, h6, mesh);
    }
  }

  {
    Triangle_mesh::Face_index f28(28);
    auto h28 = halfedge(f28, mesh);
    Triangle_mesh::Face_index f29(29);
    auto h29 = halfedge(f29, mesh);

    auto h = halfedge(Triangle_mesh::Vertex_index(31), Triangle_mesh::Vertex_index(25), mesh).first;

    for(std::size_t side=0; side<2; ++side)
    {
      for(std::size_t i=0; i<3; ++i)
      {
        Surface_mesh_shortest_path shortest_paths(mesh);
        shortest_paths.add_source_point(shortest_paths.face_location(h, 0.5));

        double dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(16)).first;
        assert(dist == 0.75);

        dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(31)).first;
        assert(CGAL::abs(dist - 0.25) < 1e-7); // numerical issue (returns '0.25000000000000006')

        dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(25)).first;
        assert(dist == 0.25);

        dist = shortest_paths.shortest_distance_to_source_points(Triangle_mesh::Vertex_index(23)).first;
        assert(dist == 1.25);

        auto loc = shortest_paths.locate<AABB_face_graph_traits>(Kernel::Point_3(1.25, 0, 0));
        dist = shortest_paths.shortest_distance_to_source_points(loc.first, loc.second).first;
        assert(dist == 0.5);

        loc = shortest_paths.locate<AABB_face_graph_traits>(Kernel::Point_3(1.25, 1, 0));
        dist = shortest_paths.shortest_distance_to_source_points(loc.first, loc.second).first;
        assert(dist == 0.5);

        h28 = next(h28, mesh);
        set_halfedge(f28, h28, mesh);
        h29 = next(h29, mesh);
        set_halfedge(f29, h29, mesh);
      }

      h = opposite(h, mesh);
    }
  }
}

int main(int, char**)
{
  std::cerr.precision(17);
  std::cout.precision(17);

//  test_all_pairs();
  test_flat_donut();

  std::cout << "Done!" << std::endl;

  return EXIT_SUCCESS;
}
