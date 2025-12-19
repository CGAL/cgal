#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>
#include <sstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = CGAL::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;

template <class TriangleMesh>
void test()
{
  // test with a clipper mesh
  {
    TriangleMesh tm1, tm2;
    std::ifstream(CGAL::data_file_path("meshes/elephant.off")) >> tm1;
    std::ifstream(CGAL::data_file_path("meshes/sphere.off")) >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false).face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(!CGAL::is_closed(tm1));
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream(CGAL::data_file_path("meshes/elephant.off")) >> tm1;
    std::ifstream(CGAL::data_file_path("meshes/sphere.off")) >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(true).face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(CGAL::is_closed(tm1));
  }
  //test with SI
  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-clip/clipper_for_tet.off") >> tm2;
    std::ifstream("data-clip/tet_with_si_to_clip.off") >> tm1;

    PMP::clip(tm1, tm2,
              params::throw_on_self_intersection(true),
              params::do_not_modify(true));
    std::vector<TriangleMesh> meshes;
    PMP::split_connected_components(tm1, meshes, params::default_values());
    assert(meshes.size() == 2);
    //if the order is not deterministc, put the num_vertices in a list and check
    //if the list does contain all those numbers.
    assert(num_vertices(meshes[0]) == 21);
    assert(num_vertices(meshes[1]) == 6);
  }

  // test with a iso-cuboid
  {
    TriangleMesh tm1;
    std::ifstream(CGAL::data_file_path("meshes/elephant.off")) >> tm1;
    K::Iso_cuboid_3 iso_cuboid(K::Point_3(0,0,0), K::Point_3(0.4, 0.6, 0.4));

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);

    PMP::clip(tm1, iso_cuboid, params::clip_volume(true).face_index_map(custom_face_index_map_1));
    assert(CGAL::is_closed(tm1));
  }

  // test with a plane
  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;
    K::Plane_3 plane(0, 0, 1, -1);

    PMP::clip(tm1, plane, params::clip_volume(true));
    assert(CGAL::is_closed(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;
    K::Plane_3 plane(0, 0, 1, -1);

    PMP::clip(tm1, plane, params::clip_volume(false).use_compact_clipper(false));
    assert(!CGAL::is_closed(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;
    K::Plane_3 plane(0, 0, 1, -1);

    PMP::clip(tm1, plane, params::clip_volume(false).use_compact_clipper(true));
    assert(CGAL::is_closed(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;

    PMP::clip(tm1, K::Plane_3(-0.236474, 0.437732, 0.867451, -0.838791), params::clip_volume(true));
    assert(CGAL::is_closed(tm1));
    assert(!CGAL::is_empty(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;

    PMP::clip(tm1, K::Plane_3(0, 0, 1, 2));
    assert(CGAL::is_empty(tm1));
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;

    PMP::clip(tm1, K::Plane_3(0, 0, 1, -2));
    assert(!CGAL::is_empty(tm1));
  }

  //test with SI
  {
    TriangleMesh tm1;
    std::ifstream("data-clip/tet_with_si_to_clip.off") >> tm1;
    if(num_vertices(tm1) == 0)
    {
      std::cerr<<"File not found. Aborting."<<std::endl;
      assert(false);
      return ;
    }
    PMP::clip(tm1, K::Plane_3(0,0,1,-4.5),
               params::throw_on_self_intersection(true)
               .allow_self_intersections(true));
    assert(vertices(tm1).size() == 13);
  }

  // clipping with identity
  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(true)
                     .use_compact_clipper(true)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(num_vertices(tm1) == 8);
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .use_compact_clipper(false)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(CGAL::is_empty(tm1));
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .use_compact_clipper(true)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(num_vertices(tm1) == 8);
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION, K::Vector_3(1,0,0)), tm2);
    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .use_compact_clipper(false)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(CGAL::is_empty(tm1));
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION, K::Vector_3(1,0,0)), tm2);
    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .use_compact_clipper(true)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 4);
  }

  // test orientation + patch without input vertex
  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;

    CGAL::make_tetrahedron(K::Point_3(0.53, -1.3, 0.2),
                           K::Point_3(0.53, 1.1, 0.2),
                           K::Point_3(0.53, -1.3, 0.4),
                           K::Point_3(0.73, -1.3, 0.2),
                           tm2);

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 6);
  }

  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-coref/cube.off") >> tm1;

    CGAL::make_tetrahedron(K::Point_3(0.53, -1.3, 0.2),
                           K::Point_3(0.53, 1.1, 0.2),
                           K::Point_3(0.53, -1.3, 0.4),
                           K::Point_3(0.73, -1.3, 0.2),
                           tm2);
    PMP::reverse_face_orientations(tm2);

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::clip_volume(false)
                     .face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 6+8);
  }

  // clip meshes with intersection polyline opened
  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(1, 0, 0, -2));
    assert(vertices(tm1).size() == 4);
  }

  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(-1, 0, 0, 2));
    assert(vertices(tm1).size() == 3);
  }

  // test with clipper on border edge
  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 1, 0), K::Point_3(1, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(0, 1, 0 , 0));
    assert(vertices(tm1).size() == 0);
  }

  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 1, 0), K::Point_3(1, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(0, -1, 0 , 0));
    assert(vertices(tm1).size() == 3);
  }

  // test with clipper on border edge: full triangle
  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(0, 0, 1, 0), params::use_compact_clipper(true));
    assert(vertices(tm1).size()!=0);
  }

  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(0, 0, 1, 0), params::use_compact_clipper(false));
    assert(vertices(tm1).size() == 0);
  }

  // test tangencies
  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 2, 0), K::Point_3(1, 1, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(1, 0, 0, -1));
    assert(vertices(tm1).size() == 3);
  }

  {
    TriangleMesh tm1;
    make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 2, 0), K::Point_3(1, 1, 0), tm1 );
    PMP::clip(tm1, K::Plane_3(-1, 0, 0, 1));
    assert(vertices(tm1).size() == 0);
  }

  {
    TriangleMesh tm1, tm2;
    make_triangle( K::Point_3(0.5, 0, 0.5), K::Point_3(1, 0.5, 0.5), K::Point_3(0.5, 1, 0.5), tm1 );
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 3);
  }

  {
    TriangleMesh tm1, tm2;
    make_triangle( K::Point_3(0.5, 0, 0.5), K::Point_3(1, 0.5, 0.5), K::Point_3(0.5, 1, 0.5), tm1 );
    std::ifstream("data-coref/cube.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::reverse_face_orientations(tm2);
    PMP::clip(tm1, tm2,
              params::face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(vertices(tm1).size() == 0);
  }

  // test combinations of use_compact_clipper and clip_volume
  {
    TriangleMesh tm1;
    std::ifstream("data-coref/cube.off") >> tm1;

    //  -> closed mesh, true/true
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(true));
    assert(faces(tm1).size() == 12);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, false/true
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(false).clip_volume(true));
    assert(faces(tm1).size() == 12);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, true/false
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(false));
    assert(faces(tm1).size() == 12);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, false/false
    PMP::clip(tm1, K::Plane_3(1,0,0,-1), params::use_compact_clipper(false).clip_volume(false));
    assert(faces(tm1).size() == 10);
    assert(!CGAL::is_closed(tm1));

    // -> open mesh true/true
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(true));
    assert(faces(tm1).size() == 10);

    // -> open mesh true/false
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(false));
    assert(faces(tm1).size() == 10);

    // -> open mesh false/false
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(false).clip_volume(false));
    assert(faces(tm1).size() == 8);

    // -> open mesh false/true
    PMP::clip(tm1, K::Plane_3(0,-1,0,0), params::use_compact_clipper(false).clip_volume(true));
    assert(faces(tm1).size() == 6);
  }

  // test special case
  {
    TriangleMesh tm1, tm2;
    std::ifstream("data-clip/tm_1.off") >> tm2;
    std::ifstream("data-clip/clipper_1.off") >> tm2;

    auto custom_face_index_map_1 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm1);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_1, tm1);
    auto custom_face_index_map_2 = get(CGAL::dynamic_face_property_t<std::size_t>(), tm2);
    CGAL::BGL::internal::initialize_face_index_map(custom_face_index_map_2, tm2);

    PMP::clip(tm1, tm2,
              params::face_index_map(custom_face_index_map_1),
              params::face_index_map(custom_face_index_map_2));
    assert(is_valid_polygon_mesh(tm1));
  }

  // non-manifold border vertices
  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 5 2 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 0 1 4\n3 1 2 3\n";
    ss >> tm1;
    PMP::clip(tm1, K::Plane_3(-1,0,0,2));
    assert(vertices(tm1).size() == 3);
  }

  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 7 4 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 3 5\n";
    ss >> tm1;
    CGAL::Euler::remove_face(halfedge(*std::prev(faces(tm1).end()),tm1),tm1);
    PMP::clip(tm1, K::Plane_3(-1,0,0,2));
    assert(vertices(tm1).size() == 5);
  }

  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
    ss >> tm1;
    for (int i=0;i<3;++i)
      CGAL::Euler::remove_face(halfedge(*std::prev(faces(tm1).end()),tm1),tm1);
    PMP::clip(tm1, K::Plane_3(-1,0,0,2));
    assert(vertices(tm1).size() == 7);
  }

  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
    ss >> tm1;
    for (int i=0;i<3;++i)
      CGAL::Euler::remove_face(halfedge(*std::prev(faces(tm1).end()),tm1),tm1);
    PMP::clip(tm1, K::Plane_3(0,1,0,0));
    assert(vertices(tm1).size() == 3);
  }

  {
    TriangleMesh tm1;
    std::stringstream ss;
    ss << "OFF\n 9 7 0\n 0 0 0\n2 0 0\n4 0 0\n4 1 0\n0 1 0\n3 1 0\n 1 1 0\n3 -1 0\n1 -1 0\n3 0 1 4\n3 1 2 3\n3 1 5 6\n3 1 8 7\n3 1 3 5\n3 1 6 4\n3 1 0 8\n";
    ss >> tm1;
    for (int i=0;i<3;++i)
      CGAL::Euler::remove_face(halfedge(*std::prev(faces(tm1).end()),tm1),tm1);
    PMP::clip(tm1, K::Plane_3(0,-1,0,0));
    assert(vertices(tm1).size() == 7);
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/open_large_cube.off") >> tm1;
    PMP::clip(tm1, K::Plane_3(0,0,1,-1), CGAL::parameters::use_compact_clipper(false));
    assert(vertices(tm1).size()==753);
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/open_large_cube.off") >> tm1;
    std::size_t nbv = vertices(tm1).size();
    PMP::clip(tm1, K::Plane_3(0,0,1,-1), CGAL::parameters::use_compact_clipper(true));
    assert(vertices(tm1).size()==nbv);
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/open_large_cube.off") >> tm1;
    PMP::clip(tm1, K::Plane_3(0,0,1,-1), CGAL::parameters::use_compact_clipper(false).allow_self_intersections(true));
    assert(vertices(tm1).size()==753);
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/open_large_cube.off") >> tm1;
    std::size_t nbv = vertices(tm1).size();
    PMP::clip(tm1, K::Plane_3(0,0,1,-1), CGAL::parameters::use_compact_clipper(true).allow_self_intersections(true));
    assert(vertices(tm1).size()==nbv);
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/open_large_cube.off") >> tm1;
    PMP::clip(tm1, K::Plane_3(0,0,-1,1), CGAL::parameters::use_compact_clipper(false));
    assert(vertices(tm1).size()==0);
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/open_large_cube.off") >> tm1;
    PMP::clip(tm1, K::Plane_3(0,0,-1,1), CGAL::parameters::use_compact_clipper(true));
    assert(vertices(tm1).size()==176);
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/open_large_cube.off") >> tm1;
    PMP::clip(tm1, K::Plane_3(0,0,-1,1), CGAL::parameters::use_compact_clipper(false).allow_self_intersections(true));
    assert(vertices(tm1).size()==0);
  }

  {
    TriangleMesh tm1;
    std::ifstream("data-coref/open_large_cube.off") >> tm1;
    PMP::clip(tm1, K::Plane_3(0,0,-1,1), CGAL::parameters::use_compact_clipper(true).allow_self_intersections(true));
    assert(vertices(tm1).size()==176);
  }

  // non-simply connected output faces
  {
    TriangleMesh tm1;

    CGAL::make_hexahedron(CGAL::Bbox_3(-3.5,-0.5,-0.5, -2.5,0.5,0.5), tm1);
    PMP::reverse_face_orientations(tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-4,-1,-1, -2,1,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-1,-1,-1, 1,1,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(2,-1,-1, 4,1,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(5,-1,-1, 7,1,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-1,2,-1, 1,4,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(2,2,-1, 4,4,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(5,2,-1, 7,4,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-1,5,-1, 1,7,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(2,5,-1, 4,7,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(5,5,-1, 7,7,1), tm1);

    PMP::clip(tm1, K::Plane_3(0,0,1,0), CGAL::parameters::do_not_triangulate_faces(true).clip_volume(true));

    assert(vertices(tm1).size()==88);
    assert(faces(tm1).size()==72);
  }
  {
    TriangleMesh tm1;

    CGAL::make_hexahedron(CGAL::Bbox_3(-1,-1,-1, 1,1,1), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-3,-3,-3, 3,3,3), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-5,-5,-5, 5,5,5), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-7,-7,-7, 7,7,7), tm1);
    PMP::reverse_face_orientations(tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-0.5,-0.5,-0.5, 0.5,0.5,0.5), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-2,-2,-2, 2,2,2), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-4,-4,-4, 4,4,4), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-6,-6,-6, 6,6,6), tm1);
    CGAL::make_hexahedron(CGAL::Bbox_3(-8,-8,-8, 8,8,8), tm1);

    PMP::clip(tm1, K::Plane_3(0,0,1,0), CGAL::parameters::do_not_triangulate_faces(true).clip_volume(true));
    assert(vertices(tm1).size()==72);
    assert(faces(tm1).size()==79);
  }


}

template <class Mesh>
void test_split_plane()
{
//test with a splitter mesh
  Mesh tm1;
  std::ifstream input(CGAL::data_file_path("meshes/elephant.off"));
  input >> tm1;

  if(!input)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }
  input.close();

  PMP::split(tm1,K::Plane_3(0,0,1,0));

  std::vector<Mesh> meshes;
  PMP::split_connected_components(tm1, meshes, params::default_values());
  assert(meshes.size() == 3);
  //if the order is not deterministc, put the num_vertices in a list and check
  //if the list does contain all those numbers.

  assert(num_vertices(meshes[2]) == 46);
  assert(num_vertices(meshes[0]) == 1523);
  assert(num_vertices(meshes[1]) == 1668);

  CGAL::clear(tm1);
  meshes.clear();

//test with a non-closed splitter mesh (border edges in the plane)
  input.open("data-coref/open_large_cube.off");
  input >> tm1;

  if(!input)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }
  input.close();

  PMP::split(tm1,K::Plane_3(0,0,1,-1));
  PMP::split_connected_components(tm1, meshes, params::default_values());
  assert(meshes.size() == 2);

  CGAL::clear(tm1);
  meshes.clear();

//test with a non-closed splitter mesh (border edges in the plane)
  input.open("data-coref/open_large_cube.off");
  input >> tm1;

  if(!input)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }
  input.close();

  PMP::split(tm1,K::Plane_3(0,-1,0,0.3));
  PMP::split_connected_components(tm1, meshes, params::default_values());
  assert(meshes.size() == 2);

  CGAL::clear(tm1);
  meshes.clear();

//test with SI
  std::ifstream("data-clip/tet_si_to_split.off") >> tm1;
  if(num_vertices(tm1) == 0)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }

  PMP::split(tm1, K::Plane_3(0,0,1,-0.5),
             params::throw_on_self_intersection(true)
             .allow_self_intersections(true));
  PMP::split_connected_components(tm1, meshes, params::default_values());
  assert(meshes.size() == 2);
  //if the order is not deterministc, put the num_vertices in a list and check
  //if the list does contain all those numbers.
  assert(num_vertices(meshes[0]) == 12);
  assert(num_vertices(meshes[1]) == 12);

  CGAL::clear(tm1);
  meshes.clear();
}

template <class TriangleMesh>
void test_split()
{
  // test with a clipper mesh
  TriangleMesh tm1, tm2;
  //closed intersection curves
  std::ifstream input(CGAL::data_file_path("meshes/elephant.off"));
  input >> tm1;

  if(!input)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }

  input.close();

  input.open(CGAL::data_file_path("meshes/sphere.off"));
  input >> tm2;

  if(!input)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }

  input.close();

  std::vector<TriangleMesh> meshes;

  PMP::split(tm1, tm2);
  //try with np
  typedef typename boost::graph_traits<TriangleMesh>::faces_size_type  faces_size_type;
  typename boost::template property_map<
  TriangleMesh, CGAL::dynamic_face_property_t<faces_size_type> >::type
      pidmap = get(CGAL::dynamic_face_property_t<faces_size_type>(), tm1);
  CGAL::Polygon_mesh_processing::connected_components(
        tm1, pidmap, CGAL::parameters::default_values());
  PMP::split_connected_components(tm1,
                                  meshes,
                                  params::face_patch_map(pidmap));

  assert(meshes.size() == 5);
  //if the order is not deterministc, put the num_vertices in a list and check
  //if the list does contain all those numbers.
  assert(num_vertices(meshes[0]) == 2641);
  assert(num_vertices(meshes[1]) == 159);
  assert(num_vertices(meshes[2]) == 142);
  assert(num_vertices(meshes[3]) == 83);
  assert(num_vertices(meshes[4]) == 104);
  assert(tm1.is_valid());

  CGAL::clear(tm1);
  CGAL::clear(tm2);
  meshes.clear();

  // open intersection curve
  input.open("data-clip/split_A.off");
  input >> tm1;

  if(!input)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }

  input.close();

  input.open("data-clip/splitter.off");
  input >> tm2;

  if(!input)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }

  input.close();

  PMP::split(tm1, tm2);
  PMP::split_connected_components(tm1,
                                  meshes,
                                  params::default_values());

  assert(meshes.size() == 2);
  //if the order is not deterministc, put the num_vertices in a list and check
  //if the list does contain all those numbers.
  assert(num_vertices(meshes[0]) == 588);
  assert(num_vertices(meshes[1]) == 50);
  assert(tm1.is_valid());

  CGAL::clear(tm1);
  CGAL::clear(tm2);
  meshes.clear();


  //test with SI
  std::ifstream("data-clip/splitter_for_tet.off") >> tm2;
  std::ifstream("data-clip/tet_si_to_split.off") >> tm1;
  //test with SI
  if(num_vertices(tm1) == 0 || num_vertices(tm2) == 0)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }
  PMP::split(tm1, tm2,
             params::throw_on_self_intersection(true),
             params::do_not_modify(true));
  PMP::split_connected_components(tm1, meshes, params::default_values());
  assert(meshes.size() == 3);
  //if the order is not deterministc, put the num_vertices in a list and check
  //if the list does contain all those numbers.
  assert(num_vertices(meshes[0]) == 29);
  assert(num_vertices(meshes[1]) == 8);
  assert(num_vertices(meshes[2]) == 17);

  CGAL::clear(tm1);
  CGAL::clear(tm2);
  meshes.clear();
}

template <class TriangleMesh>
void test_isocuboid()
{
  TriangleMesh tm;
  //closed intersection curves
  std::ifstream input(CGAL::data_file_path("meshes/elephant.off"));
  input >> tm;

  if(!input)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }

  input.close();

  std::vector<TriangleMesh> meshes;
  K::Iso_cuboid_3 splitter(K::Point_3(-0.3, -0.45, -0.25),
                           K::Point_3( 0.3, 0.45, 0.25));
  PMP::split(tm, splitter);

  PMP::split_connected_components(tm,
                                  meshes);

  assert(meshes.size() == 10);
  //if the order is not deterministc, put the num_vertices in a list and check
  //if the list does contain all those numbers.
  assert(num_vertices(meshes[0]) == 2663);
  assert(num_vertices(meshes[1]) == 131 );
  assert(num_vertices(meshes[2]) == 32  );
  assert(num_vertices(meshes[3]) == 125 );
  assert(num_vertices(meshes[4]) == 224 );
  assert(num_vertices(meshes[5]) == 107 );
  assert(num_vertices(meshes[6]) == 121 );
  assert(num_vertices(meshes[7]) == 56  );
  assert(num_vertices(meshes[8]) == 49  );
  assert(num_vertices(meshes[9]) == 13  );
  assert(tm.is_valid());
  CGAL::clear(tm);
  meshes.clear();


  std::ifstream("data-clip/tet_with_si_to_clip.off") >> tm;
  if(num_vertices(tm) == 0)
  {
    std::cerr<<"File not found. Aborting."<<std::endl;
    assert(false);
    return ;
  }
  splitter = K::Iso_cuboid_3(K::Point_3(-2, 7, 4),
                           K::Point_3(-1, 8, 5));
  PMP::split(tm, splitter,
             params::throw_on_self_intersection(true)
             .allow_self_intersections(true));
  PMP::split_connected_components(tm, meshes, params::default_values());
  assert(meshes.size() == 4);

  std::set<std::size_t> sizes;
  for (int i=0; i<4; ++i)
    sizes.insert(vertices(meshes[i]).size());

  assert(sizes.count(20)==1);
  assert(sizes.count(21)==1);
  assert(sizes.count(7)==1);
  assert(sizes.count(4)==1);

  CGAL::clear(tm);
  meshes.clear();

  std::ifstream("data-clip/tet_with_si_to_clip.off") >> tm;
  PMP::clip(tm, splitter,
             params::throw_on_self_intersection(true)
             .allow_self_intersections(true));
  PMP::split_connected_components(tm, meshes, params::default_values());
  assert(meshes.size() == 2);
  //if the order is not deterministc, put the num_vertices in a list and check
  //if the list does contain all those numbers.
  assert(vertices(meshes[0]).size() == 20);
  assert(vertices(meshes[1]).size() == 4);
}

template <class TriangleMesh>
void test_new_clip()
{
  {
    TriangleMesh e;
    std::ifstream("data-clip/ee.off") >> e;
    PMP::refine_with_plane(e, K::Plane_3(1,0,0,-2));
    assert(faces(e).size()==5);
    assert(vertices(e).size()==24);
  }

  {
    TriangleMesh c;
    std::ifstream("data-clip/c.off") >> c;
    PMP::refine_with_plane(c, K::Plane_3(1,0,0,-2));
    assert(faces(c).size()==2);
    assert(vertices(c).size()==8);
  }

  {
    TriangleMesh e;
    std::ifstream("data-clip/ee.off") >> e;
    PMP::triangulate_faces(e);
    PMP::refine_with_plane(e, K::Plane_3(1,0,0,-2), CGAL::parameters::do_not_triangulate_faces(false));
    assert(faces(e).size()==30);
    assert(vertices(e).size()==28);
  }

  {
    TriangleMesh c;
    std::ifstream("data-clip/c.off") >> c;
    PMP::triangulate_faces(c);
    PMP::refine_with_plane(c, K::Plane_3(1,0,0,-2), CGAL::parameters::do_not_triangulate_faces(false));
    assert(faces(c).size()==8);
    assert(vertices(c).size()==9);
  }

  {
    TriangleMesh ele;
    std::ifstream(CGAL::data_file_path("meshes/elephant.off")) >> ele;
    PMP::clip(ele, K::Plane_3(1,0,0,0), CGAL::parameters::do_not_triangulate_faces(true).clip_volume(true));
    PMP::clip(ele, K::Plane_3(0,1,0,0), CGAL::parameters::do_not_triangulate_faces(true).clip_volume(true));
    PMP::clip(ele, K::Plane_3(0,0,1,0), CGAL::parameters::do_not_triangulate_faces(true).clip_volume(true));
    assert(faces(ele).size()==1220);
    assert(vertices(ele).size()==691);
  }
}

struct Clip_and_split_visitor
{
  using halfedge_descriptor = typename boost::graph_traits<Surface_mesh>::halfedge_descriptor;
  using face_descriptor = typename boost::graph_traits<Surface_mesh>::face_descriptor;
  using vertex_descriptor = typename boost::graph_traits<Surface_mesh>::vertex_descriptor;

  Surface_mesh& sm;
  Surface_mesh::Property_map<Surface_mesh::Face_index, int> fid_map;
  Surface_mesh::Property_map<Surface_mesh::Vertex_index, int> vid_map;
  Surface_mesh::Property_map<Surface_mesh::Halfedge_index, int> hid_map;

  int fid=-1;
  int hid=-1;
  std::size_t nbf=0;
  std::size_t nbe=0;
  std::size_t nb_input_v=sm.number_of_vertices();
  std::size_t nb_input_v_on=0;
  std::size_t nb_new_v_on=0;

  Clip_and_split_visitor(Surface_mesh& sm)
    : sm(sm)
  {
    bool is_new=false;
    std::tie(fid_map, is_new)=sm.add_property_map<Surface_mesh::Face_index,int>("f:id", -1);
    assert(is_new);
    int i=0; for (auto f : faces(sm)) put(fid_map, f, i++);
    std::tie(vid_map, is_new)=sm.add_property_map<Surface_mesh::Vertex_index,int>("v:id", -1);
    assert(is_new);
    i=0; for (auto v : vertices(sm)) put(vid_map, v, i++);
    std::tie(hid_map, is_new)=sm.add_property_map<Surface_mesh::Halfedge_index,int>("h:id", -1);
    assert(is_new);
    i=0;
    for (auto e : edges(sm))
    {
      auto h = halfedge(e, sm);
      put(hid_map, opposite(h, sm), i);
      put(hid_map, h, i++);
    }
  }

  void before_subface_creations(face_descriptor f_split, const Surface_mesh& pm)
  {
    fid=get(fid_map, f_split);
    assert(fid!=-1);
    assert(&pm==&sm);
  }
  void after_subface_creations(const Surface_mesh& pm)
  {
    assert(&pm==&sm);
  }
  void before_subface_created(const Surface_mesh& pm)
  {
    nbf=sm.number_of_faces();
    assert(&pm==&sm);
  }
  void after_subface_created(face_descriptor f_new, const Surface_mesh& pm)
  {
    assert(&pm==&sm);
    assert(get(fid_map, f_new)==-1);
    put(fid_map, f_new, fid);
    assert(nbf+1==sm.number_of_faces());
  }

  void before_edge_split(halfedge_descriptor h, Surface_mesh& pm)
  {
    hid=get(hid_map, h);
    assert(hid!=-1);
    assert(&pm==&sm);
    nbe=sm.number_of_edges();
  }
  void edge_split(halfedge_descriptor hnew, Surface_mesh& pm)
  {
    assert(&pm==&sm);
    assert(get(hid_map, hnew)==-1);
    assert(hid!=-1);
    put(hid_map, hnew, hid);
    put(hid_map, opposite(hnew, sm), hid);
    assert(nbe+1==sm.number_of_edges());
  }
  void after_edge_split()
  {
    assert(nbe+1==sm.number_of_edges());
  }
  void add_retriangulation_edge(halfedge_descriptor hnew, const Surface_mesh& pm)
  {
    assert(get(hid_map, hnew)==-1);
    put(hid_map, hnew, -2);
    put(hid_map, opposite(hnew, sm), -2);
    assert(&pm==&sm);
  }

  void intersection_point_detected(std::size_t /* i_id */,
                                   int /* sdim */,
                                   halfedge_descriptor h_e,
                                   halfedge_descriptor h_f,
                                   const Surface_mesh& tm_e,
                                   const Surface_mesh& tm_f,
                                   bool is_target_coplanar,
                                   bool is_source_coplanar)
  {
    assert(is_source_coplanar==false);
    assert(h_f==boost::graph_traits<Surface_mesh>::null_halfedge());
    assert(h_e!=boost::graph_traits<Surface_mesh>::null_halfedge());
    if (!is_target_coplanar)
      ++nb_new_v_on;
    else
      ++nb_input_v_on;
    assert(&tm_e==&sm);
    assert(&tm_f==&sm);
  }

  void new_vertex_added(std::size_t id, vertex_descriptor v, const Surface_mesh& pm)
  {
    assert(&pm==&sm);
    assert(get(vid_map, v)==-1);
    put(vid_map, v, nb_input_v+id-nb_input_v_on);
  }

  void before_face_copy(face_descriptor f, const Surface_mesh& src, const Surface_mesh& tgt)
  {
    assert(f==boost::graph_traits<Surface_mesh>::null_face());
    assert(&src==&sm);
    assert(&tgt==&sm);
  }

  void after_face_copy(face_descriptor fsrc, const Surface_mesh& src, face_descriptor ftgt, const Surface_mesh& tgt)
  {
    assert(fsrc==boost::graph_traits<Surface_mesh>::null_face());
    assert(ftgt!=boost::graph_traits<Surface_mesh>::null_face());
    //assert(get(fid_map, ftgt)==-1);
    put(fid_map, ftgt, -2);
    assert(&src==&sm);
    assert(&tgt==&sm);
    for (auto h : halfedges_around_face(halfedge(ftgt, sm), sm))
    {
      if (get(hid_map, h)==-1)
        put(hid_map, h, -2);
    }
  }

  void before_edge_duplicated(halfedge_descriptor h, Surface_mesh& tm)
  {
    hid = get(hid_map, h);
    assert(hid!=-1);
    assert(&tm==&sm);
  }

  void after_edge_duplicated(halfedge_descriptor h, halfedge_descriptor new_hedge, Surface_mesh& tm)
  {
    assert(hid==get(hid_map, h));
    assert(&tm==&sm);
    put(hid_map, new_hedge, hid);
    put(hid_map, opposite(new_hedge, sm), hid);
  }

  void before_vertex_copy(vertex_descriptor v, Surface_mesh& src, Surface_mesh& tgt)
  {
    assert(&src==&sm);
    assert(&tgt==&sm);
    assert(get(vid_map, v)!=-1);
  }
  void after_vertex_copy(vertex_descriptor v, Surface_mesh& src, vertex_descriptor nv, Surface_mesh& tgt)
  {
    assert(&src==&sm);
    assert(&tgt==&sm);
    assert(get(vid_map, v)!=-1);
    put(vid_map, nv, get(vid_map, v));
  }



  void check()
  {
    for (auto f :faces(sm))
    {
      if (get(fid_map, f)==-1) std::cout << sm.point(source(halfedge(f, sm), sm)) << " " << sm.point(target(halfedge(f, sm), sm)) << " " << sm.point(target(next(halfedge(f, sm), sm), sm)) << "\n";
      assert(get(fid_map, f)!=-1);
    }
    for (auto h :halfedges(sm))
    {
      if (get(hid_map, h)==-1)
        std::cout << sm.point(source(h, sm)) << " " << sm.point(target(h,sm)) << "\n";
      assert(get(hid_map, h)!=-1);
    }
    std::size_t nbv_max=nb_input_v+nb_new_v_on;
    for (auto v :vertices(sm))
      assert(get(vid_map, v)!=-1 && get(vid_map, v)<(int)nbv_max);
  }

};

void test_clip_and_split_with_plane_visitor()
{
  auto test_clip =[](std::string fname, const K::Plane_3& plane, bool triangulate, bool clip_volume)
  {
    std::cout << "   testing with clip " << fname << " vs " << plane << " (" << triangulate << "," << clip_volume << ")\n";
    Surface_mesh sm;
    std::ifstream(fname) >> sm;
    Clip_and_split_visitor visitor(sm);
    visitor.check();
    PMP::clip(sm, plane, params::visitor(std::ref(visitor)).do_not_triangulate_faces(!triangulate).clip_volume(clip_volume));
    std::ofstream ("/tmp/out.off") << sm;
    visitor.check();
  };

  auto test_split =[](std::string fname, const K::Plane_3& plane, bool triangulate)
  {
    std::cout << "   testing with split" << fname << " vs " << plane << " (" << triangulate << ")\n";
    Surface_mesh sm;
    std::ifstream(fname) >> sm;
    Clip_and_split_visitor visitor(sm);
    visitor.check();
    PMP::split(sm, plane, params::visitor(std::ref(visitor)).do_not_triangulate_faces(!triangulate));
    visitor.check();
  };

  auto test = [&](std::string fname, const K::Plane_3& plane)
  {
    test_clip(fname, plane, true, true);
    test_clip(fname, plane, true, false);
    test_clip(fname, plane, false, false);
    test_clip(fname, plane, false, true);
    test_split(fname, plane, false);
    test_split(fname, plane, true);
  };

  test(CGAL::data_file_path("meshes/torus_quad.off"), K::Plane_3(0,0,1,0));
  test(CGAL::data_file_path("meshes/elephant.off"), K::Plane_3(0.137304, -0.293668, 0.945995, 0));

}

int main()
{
  std::cout << "Surface Mesh" << std::endl;
  test<Surface_mesh>();

  std::cout << "Polyhedron" << std::endl;
  test<Polyhedron>();

  std::cout << "running test_split with Surface_mesh\n";
  test_split<Surface_mesh>();

  std::cout << "running test_iso_cuboid with Surface_mesh\n";
  test_isocuboid<Surface_mesh>();
  std::cout << "running test_split_plane with Surface_mesh\n";
  test_split_plane<Surface_mesh>();
  std::cout << "running test_split with Polyhedron\n";
  test_split<Polyhedron>();
  std::cout << "running test_split_plane with Polyhedron\n";
  test_split_plane<Polyhedron>();
  std::cout << "running test_iso_cuboid with Polyhedron\n";
  test_isocuboid<Polyhedron>();
  std::cout << "Done!" << std::endl;
  std::cout << "running test_new_clip with Surface_mesh\n";
  test_new_clip<Surface_mesh>();
  std::cout << "Done!" << std::endl;
  std::cout << "running test_new_clip with Polyhedron\n";
  test_new_clip<Polyhedron>();
  std::cout << "Done!" << std::endl;
  std::cout << "running test_clip_and_split_with_plane_visitor\n";
  test_clip_and_split_with_plane_visitor();
  std::cout << "Done!" << std::endl;
  return EXIT_SUCCESS;
}
