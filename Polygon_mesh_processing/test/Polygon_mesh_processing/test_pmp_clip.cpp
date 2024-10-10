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
    assert(vertices(tm1).size() == 16);
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
    assert(vertices(tm1).size() == 4);
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
    assert(faces(tm1).size() == 14);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, false/true
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(false).clip_volume(true));
    assert(faces(tm1).size() == 14);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, true/false
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(false));
    assert(faces(tm1).size() == 14);
    assert(CGAL::is_closed(tm1));

    //  -> closed mesh, false/false
    PMP::clip(tm1, K::Plane_3(1,0,0,-1), params::use_compact_clipper(false).clip_volume(false));
    assert(faces(tm1).size() == 12);
    assert(!CGAL::is_closed(tm1));

    // -> open mesh true/true
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(true));
    assert(faces(tm1).size() == 12);

    // -> open mesh true/false
    PMP::clip(tm1, K::Plane_3(-1,0,0,0), params::use_compact_clipper(true).clip_volume(false));
    assert(faces(tm1).size() == 12);

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
    assert(vertices(tm1).size() == 6);
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
    assert(vertices(tm1).size()==nbv+2); // because of the plane diagonal
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
    assert(vertices(tm1).size()==nbv+2); // because of the plane diagonal
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
    assert(vertices(tm1).size()==178);
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
    assert(vertices(tm1).size()==178);
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
  assert(num_vertices(meshes[2]) == 48);
  assert(num_vertices(meshes[0]) == 1527);
  assert(num_vertices(meshes[1]) == 1674);

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
  assert(meshes.size() == 281);

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
  assert(num_vertices(meshes[0]) == 16);
  assert(num_vertices(meshes[1]) == 16);

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
  return EXIT_SUCCESS;
}
