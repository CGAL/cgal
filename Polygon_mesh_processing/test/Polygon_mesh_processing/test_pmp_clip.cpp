#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polygon_mesh_processing/clip.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polygon_mesh_processing/transform.h>

#include <CGAL/boost/graph/Face_filtered_graph.h>
#include <boost/property_map/property_map.hpp>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Surface_mesh;
typedef CGAL::Polyhedron_3<K> Polyhedron;

template <class TriangleMesh>
void test()
{
  // test with a clipper mesh
  TriangleMesh tm1, tm2;

  std::ifstream input("data-coref/elephant.off");
  input >> tm1;
  input.close();
  input.open("data-coref/sphere.off");
  input >> tm2;
  input.close();

  PMP::clip(tm1, tm2,
            params::clip_volume(false)
              .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
            params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2))
  );
  assert(!CGAL::is_closed(tm1));
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/elephant.off");
  input >> tm1;
  input.close();
  input.open("data-coref/sphere.off");
  input >> tm2;
  input.close();

  PMP::clip(tm1, tm2, params::clip_volume(true)
              .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
            params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(CGAL::is_closed(tm1));
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // test with a plane
  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();

  K::Plane_3 plane(0, 0, 1, -1);

  PMP::clip(tm1, plane, params::clip_volume(true));
  assert(CGAL::is_closed(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, plane, params::clip_volume(false)
              .use_compact_clipper(false));
  assert(!CGAL::is_closed(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, plane, params::clip_volume(false)
                        .use_compact_clipper(true));
  assert(CGAL::is_closed(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, K::Plane_3(-0.236474, 0.437732, 0.867451, -0.838791), params::clip_volume(true));
  assert(CGAL::is_closed(tm1));
  assert(!CGAL::is_empty(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, K::Plane_3(0, 0, 1, 2));
  assert(CGAL::is_empty(tm1));
  CGAL::clear(tm1);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, K::Plane_3(0, 0, 1, -2));
  assert(!CGAL::is_empty(tm1));
  CGAL::clear(tm1);

  // clipping with identity
  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(true)
                      .use_compact_clipper(true)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(num_vertices(tm1)==8);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .use_compact_clipper(false)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(CGAL::is_empty(tm1));
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .use_compact_clipper(true)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(num_vertices(tm1)==8);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION, K::Vector_3(1,0,0)), tm2);
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .use_compact_clipper(false)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(CGAL::is_empty(tm1));
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::transform(K::Aff_transformation_3(CGAL::TRANSLATION, K::Vector_3(1,0,0)), tm2);
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .use_compact_clipper(true)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==4);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // test orientation + patch without input vertex
  CGAL::make_tetrahedron(
    K::Point_3(0.53, -1.3, 0.2),
    K::Point_3(0.53, 1.1, 0.2),
    K::Point_3(0.53, -1.3, 0.4),
    K::Point_3(0.73, -1.3, 0.2),
    tm2);
  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==6);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  CGAL::make_tetrahedron(
    K::Point_3(0.53, -1.3, 0.2),
    K::Point_3(0.53, 1.1, 0.2),
    K::Point_3(0.53, -1.3, 0.4),
    K::Point_3(0.73, -1.3, 0.2),
    tm2);
  PMP::reverse_face_orientations(tm2);
  input.open("data-coref/cube.off");
  input >> tm1;
  input.close();
  PMP::clip(tm1, tm2,params::clip_volume(false)
                      .face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                     params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==6+8);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // clip meshes with intersection polyline opened
  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(1, 0, 0, -2));
  assert(vertices(tm1).size()==4);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(-1, 0, 0, 2));
  assert(vertices(tm1).size()==3);
  CGAL::clear(tm1);

  // test with clipper on border edge
  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 1, 0), K::Point_3(1, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(0, 1, 0 , 0));
  assert(vertices(tm1).size()==0);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 1, 0), K::Point_3(1, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(0, -1, 0 , 0));
  assert(vertices(tm1).size()==4);
  CGAL::clear(tm1);

  // test with clipper on border edge: full triangle
  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(0, 0, 1, 0), params::use_compact_clipper(true));
  assert(vertices(tm1).size()!=0);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 4, 0), K::Point_3(4, 0, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(0, 0, 1, 0), params::use_compact_clipper(false));
  assert(vertices(tm1).size()==0);
  CGAL::clear(tm1);

  // test tangencies
  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 2, 0), K::Point_3(1, 1, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(1, 0, 0, -1));
  assert(vertices(tm1).size()==3);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0, 0, 0), K::Point_3(0, 2, 0), K::Point_3(1, 1, 0), tm1 );
  PMP::clip(tm1, K::Plane_3(-1, 0, 0, 1));
  assert(vertices(tm1).size()==0);
  CGAL::clear(tm1);

  make_triangle( K::Point_3(0.5, 0, 0.5), K::Point_3(1, 0.5, 0.5), K::Point_3(0.5, 1, 0.5), tm1 );
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2, params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                      params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==3);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  make_triangle( K::Point_3(0.5, 0, 0.5), K::Point_3(1, 0.5, 0.5), K::Point_3(0.5, 1, 0.5), tm1 );
  input.open("data-coref/cube.off");
  input >> tm2;
  input.close();
  PMP::reverse_face_orientations(tm2);
  PMP::clip(tm1, tm2, params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                      params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(vertices(tm1).size()==0);
  CGAL::clear(tm1);
  CGAL::clear(tm2);

  // test special case
  input.open("data-clip/tm_1.off");
  input >> tm1;
  input.close();
  input.open("data-clip/clipper_1.off");
  input >> tm2;
  input.close();
  PMP::clip(tm1, tm2, params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)),
                      params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm2)));
  assert(is_valid_polygon_mesh(tm1));
  CGAL::clear(tm1);
  CGAL::clear(tm2);
}

template <class Mesh>
void test_split_plane()
{

  typedef K::Point_3                                     Point;

  typedef typename boost::graph_traits<Mesh>::face_descriptor          face_descriptor;
  typedef typename boost::graph_traits<Mesh>::edge_descriptor          edge_descriptor;

  typedef typename Mesh::template Property_map<face_descriptor, std::size_t> FCCmap;
  typedef typename Mesh::template Property_map<edge_descriptor, bool> Cst_edge_map;

  typedef CGAL::Face_filtered_graph<Mesh> Filtered_graph;

  Mesh tm;
  std::ifstream input("data/blobby_3cc.off");
  if (!input || !(input >> tm) || tm.is_empty()) {
    std::cerr << "Not a valid off file." << std::endl;
    return ;
  }


  // create a splitter mesh for the splitting plane using an internal CGAL function
  CGAL::Bbox_3 bbox = ::CGAL::Polygon_mesh_processing::bbox(tm);
  double xd=(std::max)(1.,(bbox.xmax()-bbox.xmin())/100);
  double yd=(std::max)(1.,(bbox.ymax()-bbox.ymin())/100);
  double zd=(std::max)(1.,(bbox.zmax()-bbox.zmin())/100);
  bbox=CGAL::Bbox_3(bbox.xmin()-xd, bbox.ymin()-yd, bbox.zmin()-zd,
                    bbox.xmax()+xd, bbox.ymax()+yd, bbox.zmax()+zd);

  typename K::Plane_3 plane(0, 0, 1, (bbox.zmin()+bbox.zmax())/2); // arbitratry plane
  Mesh splitter;
  CGAL::Oriented_side os = PMP::internal::clip_to_bbox(plane, bbox, splitter, PMP::parameters::all_default());


  if (os == CGAL::ON_ORIENTED_BOUNDARY)
  {
    // create a constrained edge map and corefine input mesh with the plane
    Cst_edge_map ecm = tm.template add_property_map<edge_descriptor, bool>("e:cst", false).first;
    PMP::corefine(tm, splitter, CGAL::parameters::edge_is_constrained_map(ecm));

    FCCmap fccmap = tm.template add_property_map<face_descriptor, std::size_t>("f:CC").first;
    std::size_t num = PMP::connected_components(tm, fccmap, PMP::parameters::edge_is_constrained_map(ecm));

    Filtered_graph ffg(tm, 0, fccmap);
    for (std::size_t i=0; i <num; ++i)
    {
      if (i!=0)
        ffg.set_selected_faces(i, fccmap);

      // create a new mesh for the component
      Mesh cc_mesh;
      CGAL::copy_face_graph(ffg, cc_mesh);

      // write the copy into a file
      std::ofstream(std::string("out_cc-")+std::to_string(i)+std::string(".off")) << std::setprecision(17) << cc_mesh;
    }
  }
  else
  {
    // nothing to do, no intersection.
    std::ofstream("out_cc-0.off") << std::setprecision(17) << tm;
  }

  return;
}

template <class TriangleMesh>
void test_split()
{

  // test with a clipper mesh
  TriangleMesh tm1, tm2;

  std::ifstream input("data-coref/elephant.off");
  input >> tm1;
  input.close();
  input.open("data-coref/sphere.off");
  input >> tm2;
  input.close();
  std::vector<TriangleMesh> output;
  PMP::split(tm1, tm2, std::back_inserter(output),
             params::face_index_map(get(CGAL::dynamic_face_property_t<std::size_t>(), tm1)));
int i=0;
  for(const auto& m : output)
  {
    std::ofstream(std::string("out_cc-")+std::to_string(i++)+std::string(".off")) << std::setprecision(17) << m;
    //assert(!CGAL::is_closed(m));
  }
  CGAL::clear(tm1);
  CGAL::clear(tm2);

}

int main()
{
  //test<Surface_mesh>();
  //test<Polyhedron>();
  test_split<Surface_mesh>();
  //test_split<Polyhedron>();

  return 0;
}


