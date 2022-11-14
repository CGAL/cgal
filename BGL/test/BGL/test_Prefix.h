#ifndef CGAL_TEST_PREFIX_H
#define CGAL_TEST_PREFIX_H

#include <vector>
#include <fstream>

#include <CGAL/boost/graph/properties.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>

#include <CGAL/Linear_cell_complex_for_bgl_combinatorial_map_helper.h>
#include <CGAL/boost/graph/graph_traits_Linear_cell_complex_for_combinatorial_map.h>
#include <CGAL/boost/graph/properties_Linear_cell_complex_for_combinatorial_map.h>

#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>

#include <CGAL/Constrained_triangulation_face_base_2.h>
#include <CGAL/Constrained_triangulation_plus_2.h>
#include <CGAL/Triangulation_data_structure_2.h>
#include <CGAL/Triangulation_vertex_base_with_id_2.h>
#include <CGAL/Triangulation_face_base_with_id_2.h>
#include <CGAL/Triangulation_hierarchy_vertex_base_2.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_2.h>
#include <CGAL/boost/graph/properties_Triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Triangulation_hierarchy_2.h>
#include <CGAL/boost/graph/properties_Triangulation_hierarchy_2.h>
#include <CGAL/boost/graph/graph_traits_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/properties_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Regular_triangulation_2.h>
#include <CGAL/boost/graph/properties_Regular_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_2.h>
#include <CGAL/boost/graph/properties_Constrained_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/properties_Constrained_Delaunay_triangulation_2.h>
#include <CGAL/boost/graph/graph_traits_Constrained_triangulation_plus_2.h>
#include <CGAL/boost/graph/properties_Constrained_triangulation_plus_2.h>
#include <CGAL/boost/graph/Seam_mesh.h>

#include <CGAL/boost/graph/io.h>

// ATTN: If you change this kernel remember to also hack
// properties_PolyMesh_ArrayKernelT.h accordingly
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef Kernel::Point_3  Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Triangle_3 Triangle_3;
typedef Kernel::Iso_cuboid_3 Iso_cuboid_3;

typedef CGAL::Linear_cell_complex_traits<3, Kernel> MyTraits;
typedef CGAL::Linear_cell_complex_for_bgl_combinatorial_map_helper
          <2, 3, MyTraits>::type LCC;

typedef CGAL::Surface_mesh<Point_3> SM;

typedef SM::Property_map<SM::Edge_index, bool>                  Seam_edge_pmap;
typedef SM::Property_map<SM::Vertex_index, bool>                Seam_vertex_pmap;
typedef CGAL::Seam_mesh<SM, Seam_edge_pmap, Seam_vertex_pmap>   Seam_mesh;

#if defined(CGAL_USE_OPENMESH)

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#define CGAL_BGL_TESTSUITE 1

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> OMesh;
#endif

typedef CGAL::Triangulation_vertex_base_with_id_2<Kernel>        Vbb;
typedef CGAL::Triangulation_face_base_with_id_2<Kernel>          Fbb;

typedef CGAL::Triangulation_2<Kernel>                            Triangulation_no_id_2;

typedef CGAL::Triangulation_2<Kernel,
          CGAL::Triangulation_data_structure_2<Vbb, Fbb> >       Triangulation_2;
typedef CGAL::Delaunay_triangulation_2<Kernel,
          CGAL::Triangulation_data_structure_2<Vbb, Fbb> >       Delaunay_triangulation_2;

typedef CGAL::Regular_triangulation_vertex_base_2<Kernel, Vbb>   RVb;
typedef CGAL::Regular_triangulation_face_base_2<Kernel, Fbb>     RFb;
typedef CGAL::Regular_triangulation_2<Kernel,
          CGAL::Triangulation_data_structure_2<RVb, RFb> >       Regular_triangulation_2;

typedef CGAL::Constrained_triangulation_face_base_2<Kernel, Fbb> CDFb;
typedef CGAL::Triangulation_hierarchy_vertex_base_2<Vbb>         CDVb;
typedef CGAL::Constrained_triangulation_2<Kernel,
          CGAL::Triangulation_data_structure_2<CDVb, CDFb> >     Constrained_triangulation_2;
typedef CGAL::Constrained_Delaunay_triangulation_2<Kernel,
          CGAL::Triangulation_data_structure_2<CDVb, CDFb> >     Constrained_Delaunay_triangulation_2;

typedef CGAL::Constrained_triangulation_plus_2<
          Constrained_Delaunay_triangulation_2>                  CDT_P2;

typedef CGAL::Triangulation_hierarchy_2<CDT_P2>                  Triangulation_hierarchy_2;

#include <CGAL/boost/graph/helpers.h>

// helper to easily define all graph_traits members
#define CGAL_GRAPH_TRAITS_MEMBERS(T)                                    \
  typedef boost::graph_traits< T >                         Traits;               \
  typedef typename Traits::vertex_descriptor               vertex_descriptor;    \
  typedef typename Traits::edge_descriptor                 edge_descriptor;      \
  typedef typename Traits::halfedge_descriptor             halfedge_descriptor;  \
  typedef typename Traits::face_descriptor                 face_descriptor;      \
  typedef typename Traits::face_iterator                   face_iterator;        \
  typedef typename Traits::faces_size_type                 faces_size_type;      \
  typedef typename Traits::out_edge_iterator               out_edge_iterator;    \
  typedef typename Traits::in_edge_iterator                in_edge_iterator;     \
  typedef typename Traits::vertex_iterator                 vertex_iterator;      \
  typedef typename Traits::vertices_size_type              vertices_size_type;   \
  typedef typename Traits::edge_iterator                   edge_iterator;        \
  typedef typename Traits::edges_size_type                 edges_size_type;      \
  typedef typename Traits::halfedge_iterator               halfedge_iterator;    \
  typedef typename Traits::halfedges_size_type             halfedges_size_type;  \
  typedef CGAL::Halfedge_around_face_iterator<T>           halfedge_around_face_iterator; \
  typedef CGAL::Halfedge_around_target_iterator<T>           halfedge_around_target_iterator; \
  CGAL_USE_TYPE(vertex_descriptor); \
  CGAL_USE_TYPE(halfedge_descriptor); \
  CGAL_USE_TYPE(edge_descriptor); \
  CGAL_USE_TYPE(face_descriptor); \
  CGAL_USE_TYPE(vertices_size_type); \
  CGAL_USE_TYPE(halfedges_size_type); \
  CGAL_USE_TYPE(edges_size_type); \
  CGAL_USE_TYPE(faces_size_type); \
  CGAL_USE_TYPE(vertex_iterator); \
  CGAL_USE_TYPE(halfedge_iterator); \
  CGAL_USE_TYPE(edge_iterator); \
  CGAL_USE_TYPE(face_iterator);   \
  CGAL_USE_TYPE(out_edge_iterator); \
  CGAL_USE_TYPE(in_edge_iterator); \
  CGAL_USE_TYPE(halfedge_around_target_iterator); \
  CGAL_USE_TYPE(halfedge_around_face_iterator); \
  do { } while(0)




/*
#if defined(CGAL_USE_OPENMESH)
bool read_a_mesh(OMesh& s, const std::string& str) {
  return OpenMesh::IO::read_mesh(s, str);
}
#endif
*/

template<typename T>
bool read_a_mesh(T& m, const std::string& str)
{
  return CGAL::IO::read_OFF(str, m);
}

bool read_a_mesh(Polyhedron& p, const std::string& str)
{
  std::ifstream in(str.c_str());
  in >> p;
  bool success = in.good();
  if(success)
    set_halfedgeds_items_id(p);
  return success;
}

template <typename T>
std::vector<T> t_data()
{
  static const std::string data[] =
    { "data/7_faces_triangle.off", "data/genus3.off", CGAL::data_file_path("meshes/head.off"),
      CGAL::data_file_path("meshes/hedra.off"), CGAL::data_file_path("meshes/hedra_open.off"), CGAL::data_file_path("meshes/open_cube.off"),
      "data/rombus.off", "data/tetrahedron.off", "data/triangle.off",
      CGAL::data_file_path("meshes/triangular_hole.off"), CGAL::data_file_path("meshes/cube.off") };

  std::vector<T> vs;
  for(unsigned int i = 0; i < sizeof(data) / sizeof(data[0]); ++i) {
    vs.push_back(T());
    T& s = vs.back();
    if(!read_a_mesh(s, std::string(data[i])))
      throw std::runtime_error(std::string("Failed to read test data: ") + data[i]);
  }

  return vs;
}

std::vector<Polyhedron> poly_data() { return t_data<Polyhedron>(); }
std::vector<SM> sm_data() { return t_data<SM>(); }
std::vector<LCC> lcc_data() { return t_data<LCC>(); }

#if defined(CGAL_USE_OPENMESH)
std::vector<OMesh> omesh_data() { return t_data<OMesh>(); }
#endif

template <typename Tr>
Tr build_dummy_triangulation()
{
  typedef typename Tr::Point                                       Point;

  Tr t;
  t.insert(Point(0.1,0));
  t.insert(Point(1,0));
  t.insert(Point(0.2,0.2));
  t.insert(Point(0,1));
  t.insert(Point(0,2));

  return t;
}

template <typename Tr>
Tr build_dummy_triangulation_with_ids()
{
  Tr t = build_dummy_triangulation<Tr>();
  CGAL::set_triangulation_ids(t);
  return t;
}

Triangulation_no_id_2 t2_no_id_data() { return build_dummy_triangulation<Triangulation_no_id_2>(); }
Triangulation_2 t2_data() { return build_dummy_triangulation_with_ids<Triangulation_2>(); }
Delaunay_triangulation_2 dt2_data() { return build_dummy_triangulation_with_ids<Delaunay_triangulation_2>(); }
Regular_triangulation_2 rt2_data() { return build_dummy_triangulation_with_ids<Regular_triangulation_2>(); }
Constrained_triangulation_2 ct2_data() { return build_dummy_triangulation_with_ids<Constrained_triangulation_2>(); }
Constrained_Delaunay_triangulation_2 cdt2_data() { return build_dummy_triangulation_with_ids<Constrained_Delaunay_triangulation_2>(); }
CDT_P2 cdtp2_data() { return build_dummy_triangulation_with_ids<CDT_P2>(); }
Triangulation_hierarchy_2 t2h_data() { return build_dummy_triangulation_with_ids<Triangulation_hierarchy_2>(); }

template <typename Graph>
struct Surface_fixture_1 {
  Surface_fixture_1() {
    const bool is_reading_successful = read_a_mesh(m, "data/fixture1.off");
    assert(is_reading_successful);
    assert(CGAL::is_valid_polygon_mesh(m));
    typename boost::property_map<Graph, CGAL::vertex_point_t>::const_type
      pm = get(CGAL::vertex_point, const_cast<const Graph&>(m));

    typename boost::graph_traits<Graph>::vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = vertices(m); vb != ve; ++vb) {
      if     (get(pm, *vb) == Point_3(0, 0, 0))
        u = *vb;
      else if(get(pm, *vb) == Point_3(1, 0, 0))
        v = *vb;
      else if(get(pm, *vb) == Point_3(0, 1, 0))
        w = *vb;
      else if(get(pm, *vb) == Point_3(1, 1, 0))
        x = *vb;
      else if(get(pm, *vb) == Point_3(2, 0, 0))
        y = *vb;
    }
    assert(u != boost::graph_traits<Graph>::null_vertex());
    assert(v != boost::graph_traits<Graph>::null_vertex());
    assert(w != boost::graph_traits<Graph>::null_vertex());
    assert(x != boost::graph_traits<Graph>::null_vertex());
    assert(y != boost::graph_traits<Graph>::null_vertex());

    f1 = CGAL::is_border(halfedge(u, m),m) ? face(opposite(halfedge(u, m), m), m) : face(halfedge(u, m), m);
    assert(f1 != boost::graph_traits<Graph>::null_face());
    CGAL::Halfedge_around_face_iterator<Graph> hafib, hafie;
    for(boost::tie(hafib, hafie) = CGAL::halfedges_around_face(halfedge(f1, m), m); hafib != hafie; ++hafib)
    {
      if(! CGAL::is_border(opposite(*hafib, m), m))
        f2 = face(opposite(*hafib, m), m);
    }
    typename boost::graph_traits<Graph>::face_iterator fb, fe;
    for(boost::tie(fb, fe) = faces(m); fb != fe; ++fb) {
      if(*fb != f1 && *fb != f2)
        f3 = *fb;
    }
    assert(f2 != boost::graph_traits<Graph>::null_face());
    assert(f3 != boost::graph_traits<Graph>::null_face());
  }

  Graph m;
  typename boost::graph_traits<Graph>::vertex_descriptor u, v, w, x, y;
  typename boost::graph_traits<Graph>::face_descriptor f1, f2, f3;
};

template <typename Graph>
struct Surface_fixture_2 {
  Surface_fixture_2() {
    const bool is_reading_successful = read_a_mesh(m, "data/fixture2.off");
    assert(is_reading_successful);
    assert(CGAL::is_valid_polygon_mesh(m));

    typename boost::property_map<Graph, CGAL::vertex_point_t>::const_type
      pm = get(CGAL::vertex_point, const_cast<const Graph&>(m));

    typename boost::graph_traits<Graph>::vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = vertices(m); vb != ve; ++vb) {
      if     (get(pm, *vb) == Point_3(0, 2, 0))
        u = *vb;
      else if(get(pm, *vb) == Point_3(2, 2, 0))
        v = *vb;
      else if(get(pm, *vb) == Point_3(0, 0, 0))
        w = *vb;
      else if(get(pm, *vb) == Point_3(2, 0, 0))
        x = *vb;
      else if(get(pm, *vb) == Point_3(1, 1, 0))
        y = *vb;
    }
    assert(u != boost::graph_traits<Graph>::null_vertex());
    assert(v != boost::graph_traits<Graph>::null_vertex());
    assert(w != boost::graph_traits<Graph>::null_vertex());
    assert(x != boost::graph_traits<Graph>::null_vertex());
    assert(y != boost::graph_traits<Graph>::null_vertex());
    typename boost::graph_traits<Graph>::halfedge_descriptor h;
    bool found;
    boost::tie(h, found) = halfedge(x, v, m);
    assert(found);
    assert(! CGAL::is_border(h,m));
    f1 = face(h, m);
    assert(f1 != boost::graph_traits<Graph>::null_face());

    boost::tie(h, found) = halfedge(v, u, m);
    assert(found);
    assert(!CGAL::is_border(h,m));
    f2 = face(h, m);
    assert(f2 != boost::graph_traits<Graph>::null_face());

    boost::tie(h, found) = halfedge(u, w, m);
    assert(found);
    assert(!CGAL::is_border(h,m));
    f3 = face(h, m);
    assert(f3 != boost::graph_traits<Graph>::null_face());

    boost::tie(h, found) = halfedge(w, x, m);
    assert(found);
    assert(!CGAL::is_border(h,m));
    f4 = face(h, m);
    assert(f4 != boost::graph_traits<Graph>::null_face());
  }

  Graph m;
  typename boost::graph_traits<Graph>::vertex_descriptor u, v, w, x, y;
  typename boost::graph_traits<Graph>::face_descriptor f1, f2, f3, f4;

  ~Surface_fixture_2() {}
};

template <typename Graph>
struct Surface_fixture_3 {
  Surface_fixture_3() {
    const bool is_reading_successful = read_a_mesh(m, "data/fixture3.off");
    assert(is_reading_successful);
    assert(CGAL::is_valid_polygon_mesh(m));

    typename boost::property_map<Graph, CGAL::vertex_point_t>::const_type
      pm = get(CGAL::vertex_point, const_cast<const Graph&>(m));

    typename boost::graph_traits<Graph>::vertex_iterator vb, ve;
    for(boost::tie(vb, ve) = vertices(m); vb != ve; ++vb) {
      if     (get(pm, *vb) == Point_3(0, 1, 0))
        u = *vb;
      else if(get(pm, *vb) == Point_3(0, 0, 0))
        v = *vb;
      else if(get(pm, *vb) == Point_3(1, 0, 0))
        w = *vb;
      else if(get(pm, *vb) == Point_3(1, 1, 0))
        x = *vb;
      else if(get(pm, *vb) == Point_3(2, 0, 0))
        y = *vb;
      else if(get(pm, *vb) == Point_3(2, 1, 0))
        z = *vb;
    }
    assert(u != boost::graph_traits<Graph>::null_vertex());
    assert(v != boost::graph_traits<Graph>::null_vertex());
    assert(w != boost::graph_traits<Graph>::null_vertex());
    assert(x != boost::graph_traits<Graph>::null_vertex());
    assert(y != boost::graph_traits<Graph>::null_vertex());
    assert(z != boost::graph_traits<Graph>::null_vertex());

    f1 = CGAL::is_border(halfedge(u, m),m) ? face(opposite(halfedge(u, m), m), m) : face(halfedge(u, m), m);
    f2 = CGAL::is_border(halfedge(z, m),m) ? face(opposite(halfedge(z, m), m), m) : face(halfedge(z, m), m);

    assert(f1 != boost::graph_traits<Graph>::null_face());
    assert(f2 != boost::graph_traits<Graph>::null_face());

  }

  Graph m;
  typename boost::graph_traits<Graph>::vertex_descriptor u, v, w, x, y, z;
  typename boost::graph_traits<Graph>::face_descriptor f1, f2;

  ~Surface_fixture_3() {}
};

template <typename Graph>
struct Surface_fixture_4 {
  Surface_fixture_4() {
    const bool is_reading_successful = read_a_mesh(m, "data/fixture4.off");
    assert(is_reading_successful);
    assert(CGAL::is_valid_polygon_mesh(m));

   typename boost::property_map<Graph, CGAL::vertex_point_t>::const_type
      pm = get(CGAL::vertex_point, const_cast<const Graph&>(m));

    int found = 0;
    typename boost::graph_traits<Graph>::halfedge_iterator hb, he;
    for(boost::tie(hb, he) = halfedges(m); hb != he; ++hb) {
      if(CGAL::is_border(*hb,m)){
        if(get(pm, target(*hb,m)) == Point_3(0,0,0)){
          if(found == 0){
            h1 = *hb;
            ++found;
          } else if(found == 1){
            h2 = *hb;
            ++found;
          }
        }
      }
    }
    assert(found == 2);
  }

  Graph m;
  typename boost::graph_traits<Graph>::halfedge_descriptor h1, h2;



};


template <typename Graph>
struct Surface_fixture_5 {
  Surface_fixture_5() {
    const bool is_reading_successful = read_a_mesh(m, "data/add_face_to_border.off");
    assert(is_reading_successful);
    assert(CGAL::is_valid_polygon_mesh(m));

   typename boost::property_map<Graph, CGAL::vertex_point_t>::const_type
      pm = get(CGAL::vertex_point, const_cast<const Graph&>(m));

    int found = 0;
    typename boost::graph_traits<Graph>::halfedge_iterator hb, he;
    for(boost::tie(hb, he) = halfedges(m); hb != he; ++hb) {
      if(CGAL::is_border(*hb,m)){
        if(get(pm, target(*hb,m)) == Point_3(2,1,0)){
          h1 = *hb;
          found++;
        } else if(get(pm, target(*hb,m)) == Point_3(2,-1,0)){
          h2 = *hb;
          found++;
        }
      }
    }
    assert(found == 2);
  }

  Graph m;
  typename boost::graph_traits<Graph>::halfedge_descriptor h1, h2;

};

template <typename Graph>
struct Surface_fixture_6 {
  Surface_fixture_6() {
    const bool is_reading_successful = read_a_mesh(m, "data/quad.off");
    assert(is_reading_successful);
    assert(CGAL::is_valid_polygon_mesh(m));

    typename boost::graph_traits<Graph>::halfedge_descriptor h;

    h1 = halfedge(*faces(m).first, m);

    h2 = next(next(h1,m),m);
  }

  Graph m;
  typename boost::graph_traits<Graph>::halfedge_descriptor h1, h2;

};


template <typename Graph>
struct Surface_fixture_7 {
  Surface_fixture_7() {
    const bool is_reading_successful = read_a_mesh(m, CGAL::data_file_path("meshes/cube.off"));
    assert(is_reading_successful);
    assert(CGAL::is_valid_polygon_mesh(m));

    h = *(halfedges(m).first);
  }

  Graph m;
  typename boost::graph_traits<Graph>::halfedge_descriptor h;

};
template <typename Graph>
struct Surface_fixture_8 {
  Surface_fixture_8() {
    const bool is_reading_successful = read_a_mesh(m, "data/fixture5.off");
    assert(is_reading_successful);
    assert(CGAL::is_valid_polygon_mesh(m));

   typename boost::property_map<Graph, CGAL::vertex_point_t>::const_type
      pm = get(CGAL::vertex_point, const_cast<const Graph&>(m));

    int found = 0;
    typename boost::graph_traits<Graph>::halfedge_iterator hb, he;
    for(boost::tie(hb, he) = halfedges(m); hb != he; ++hb) {
      if(get(pm, source(*hb,m)) == Point_3(0,0,0) &&
         get(pm, target(*hb,m)) == Point_3(1,0,0)){
          h1 = *hb;
          found++;
      } else if(get(pm, source(*hb,m)) == Point_3(1,0,0) &&
                get(pm, target(*hb,m)) == Point_3(0,1,0)){
        h2 = *hb;
        found++;
      }  else if(get(pm, source(*hb,m)) == Point_3(0,1,0) &&
                get(pm, target(*hb,m)) == Point_3(0,0,0)){
        h3 = *hb;
        found++;
      }
    }

    assert(found == 3);
  }

  Graph m;
  typename boost::graph_traits<Graph>::halfedge_descriptor h1, h2, h3;

};


#endif /* CGAL_TEST_PREFIX_H */
