#ifndef CGAL_TEST_PREFIX_H
#define CGAL_TEST_PREFIX_H

#define CGAL_USE_SURFACE_MESH
#define CGAL_USE_OPENMESH

#include <vector>
#include <fstream>

#include <CGAL/boost/graph/properties.h>
#include <boost/assign.hpp>
#include <boost/mpl/list.hpp>

#include <CGAL/Simple_cartesian.h>

#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>

#ifdef CGAL_USE_SURFACE_MESH
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh/IO.h>

#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>
#include <CGAL/boost/graph/properties_Surface_mesh.h>
#endif

// ATTN: If you change this kernel remember to also hack
// properties_PolyMesh_ArrayKernelT.h accordingly
typedef CGAL::Simple_cartesian<double> Kernel;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3> Polyhedron;

typedef Kernel::Point_3  Point_3;
typedef Kernel::Vector_3 Vector_3;

#ifdef CGAL_USE_SURFACE_MESH
typedef CGAL::Surface_mesh<Point_3> SM;
#endif

#if defined(CGAL_USE_OPENMESH)

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#define CGAL_BGL_TESTSUITE 1

#include <CGAL/boost/graph/graph_traits_PolyMesh_ArrayKernelT.h>
#include <CGAL/boost/graph/properties_PolyMesh_ArrayKernelT.h>

typedef OpenMesh::PolyMesh_ArrayKernelT</* MyTraits*/> OMesh;
#endif

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



static const char* data[] =
{ "data/7_faces_triangle.off", "data/genus3.off", "data/head.off",
  "data/hedra.off", "data/hedra_open.off",   "data/open_cube.off",
  "data/rombus.off", "data/tetrahedron.off", "data/triangle.off",
  "data/triangular_hole.off", "data/cube.off" };

#if defined(CGAL_USE_OPENMESH)
bool read_a_mesh(OMesh& s, const std::string& str) {
  return OpenMesh::IO::read_mesh(s, str);
}
#endif

#if defined(CGAL_USE_SURFACE_MESH)
// quick hack to generically read a file
bool read_a_mesh(SM& s, const std::string& str) {
  return CGAL::read_off(s, str);
}
#endif

bool read_a_mesh(Polyhedron& p, const std::string& str) {
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
  std::vector<T> vs;
  for(unsigned int i = 0; i < sizeof(data) / sizeof(data[0]); ++i) {
    vs.push_back(T());
    T& s = vs.back();
    if(!read_a_mesh(s, std::string(data[i])))
      throw std::runtime_error(std::string("Failed to read test data: ") + data[i]);
  }

  return vs;
}

std::vector<Polyhedron> poly_data() 
{ return t_data<Polyhedron>(); }

#if defined(CGAL_USE_SURFACE_MESH)
std::vector<SM> sm_data() 
{ return t_data<SM>(); }
#endif

#if defined(CGAL_USE_OPENMESH)
std::vector<OMesh> omesh_data() 
{ return t_data<OMesh>(); }
#endif

template <typename Graph>
struct Surface_fixture_1 {
  Surface_fixture_1() {
    assert(read_a_mesh(m, "data/fixture1.off"));
    assert(CGAL::is_valid(m));
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

    f1 = is_border(halfedge(u, m),m) ? face(opposite(halfedge(u, m), m), m) : face(halfedge(u, m), m);
    assert(f1 != boost::graph_traits<Graph>::null_face());
    CGAL::Halfedge_around_face_iterator<Graph> hafib, hafie;
    for(boost::tie(hafib, hafie) = halfedges_around_face(halfedge(f1, m), m); hafib != hafie; ++hafib) 
    {
      if(! is_border(opposite(*hafib, m), m))
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
    assert(read_a_mesh(m, "data/fixture2.off"));
    assert(CGAL::is_valid(m));

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
    assert(! is_border(h,m));
    f1 = face(h, m);
    assert(f1 != boost::graph_traits<Graph>::null_face());

    boost::tie(h, found) = halfedge(v, u, m);
    assert(found);
    assert(!is_border(h,m));
    f2 = face(h, m);
    assert(f2 != boost::graph_traits<Graph>::null_face());

    boost::tie(h, found) = halfedge(u, w, m);
    assert(found);
    assert(!is_border(h,m));
    f3 = face(h, m);
    assert(f3 != boost::graph_traits<Graph>::null_face());

    boost::tie(h, found) = halfedge(w, x, m);
    assert(found);
    assert(!is_border(h,m));
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
    assert(read_a_mesh(m, "data/fixture3.off"));
    assert(CGAL::is_valid(m));

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

    f1 = is_border(halfedge(u, m),m) ? face(opposite(halfedge(u, m), m), m) : face(halfedge(u, m), m);
    f2 = is_border(halfedge(u, m),m) ? face(opposite(halfedge(z, m), m), m) : face(halfedge(z, m), m);

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
    assert(read_a_mesh(m, "data/fixture4.off"));
    assert(CGAL::is_valid(m));

   typename boost::property_map<Graph, CGAL::vertex_point_t>::const_type
      pm = get(CGAL::vertex_point, const_cast<const Graph&>(m));

    int found = 0;
    typename boost::graph_traits<Graph>::halfedge_iterator hb, he;
    for(boost::tie(hb, he) = halfedges(m); hb != he; ++hb) {
      if(is_border(*hb,m)){
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
    assert(read_a_mesh(m, "data/add_face_to_border.off"));
    assert(CGAL::is_valid(m));

   typename boost::property_map<Graph, CGAL::vertex_point_t>::const_type
      pm = get(CGAL::vertex_point, const_cast<const Graph&>(m));

    int found = 0;
    typename boost::graph_traits<Graph>::halfedge_iterator hb, he;
    for(boost::tie(hb, he) = halfedges(m); hb != he; ++hb) {
      if(is_border(*hb,m)){
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
    assert(read_a_mesh(m, "data/quad.off"));
    assert(CGAL::is_valid(m));

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
    assert(read_a_mesh(m, "data/cube.off"));
    assert(CGAL::is_valid(m));

    h = *(halfedges(m).first);    
  }

  Graph m;
  typename boost::graph_traits<Graph>::halfedge_descriptor h;

};
template <typename Graph>
struct Surface_fixture_8 {
  Surface_fixture_8() {
    assert(read_a_mesh(m, "data/fixture5.off"));
    assert(CGAL::is_valid(m));

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
