
#include "test_Prefix.h"
#include <boost/range/distance.hpp>
#include <CGAL/boost/graph/Euler_operations.h>

#include <CGAL/IO/OFF.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <CGAL/boost/graph/copy_face_graph.h>

template <typename T>
void
test_copy_face_graph_nm_umbrella()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  T g;
  Kernel::Point_3 p(0,0,0);

  CGAL::make_tetrahedron(p, p, p, p, g);
  CGAL::make_tetrahedron(p, p, p, p, g);

  std::vector<vertex_descriptor> verts(vertices(g).begin(), vertices(g).end());

  //merge verts[0] and verts[4]
  for (halfedge_descriptor h : CGAL::halfedges_around_target(verts[4], g))
    set_target(h, verts[0], g);
  remove_vertex(verts[4], g);

  T g_copy;
  CGAL::copy_face_graph(g, g_copy);

  for (halfedge_descriptor h : halfedges(g_copy))
  {
    assert( target(h, g_copy) != Traits::null_vertex() );
  }
}

template <typename T>
void
test_copy_face_graph_isolated_vertices()
{
  Kernel::Point_3 p(0,0,0);
  {
    T s, t;
    add_vertex(s);
    CGAL::copy_face_graph(s, t);
  }

  {
    T s, t;
    add_vertex(t);
    CGAL::copy_face_graph(s, t);
  }

  {
    T s, t;
    CGAL::make_triangle(p, p, p, s);
    add_vertex(s);
    t=s;
    CGAL::copy_face_graph(s, t);
  }

  {
    T s, t;
    CGAL::make_triangle(p, p, p, s);
    add_vertex(s);
    add_vertex(t);
    CGAL::copy_face_graph(s, t);
  }

  {
    T s, t;
    CGAL::make_tetrahedron(p, p, p, p, s);
    add_vertex(s);
    t=s;
    CGAL::copy_face_graph(s, t);
  }

  {
    T s, t;
    CGAL::make_tetrahedron(p, p, p, p, s);
    add_vertex(s);
    add_vertex(t);
    CGAL::copy_face_graph(s, t);
  }
}

template <typename T>
void
join_face_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_1<T> f;

  bool found;
  halfedge_descriptor e;
  boost::tie(e, found) = halfedge(f.w, f.v, f.m);
  assert(found);
  // manually set the halfedge of f.f1 to the edge that is to be
  // removed to provoke a special case
  set_halfedge(f.f1, e, f.m);
  CGAL::Euler::join_face(e,f.m);

  assert(CGAL::internal::exact_num_faces(f.m) == 2);
  assert(CGAL::internal::exact_num_edges(f.m) == 6);

  CGAL::Halfedge_around_face_iterator<T> begin, end;
  boost::tie(begin, end) = CGAL::halfedges_around_face(halfedge(f.f1, f.m), f.m);
  assert(std::distance(begin, end) == 4);
  for(; begin != end; ++begin)
  {

    halfedge_descriptor hd = *begin;
    assert(face(hd, f.m) == f.f1);

  }

  face_iterator fit, fend;
  for(boost::tie(fit, fend) = faces(f.m); fit != fend; ++fit) {
    assert(*fit == f.f1 || *fit == f.f3);
  }

  assert(degree(f.w, f.m) == 2);
  assert(degree(f.v, f.m) == 3);
  assert(CGAL::is_valid_polygon_mesh(f.m));

}



template <typename T>
void
remove_face_test_1()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
  Surface_fixture_1<T> f;

  // find the edge between x and y
  bool found;
  halfedge_descriptor e;
  boost::tie(e, found) = halfedge(f.x, f.y, f.m);
  assert(found);
  assert(face(e, f.m) == f.f3);

  CGAL::Euler::remove_face(e,f.m);

  assert(CGAL::is_valid_polygon_mesh(f.m));

  assert_EQUAL(degree(f.v, f.m) == 3);
  assert_EQUAL(degree(f.x, f.m) == 2);
  assert_EQUAL(CGAL::internal::exact_num_faces(f.m) == 2);
  assert_EQUAL(CGAL::internal::exact_num_edges(f.m) == 5);
  assert_EQUAL(CGAL::internal::exact_num_vertices(f.m) == 4);
  halfedge_iterator eb, ee;
  int count = 0;
  for(boost::tie(eb, ee) = halfedges(f.m); eb != ee; ++eb) {
    if(face(*eb,f.m) == boost::graph_traits<T>::null_face())
      ++count;
  }
  assert(count == 4);
}



template <typename T>
void
remove_face_test_2()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_2<T> f;

  // find the edge between x and v
  bool found;
  halfedge_descriptor e;

  boost::tie(e, found) = halfedge(f.x, f.w, f.m);
  assert(found);
  boost::tie(e, found) = halfedge(f.x, f.v, f.m);
  assert(found);
  assert(face(e, f.m) == f.f1);
  CGAL::Euler::remove_face(e,f.m);
  assert(CGAL::is_valid_polygon_mesh(f.m));

  assert(CGAL::internal::exact_num_faces(f.m) == 3);
  assert(CGAL::internal::exact_num_edges(f.m) == 7);
  assert(CGAL::internal::exact_num_vertices(f.m) == 5);

  boost::tie(e, found) = halfedge(f.x, f.w, f.m);
  assert(found);
  assert(face(e,f.m) == boost::graph_traits<T>::null_face());

  // check the boundary

  halfedge_descriptor n = next(e, f.m);
  while(n != e) {
    assert(face(n,f.m) == boost::graph_traits<T>::null_face() );
    n = next(n, f.m);
  }
}

template <typename T>
void
add_face_to_border_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_5<T> f;

  CGAL::Euler::add_face_to_border(f.h1, f.h2, f.m);

  assert(CGAL::is_valid_polygon_mesh(f.m));

}

template <typename T>
void
add_vertex_and_face_to_border_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_5<T> f;
  typedef typename boost::graph_traits<T>::halfedge_descriptor  halfedge_descriptor;
  halfedge_descriptor h1 = f.h1, h2 = f.h2;
  const T& m = f.m;

  int dist=0;
  halfedge_descriptor hd = h1;
  while(hd != h2){
    ++dist;
    hd = next(hd,m);
  }
  assert(dist == 2);

  int blength = 0;
  for(halfedge_descriptor hd : CGAL::halfedges_around_face(h1,m)){
    CGAL_USE(hd);
    blength++;
  }

  halfedge_descriptor res = CGAL::Euler::add_vertex_and_face_to_border(f.h1, f.h2, f.m);
  assert(CGAL::is_valid_polygon_mesh(f.m));

  assert(! CGAL::is_border(res,m));
  assert(CGAL::is_border(opposite(res,m),m));
  res = opposite(res,m);
  for(halfedge_descriptor hd : CGAL::halfedges_around_face(res,m)){
    CGAL_USE(hd);
    blength--;
  }
  assert(blength == 0);

}


template <typename T>
void
join_vertex_interior_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_3<T> f;
  halfedge_descriptor e;

  bool found;
  boost::tie(e, found) = halfedge(f.w, f.x, f.m);
  assert(found);
  CGAL::Euler::join_vertex(e,f.m);
  assert(CGAL::internal::exact_num_faces(f.m) == 2);
  assert(CGAL::internal::exact_num_vertices(f.m) == 5);
  assert(CGAL::internal::exact_num_edges(f.m) == 6);
  assert(CGAL::halfedges_around_face(halfedge(f.f1, f.m), f.m).size() == 3);
  assert(CGAL::halfedges_around_face(halfedge(f.f2, f.m), f.m).size() == 3);
  assert(degree(f.x, f.m) == 4);
  assert(CGAL::is_valid_polygon_mesh(f.m));
}

template <typename T>
void
join_vertex_exterior_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  {
    // exterior edge is collapsed
    Surface_fixture_3<T> f;
    halfedge_descriptor e;
    bool found;
    boost::tie(e, found) = halfedge(f.w, f.y, f.m);
    assert(source(e,f.m) == f.w);
    assert(target(e,f.m) == f.y);
    assert(found);
    CGAL::Euler::join_vertex(e,f.m);
    assert(CGAL::internal::exact_num_faces(f.m) == 2);
    assert(CGAL::internal::exact_num_vertices(f.m) == 5);
    assert(CGAL::internal::exact_num_edges(f.m) == 6);
    assert(CGAL::halfedges_around_face(halfedge(f.f1, f.m), f.m).size() == 4);
    assert(CGAL::halfedges_around_face(halfedge(f.f2, f.m), f.m).size() == 3);
    assert(degree(f.y, f.m) == 3);
    assert(CGAL::is_valid_polygon_mesh(f.m));
  }

  {
    Surface_fixture_3<T> f;
    halfedge_descriptor e;
    bool found;
    boost::tie(e, found) = halfedge(f.y, f.w, f.m);

    assert(source(e,f.m) == f.y);
    assert(target(e,f.m) == f.w);
    assert(found);
    CGAL::Euler::join_vertex(e,f.m);
    assert(CGAL::internal::exact_num_faces(f.m) == 2);
    assert(CGAL::internal::exact_num_vertices(f.m) == 5);
    assert(CGAL::internal::exact_num_edges(f.m) == 6);
    assert(CGAL::halfedges_around_face(halfedge(f.f1, f.m), f.m).size() == 4);
    assert(CGAL::halfedges_around_face(halfedge(f.f2, f.m), f.m).size() == 3);

    assert(CGAL::is_valid_polygon_mesh(f.m));
    assert(degree(f.w, f.m) == 3);

  }
}


template <typename T>
void
split_vertex()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);

  Surface_fixture_3<T> f;
  halfedge_descriptor h1, h2;
  bool found;
  boost::tie(h1, found) = halfedge(f.w, f.y, f.m);
  assert(found);
  boost::tie(h2, found) = halfedge(f.z, f.y, f.m);
  assert(found);
  assert(face(h2, f.m) == Traits::null_face());

  // split border vertex y
  CGAL::Euler::split_vertex(h1, h2,f.m);
  assert(CGAL::is_valid_polygon_mesh(f.m));
  assert(CGAL::internal::exact_num_vertices(f.m) == 7);
  assert(CGAL::internal::exact_num_edges(f.m) == 8);
  assert(CGAL::halfedges_around_face(h1, f.m).size() == 5);
  assert(CGAL::halfedges_around_face(h2, f.m).size() == 7);
}

template <typename T>
void
split_join_vertex_inverse()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
  Surface_fixture_3<T> f;
  halfedge_descriptor h, h1, h2;
  bool found;
  boost::tie(h, found) = halfedge(f.w, f.x, f.m);
  assert(found);
  CGAL::Euler::join_vertex(h,f.m);
  assert(CGAL::is_valid_polygon_mesh(f.m));
  boost::tie(h1, found) = halfedge(f.z, f.x, f.m);
  assert(found);
  boost::tie(h2, found) = halfedge(f.v, f.x, f.m);
  assert(found);
  CGAL::Euler::join_vertex(CGAL::Euler::split_vertex(h1, h2,f.m),f.m);
  assert(CGAL::is_valid_polygon_mesh(f.m));

  assert(CGAL::internal::exact_num_vertices(f.m)== 5);
  assert(CGAL::internal::exact_num_faces(f.m) == 2);
  assert(CGAL::internal::exact_num_edges(f.m) == 6);
  assert(CGAL::internal::exact_num_halfedges(f.m) == 12);
  assert(CGAL::halfedges_around_face(h1, f.m).size() == 3);
  assert(CGAL::halfedges_around_face(h2, f.m).size() == 3);
}


template <typename T>
void
join_loop_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
  Surface_fixture_4<T> f;

  CGAL::Euler::join_loop(f.h1, f.h2, f.m);

  assert(CGAL::is_valid_polygon_mesh(f.m));
}

template <typename T>
void
split_loop_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(T);
  Surface_fixture_8<T> f;

  CGAL::Euler::split_loop(f.h1, f.h2, f.h3, f.m);
  assert(CGAL::internal::exact_num_vertices(f.m) == 8);
  assert(CGAL::internal::exact_num_faces(f.m) == 8);
  assert(CGAL::internal::exact_num_halfedges(f.m) == 24);
  assert(CGAL::is_valid_polygon_mesh(f.m));
}

template <typename T>
void
split_face_test()
{
 Surface_fixture_6<T> f;
 CGAL::Euler::split_face(f.h1, f.h2,f.m);
 assert(num_vertices(f.m) == 4);
 assert(num_faces(f.m) == 2);
 assert(num_halfedges(f.m) == 10);
}

template <typename T>
void
 make_hole_test()
{
 Surface_fixture_7<T> f;
 std::size_t nv = num_vertices(f.m);
 std::size_t nf = num_faces(f.m);
 std::size_t nh = num_halfedges(f.m);

 CGAL::Euler::make_hole(f.h, f.m);

 assert(CGAL::internal::exact_num_vertices(f.m) == nv);
 assert(CGAL::internal::exact_num_faces(f.m) == nf-1 );
 assert(CGAL::internal::exact_num_halfedges(f.m) == nh);
}

template <typename T>
void
remove_center_vertex_test()
{
 Surface_fixture_7<T> f;
 std::size_t nv = num_vertices(f.m);
 std::size_t nf = num_faces(f.m);
 std::size_t nh = num_halfedges(f.m);

 typename boost::graph_traits<T>::degree_size_type deg = degree(target(f.h,f.m),f.m);
 CGAL::Euler::remove_center_vertex(f.h,f.m);

 assert(CGAL::internal::exact_num_vertices(f.m) == nv-1);
 assert(CGAL::internal::exact_num_faces(f.m) == (nf-deg)+1);
 assert(CGAL::internal::exact_num_halfedges(f.m) == nh-(2*deg));
}

template <typename T>
void
join_split_inverse()
{

}

template <typename T>
void
does_satisfy_link_condition()
{
  Surface_fixture_7<T> f;

  assert(CGAL::Euler::does_satisfy_link_condition(*edges(f.m).first,f.m));
}

template <typename Graph>
void
test_swap_edges()
{
  typedef typename boost::graph_traits<Graph>::halfedge_descriptor halfedge_descriptor;
  std::size_t nbh=12;
  Kernel::Point_3 pt(0,0,0);
  // test all possible pairs of halfedges
  for (std::size_t i=0; i<nbh-1; ++i)
  {
    for(std::size_t j=i+1; j<nbh; ++j)
    {
      Graph g;
      CGAL::make_tetrahedron(pt,pt,pt,pt,g);
      halfedge_descriptor h1 = *std::next(std::begin(halfedges(g)), i);
      halfedge_descriptor h2 = *std::next(std::begin(halfedges(g)), j);
      CGAL::internal::swap_edges(h1, h2, g);
      assert(CGAL::is_valid_polygon_mesh(g));
    }
  }
}

template <typename T>
void
add_face_bug()
{
  typedef boost::graph_traits<T> GT;
  typedef typename GT::vertex_descriptor vertex_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;

  T g;

  std::vector<vertex_descriptor> vs;
  vs.push_back( add_vertex(g) ); // Kernel::Point_3(0,1,0)
  vs.push_back( add_vertex(g) ); // Kernel::Point_3(4,1,0)
  vs.push_back( add_vertex(g) ); // Kernel::Point_3(5,2,0)
  vs.push_back( add_vertex(g) ); // Kernel::Point_3(4,0,0)

  CGAL::Euler::add_face(CGAL::make_array(vs[0], vs[1], vs[2]), g);
  CGAL::Euler::add_face(CGAL::make_array(vs[1], vs[3], vs[2]), g);

  // force vertex halfedge to not be a border halfedge
  for(vertex_descriptor v : vertices(g))
  {
    halfedge_descriptor h = halfedge(v, g);
    if ( CGAL::is_border(h, g) )
      set_halfedge(v, prev(opposite(h, g), g), g);
    assert(target(halfedge(v, g), g)==v);
  }

  vs.push_back( add_vertex(g) ); // Kernel::Point_3(0,0,0)
  vs.push_back( add_vertex(g) );  // Kernel::Point_3(1,0,0)
  CGAL::Euler::add_face(CGAL::make_array(vs[4],vs[5],vs[0]), g);

  vs.push_back( add_vertex(g) ); // Kernel::Point_3(2,0,0)
  vs.push_back( add_vertex(g) ); // Kernel::Point_3(3,0,0)
  CGAL::Euler::add_face(CGAL::make_array(vs[6],vs[7],vs[1]), g);
  CGAL::Euler::add_face(CGAL::make_array(vs[7],vs[3],vs[1]), g);
}

template <typename T>
void
add_faces()
{
  typedef typename boost::graph_traits<T>::vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<T>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<T>::halfedge_descriptor halfedge_descriptor;

  // read a mesh with bord + test append
  {
  T m;

  for (int i=0; i<2; ++i)
  {
    std::ifstream in(CGAL::data_file_path("meshes/head.off"));
    std::vector<Kernel::Point_3> points;
    std::vector<std::array<std::size_t, 3> > faces_ids;
    CGAL::IO::read_OFF(in, points, faces_ids);

    std::vector<vertex_descriptor> verts;
    verts.reserve(points.size());
    std::vector<std::array<vertex_descriptor, 3> > faces_vs;
    faces_vs.reserve(faces_ids.size());
    for (std::size_t k=0; k<points.size(); ++k)
      verts.push_back( add_vertex(m) );
    for (const std::array<std::size_t, 3>& f : faces_ids)
      faces_vs.push_back( CGAL::make_array(verts[f[0]], verts[f[1]], verts[f[2]]) );
    CGAL::Euler::add_faces(faces_vs, m);
  }
  assert(faces(m).size()==2*2918);
  assert(vertices(m).size()==2*1487);
  assert( CGAL::is_valid_polygon_mesh(m) );
  }

  // closing a cube no extra vertex
  {
    std::ifstream in(CGAL::data_file_path("meshes/open_cube.off"));
    T m;
    CGAL::IO::read_OFF(in, m);
    std::vector<vertex_descriptor> verts(vertices(m).begin(), vertices(m).end());
    std::list< std::vector<vertex_descriptor> > new_faces;
    new_faces.push_back({verts[1], verts[7], verts[4]});
    new_faces.push_back({verts[1], verts[0], verts[7]});
    CGAL::Euler::add_faces(new_faces, m);
    assert(CGAL::is_closed(m));
    assert(CGAL::is_valid_polygon_mesh(m));
  }
  // closing a cube with extra vertex
  {
    std::ifstream in(CGAL::data_file_path("meshes/open_cube.off"));
    T m;
    CGAL::IO::read_OFF(in, m);
    std::vector<vertex_descriptor> verts(vertices(m).begin(), vertices(m).end());
    verts.push_back(add_vertex(m));
    put(CGAL::vertex_point, m, verts.back(), Kernel::Point_3(50,0,50));
    std::list< std::vector<vertex_descriptor> > new_faces;
    new_faces.push_back({verts[1], verts[0], verts.back()});
    new_faces.push_back({verts[0], verts[7], verts.back()});
    new_faces.push_back({verts[7], verts[4], verts.back()});
    new_faces.push_back({verts[4], verts[1], verts.back()});

    CGAL::Euler::add_faces(new_faces, m);
    assert(CGAL::is_closed(m));
    assert(CGAL::is_valid_polygon_mesh(m));
  }

  // build a model with non-manifold vertices
  for (int run=0; run<4; ++run)
  {
    T m;
    std::vector<Kernel::Point_3> points;
    points.push_back( Kernel::Point_3(0,0,0) );//v0
    points.push_back( Kernel::Point_3(4,0,0) );//v1
    points.push_back( Kernel::Point_3(8,0,0) );//v2
    points.push_back( Kernel::Point_3(4,1,0) );//v3
    points.push_back( Kernel::Point_3(4,2,0) );//v4
    points.push_back( Kernel::Point_3(4,3,0) );//v5
    points.push_back( Kernel::Point_3(0,3,0) );//v6
    points.push_back( Kernel::Point_3(1,0,3) );//v7
    points.push_back( Kernel::Point_3(2,0,3) );//v8
    points.push_back( Kernel::Point_3(3,0,3) );//v9
    points.push_back( Kernel::Point_3(5,0,3) );//v10
    points.push_back( Kernel::Point_3(6,0,3) );//v11
    points.push_back( Kernel::Point_3(7,0,3) );//v12
    points.push_back( Kernel::Point_3(4,-3,0) );//v13

    std::vector<vertex_descriptor> v;
    v.reserve(points.size());

    for (const Kernel::Point_3& p : points)
    {
      v.push_back(add_vertex(m));
      put(CGAL::vertex_point, m, v.back(), p);
    }

    // used only to control which border halfedge is created first
    std::vector<std::array<vertex_descriptor, 3> > face_array;
    face_array.push_back(CGAL::make_array(v[0], v[1], v[3]));
    face_array.push_back(CGAL::make_array(v[2], v[3], v[1]));
    if (run>1)
      std::swap(face_array[0], face_array[1]);

    // add faces
    CGAL::Euler::add_face(face_array[0], m);
    CGAL::Euler::add_face(face_array[1], m);
    face_descriptor to_be_removed1 = CGAL::Euler::add_face(CGAL::make_array(v[0], v[3], v[4]), m);
    face_descriptor to_be_removed2 = CGAL::Euler::add_face(CGAL::make_array(v[2], v[4], v[3]), m);
    CGAL::Euler::add_face(CGAL::make_array(v[0], v[4], v[5], v[6]), m);
    CGAL::Euler::add_face(CGAL::make_array(v[2], v[5], v[4]), m);
    CGAL::Euler::add_face(CGAL::make_array(v[1], v[7], v[8]), m);
    CGAL::Euler::add_face(CGAL::make_array(v[1], v[8], v[9]), m);
    CGAL::Euler::add_face(CGAL::make_array(v[1], v[10], v[11]), m);
    CGAL::Euler::add_face(CGAL::make_array(v[1], v[11], v[12]), m);

    if (run%2==0)
    {
      for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(to_be_removed1, m), m))
        set_face(h, boost::graph_traits<T>::null_face(), m);
      remove_face(to_be_removed1, m);
    }
    else
    {
      for(halfedge_descriptor h : CGAL::halfedges_around_face(halfedge(to_be_removed2, m), m))
        set_face(h, boost::graph_traits<T>::null_face(), m);
      remove_face(to_be_removed2, m);
    }
    std::vector< std::vector<vertex_descriptor> > new_faces;
    new_faces.push_back( {v[0], v[13], v[1]} );
    new_faces.push_back( {v[1], v[13], v[2]} );

    std::vector<halfedge_descriptor> border_hedges;
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m, std::back_inserter(border_hedges));
    assert(border_hedges.size()==2);

    CGAL::Euler::add_faces(new_faces, m);

    border_hedges.clear();
    CGAL::Polygon_mesh_processing::extract_boundary_cycles(m, std::back_inserter(border_hedges));
    assert(border_hedges.size()==3);
  }
}

template <typename Graph>
void
test_Euler_operations()
{
  test_copy_face_graph_nm_umbrella<Graph>();
  test_copy_face_graph_isolated_vertices<Graph>();
  join_face_test<Graph>();
  add_vertex_and_face_to_border_test<Graph>();
  add_face_to_border_test<Graph>();
  join_vertex_interior_test<Graph>();
  join_vertex_exterior_test<Graph>();
  split_vertex<Graph>();
  split_join_vertex_inverse<Graph>();
  join_loop_test<Graph>();
  split_loop_test<Graph>();
  split_face_test<Graph>();
  make_hole_test<Graph>();
  remove_center_vertex_test<Graph>();
  join_split_inverse<Graph>();
  does_satisfy_link_condition<Graph>();
  test_swap_edges<Graph>();
  add_face_bug<Graph>();
  add_faces<Graph>();
}

int main()
{
  test_Euler_operations<Polyhedron>();
  test_Euler_operations<SM>();
  test_Euler_operations<LCC>();

#ifdef CGAL_USE_OPENMESH
  test_Euler_operations<OMesh>();
#endif

  std::cerr << "done\n";
  return 0;
}
