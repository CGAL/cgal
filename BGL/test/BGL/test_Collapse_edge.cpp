#include "test_Prefix.h"

#include <CGAL/boost/graph/Euler_operations.h>
#include <CGAL/boost/graph/IO/OFF.h>

#include <boost/range/distance.hpp>

#include <string>

template < typename Mesh>
typename boost::graph_traits<Mesh>::
halfedge_descriptor find_halfedge(double x1, double y1,
                                  double x2, double y2,
                                  Mesh& m,
                                  bool is_border = false)
{
  typedef typename boost::property_map<Mesh, CGAL::vertex_point_t>::type VPMAP;
  typedef typename boost::property_traits<VPMAP>::value_type Point;

  typedef typename boost::graph_traits<Mesh>::halfedge_descriptor halfedge_descriptor;
  VPMAP vpmap = get(CGAL::vertex_point, m);
  for(halfedge_descriptor h : halfedges(m))
  {
    if(get(vpmap, source(h, m)) == Point(x1,y1,0)
       && get(vpmap, target(h, m)) == Point(x2,y2,0))
    {
      if(is_border == CGAL::is_border(h, m))
        return h;
      else
        return opposite(h, m);
    }
  }
  return boost::graph_traits<Mesh>::null_halfedge();
}

template <typename Mesh>
void
collapse_edge_test()
{
  CGAL_GRAPH_TRAITS_MEMBERS(Mesh);
  typedef typename boost::graph_traits<Mesh>:: vertex_descriptor vertex_descriptor;
  typedef typename boost::graph_traits<Mesh>:: halfedge_descriptor halfedge_descriptor;

  const std::string fname = "data/flat_hexahedron.off";
  Mesh m;
  if(!CGAL::IO::read_OFF(fname, m))
    std::cout << "Error reading file: " << fname << std::endl;

  bool m_is_valid = CGAL::is_valid(m);
  assert(m_is_valid);

  Mesh test_mesh;
  CGAL::copy_face_graph(m, test_mesh);
  m_is_valid = CGAL::is_valid(m);
  assert(m_is_valid);

  //case 1: General Case.
  {
    halfedge_descriptor he = find_halfedge(-0.5,0,
                                           0.5,0,
                                           test_mesh);
    halfedge_descriptor en = next(he, test_mesh);
    halfedge_descriptor eno = opposite(en, test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(he, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(he, test_mesh);
    bool ok = CGAL::Euler::collapse_edge(edge(he, test_mesh), test_mesh) == v1;
    assert(ok);
    char found = 0;
    for(halfedge_descriptor it : CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == eno
         || it == eno_prime){
        ++found;
      }
    }
    assert(found == 2);
    CGAL::clear(test_mesh);

  }
  //case 2:  collapsing edge is not itself a border, but is incident upon a border edge that is removed.
  {
    CGAL::copy_face_graph(m, test_mesh);
    halfedge_descriptor he = find_halfedge(0,0.5,
                                           -0.75,0.5,
                                           test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);

    he = find_halfedge(-0.5,0,
                       0.5,0,
                       test_mesh);
    halfedge_descriptor en = next(he, test_mesh);
    halfedge_descriptor eno = opposite(en, test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(he, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(he, test_mesh);
    bool ok = CGAL::Euler::collapse_edge(edge(he, test_mesh), test_mesh) == v1;
    assert(ok);
    char found = 0;
    for(halfedge_descriptor it : CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == eno
         || it == eno_prime){
        ++found;
      }
    }
    assert(found == 2);
    CGAL::clear(test_mesh);
  }
  //case 3: collapsing edge is not itself a border, but is incident upon a border edge that is not removed
  {
    CGAL::copy_face_graph(m, test_mesh);
    halfedge_descriptor he = find_halfedge(1.5,0,
                                           0.75,0.5,
                                           test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);

    he = find_halfedge(-0.5,0,
                       0.5,0,
                       test_mesh);
    halfedge_descriptor en = next(he, test_mesh);
    halfedge_descriptor eno = opposite(en, test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(he, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(he, test_mesh);
    bool ok = CGAL::Euler::collapse_edge(edge(he, test_mesh), test_mesh) == v1;
    assert(ok);
    char found = 0;
    for(halfedge_descriptor it : CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == eno
         || it == eno_prime){
        ++found;
      }
    }
    assert(found == 2);
    CGAL::clear(test_mesh);
  }
  //case 4: collapsing edge is itself a border
  {
    CGAL::copy_face_graph(m, test_mesh);
    halfedge_descriptor he = find_halfedge(-0.5, 0,
                                           0, -0.5,
                                           test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);
    he = find_halfedge(0, -0.5,
                       -0.5, 0,
                       test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);
    he = find_halfedge(0, -0.5,
                       0.75, -0.5,
                       test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);


    he = find_halfedge(-0.5,0,
                       0.5,0,
                       test_mesh);
    halfedge_descriptor en = next(he, test_mesh);
    halfedge_descriptor eno = opposite(en, test_mesh);
    halfedge_descriptor ep_prime = prev(opposite(he, test_mesh), test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(he, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(he, test_mesh);
    bool ok = CGAL::Euler::collapse_edge(edge(he, test_mesh), test_mesh) == v1;
    assert(ok);
    char found = 0;
    for(halfedge_descriptor it : CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == eno
         || it == eno_prime
         || it == ep_prime){
        ++found;
      }
    }
    assert(found == 3);
    CGAL::clear(test_mesh);
  }
  //case 5 singular case.
  {
    CGAL::copy_face_graph(m, test_mesh);
    halfedge_descriptor he = find_halfedge(0.75,0.5,
                                           1.5,0,
                                           test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);
    he = find_halfedge(0.75,-0.5,
                       1.5,0,
                       test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);
    he = find_halfedge(0,0.5,
                       0.5,0,
                       test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);
    he = find_halfedge(0.5,0,
                       0,-0.5,
                       test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);

    he = find_halfedge(-0.5,0,
                       0.5,0,
                       test_mesh);
    CGAL::Euler::remove_face(he, test_mesh);
    halfedge_descriptor ep = prev(he, test_mesh);
    halfedge_descriptor eno_prime = opposite(next(opposite(he, test_mesh), test_mesh), test_mesh);
    vertex_descriptor v1 = target(he, test_mesh);
    bool ok = CGAL::Euler::collapse_edge(edge(he, test_mesh), test_mesh) == v1;
    assert(ok);
    char found = 0;
    for(halfedge_descriptor it : CGAL::halfedges_around_target(v1,test_mesh))
    {
      if(it == ep)
        ++found;
      else if( it == eno_prime){
        ++found;
      }
    }
    assert(found == 2);
    CGAL::clear(test_mesh);
  }
}


int main()
{

  collapse_edge_test<Polyhedron>();
  collapse_edge_test<SM>();

#ifdef CGAL_USE_OPENMESH
  collapse_edge_test<OMesh>();
#endif

  std::cerr << "done\n";
  return 0;
}
