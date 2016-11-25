#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Mesh_3/dihedral_angle_3.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>

#include <iostream>
#include <fstream>

namespace PMP = CGAL::Polygon_mesh_processing;
namespace params = PMP::parameters;

template <class TriangleMesh, class VertexPointMap, class EdgeIsConstrainedMap>
std::size_t
mark_sharp_edge(const TriangleMesh& tm,
                      VertexPointMap vpm,
                      EdgeIsConstrainedMap ecm,
                      double angle_threshold = 100)
{
  typedef boost::graph_traits<TriangleMesh> GT;
  typedef typename GT::edge_descriptor edge_descriptor;
  typedef typename GT::halfedge_descriptor halfedge_descriptor;
  std::size_t nb_sharp_edges = 0;

  BOOST_FOREACH(edge_descriptor ed, edges(tm))
  {
    halfedge_descriptor hd = halfedge(ed, tm);

    if ( is_border_edge(hd, tm) )
    {
      ++nb_sharp_edges;
      put(ecm, ed, true);
    }
    else{
      double angle = CGAL::Mesh_3::dihedral_angle(get(vpm, source(hd, tm)),
                                                  get(vpm, target(hd, tm)),
                                                  get(vpm, target(next(hd, tm), tm)),
                                                  get(vpm, target(next(opposite(hd, tm), tm), tm)));
      if ( CGAL::abs(angle) < angle_threshold )
      {
        ++nb_sharp_edges;
        put(ecm, ed, true);
      }
    }
  }
}


typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Surface_mesh<K::Point_3> Triangle_mesh;
typedef Triangle_mesh::Property_map<Triangle_mesh::Edge_index,bool> Constrained_edge_map;

void read_input(Triangle_mesh& tm)
{
  std::ifstream input("data/joint_refined.off");
  input >> tm;
}

void translate(Triangle_mesh& tm,
               K::Vector_3 tslt = K::Vector_3(0.2, 0.2, 0.2))
{
  BOOST_FOREACH(Triangle_mesh::Vertex_index vi, tm.vertices())
    tm.point(vi) = tm.point(vi)+tslt;
}

void dump_constrained_edges(const Triangle_mesh& tm, const Constrained_edge_map& ecm, const char* fname)
{
  std::ofstream output(fname);
  BOOST_FOREACH(Triangle_mesh::Edge_index e, tm.edges())
  {
    if ( get(ecm, e) )
      output << "2 " << tm.point( tm.vertex(e, 0) )
             <<  " " << tm.point( tm.vertex(e, 1) ) << "\n";
  }
}

std::size_t
count_constrained_edges(const Triangle_mesh& tm, const Constrained_edge_map& ecm)
{
  std::size_t n=0;
  BOOST_FOREACH(Triangle_mesh::Edge_index e, tm.edges())
  {
    if ( get(ecm, e) ) ++n;
  }
  return n;
}

void test_corefine(Triangle_mesh tm1, Triangle_mesh tm2)
{
  Constrained_edge_map ecm1 =
    tm1.property_map<Triangle_mesh::Edge_index,bool>("e:cst").first;
  Constrained_edge_map ecm2 =
    tm2.property_map<Triangle_mesh::Edge_index,bool>("e:cst").first;

  assert( count_constrained_edges(tm1, ecm1)==307 );
  assert( count_constrained_edges(tm2, ecm2)==307 );

  PMP::corefine(tm1,
                tm2,
                params::edge_is_constrained_map(ecm1),
                params::edge_is_constrained_map(ecm2) );

  assert( count_constrained_edges(tm1, ecm1)==658 );
  assert( count_constrained_edges(tm2, ecm2)==655 );
}

void test_union(Triangle_mesh tm1, Triangle_mesh tm2, Triangle_mesh tm_out,
                const char* outname, bool skip_test_1=false, bool skip_test_2=false)
{
  Constrained_edge_map ecm1 =
    tm1.property_map<Triangle_mesh::Edge_index,bool>("e:cst").first;
  Constrained_edge_map ecm2 =
    tm2.property_map<Triangle_mesh::Edge_index,bool>("e:cst").first;
  Constrained_edge_map ecm_out =
    tm_out.property_map<Triangle_mesh::Edge_index,bool>(outname).first;

  assert( count_constrained_edges(tm1, ecm1)==307 );
  assert( count_constrained_edges(tm2, ecm2)==307 );

  PMP::corefine_and_compute_union(tm1,
                                  tm2,
                                  tm_out,
                                  params::edge_is_constrained_map(ecm1),
                                  params::edge_is_constrained_map(ecm2),
                                  params::edge_is_constrained_map(ecm_out) );

  assert( skip_test_1 || count_constrained_edges(tm1, ecm1)==658 );
  assert( skip_test_2 || count_constrained_edges(tm2, ecm2)==655 );

  dump_constrained_edges(tm1, ecm1, "out_cst1.cgal");
  dump_constrained_edges(tm2, ecm2, "out_cst2.cgal");
  dump_constrained_edges(tm_out, ecm_out, "out_cst_out.cgal");

  std::cout << "nb cst " << count_constrained_edges(tm_out, ecm_out) << "\n";
  std::ofstream output("union.off");
  output << tm_out;
  output.close();

  assert( count_constrained_edges(tm2, ecm2)==978 );
}

int main()
{
  Triangle_mesh tm1, tm2;
  read_input(tm1);
  read_input(tm2);
  translate(tm2);

  Constrained_edge_map ecm1 =
    tm1.add_property_map<Triangle_mesh::Edge_index,bool>("e:cst", false).first;
  Constrained_edge_map ecm1_out =
    tm1.add_property_map<Triangle_mesh::Edge_index,bool>("e:cst_out", false).first;
  Constrained_edge_map ecm2 =
    tm2.add_property_map<Triangle_mesh::Edge_index,bool>("e:cst", false).first;
  Constrained_edge_map ecm2_out =
    tm2.add_property_map<Triangle_mesh::Edge_index,bool>("e:cst_out", false).first;

  mark_sharp_edge(tm1, get(CGAL::vertex_point, tm1), ecm1);
  mark_sharp_edge(tm1, get(CGAL::vertex_point, tm1), ecm1_out);
  mark_sharp_edge(tm2, get(CGAL::vertex_point, tm2), ecm2);
  mark_sharp_edge(tm2, get(CGAL::vertex_point, tm2), ecm2_out);

  assert( count_constrained_edges(tm1, ecm1)==307 );
  assert( count_constrained_edges(tm2, ecm2)==307 );

  Triangle_mesh tm_out;
  tm_out.add_property_map<Triangle_mesh::Edge_index,bool>("e:cst_out", false);

  std::cout << "Testing corefinement\n";
  test_corefine(tm1, tm2);
  std::cout << "Testing union out-of-place\n";
  test_union(tm1, tm2, tm_out, "e:cst_out");
  std::cout << "Testing union in-place\n";
  test_union(tm1, tm2, tm1, "e:cst_out");
  test_union(tm1, tm2, tm2, "e:cst_out");
  test_union(tm1, tm2, tm1, "e:cst", true, false);
  test_union(tm1, tm2, tm2, "e:cst", false, true);
}
