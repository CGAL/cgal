#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <iostream>
#include <fstream>

#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>

#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_items_with_id_3.h>
#include <CGAL/Simple_cartesian.h>

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Visitor base
#include <CGAL/Surface_mesh_simplification/Edge_collapse_visitor_base.h>

// Extended polyhedron items which include an id() field
#include <CGAL/Polyhedron_items_with_id_3.h>

// Stop-condition policy
#include <CGAL/internal/Mean_curvature_skeleton/Edge_minimum_length_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

typedef CGAL::Simple_cartesian<double>     Kernel;
typedef Kernel::Point_3                    Point_3;
typedef CGAL::Polyhedron_3<Kernel, CGAL::Polyhedron_items_with_id_3>  Polyhedron;

typedef Polyhedron::Vertex_iterator        Vertex_iterator;

template<class PolyhedronWithId, class KeyType>
struct Polyhedron_with_id_property_map
    : public boost::put_get_helper<std::size_t&,
             Polyhedron_with_id_property_map<PolyhedronWithId, KeyType> >
{
public:
    typedef KeyType      key_type;
    typedef std::size_t  value_type;
    typedef value_type&  reference;
    typedef boost::lvalue_property_map_tag category;

    reference operator[](key_type key) const { return key->id(); }
};

typedef boost::graph_traits<Polyhedron>::vertex_descriptor             vertex_descriptor;
typedef boost::graph_traits<Polyhedron>::vertex_iterator               vertex_iterator;
typedef boost::graph_traits<Polyhedron>::edge_descriptor               edge_descriptor;
typedef boost::graph_traits<Polyhedron>::edge_iterator                 edge_iterator;

typedef Polyhedron_with_id_property_map<Polyhedron, vertex_descriptor> Vertex_index_map; // use id field of vertices
typedef Polyhedron_with_id_property_map<Polyhedron, edge_descriptor>   Edge_index_map;   // use id field of edges

namespace SMS = CGAL::Surface_mesh_simplification;

int main()
{
    Polyhedron P;
    std::ifstream f("test_collapse.off", std::ios::in) ;
    f >> P;

    Vertex_index_map vertex_id_pmap;
    Edge_index_map edge_id_pmap;

    // initialize index maps
    vertex_iterator vb, ve;
    int vertex_id_count = 0;
    for (boost::tie(vb, ve) = boost::vertices(P); vb != ve; ++vb)
    {
      boost::put(vertex_id_pmap, *vb, vertex_id_count++);
    }

    edge_iterator eb, ee;
    int idx = 0;
    for (boost::tie(eb, ee) = boost::edges(P); eb != ee; ++eb)
    {
      boost::put(edge_id_pmap, *eb, idx++);
    }

    double edgelength_TH = 0.00236;
    // This is a stop predicate (defines when the algorithm terminates).
    // The simplification stops when the length of all edges is greater than the minimum threshold.
    CGAL::internal::Minimum_length_predicate<Polyhedron> stop(edgelength_TH);

    std::cout << "edge length threshold " << edgelength_TH << "\n";

    int cnt = 0;
    for (boost::tie(eb, ee) = boost::edges(P); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, P);
      vertex_descriptor vj = boost::target(*eb, P);
      Point_3 pi = vi->point();
      Point_3 pj = vj->point();
      double dis = sqrt(squared_distance(pi, pj));
      if (dis < edgelength_TH)
      {
        //std::cout << "edge length still shorter than threshold: " << dis << "\n";
        cnt++;
      }
    }
    std::cout << "before collapse " << cnt << " edges shorter than threshold\n";

    int r = SMS::edge_collapse
                (P
                ,stop
                ,CGAL::get_cost     (SMS::Edge_length_cost  <Polyhedron>())
                      .get_placement(SMS::Midpoint_placement<Polyhedron>())
                );

    cnt = 0;
    for (boost::tie(eb, ee) = boost::edges(P); eb != ee; ++eb)
    {
      vertex_descriptor vi = boost::source(*eb, P);
      vertex_descriptor vj = boost::target(*eb, P);
      Point_3 pi = vi->point();
      Point_3 pj = vj->point();
      double dis = sqrt(squared_distance(pi, pj));
      if (dis < edgelength_TH)
      {
        //std::cout << "edge length still shorter than threshold: " << dis << "\n";
        cnt++;
      }
    }
    std::cout << "after collapse " << cnt << " edges shorter than threshold\n";

    return 0;
}

