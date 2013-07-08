#include <iostream>
#include <fstream>
#include <map>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

// Extended polyhedron items which include an id() field
#include <CGAL/Polyhedron_items_with_id_3.h>

// Adaptor for Polyhedron_3
#include <CGAL/Surface_mesh_simplification/HalfedgeGraph_Polyhedron_3.h>

// Simplification function
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

// Stop-condition policy
#include <CGAL/internal/Mean_curvature_skeleton/Edge_minimum_length_stop_predicate.h>

// Non-default cost and placement policies
#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Midpoint_and_length.h>

// Map used to mark edges as fixed
#include <CGAL/Unique_hash_map.h>

typedef CGAL::Simple_cartesian<double> Kernel;
//
// Setup an enriched polyhedron type which stores an id() field in the items
//
typedef CGAL::Polyhedron_3<Kernel,CGAL::Polyhedron_items_with_id_3> Surface;

typedef boost::graph_traits<Surface>::vertex_descriptor	        vertex_descriptor;
typedef boost::graph_traits<Surface>::vertex_iterator            vertex_iterator;
typedef boost::graph_traits<Surface>::edge_descriptor            edge_descriptor;
typedef boost::graph_traits<Surface>::edge_iterator              edge_iterator;

namespace SMS = CGAL::Surface_mesh_simplification ;

//
// BGL property map which indicates whether an edge is border OR is marked as non-removable
//
class Constrains_map : public boost::put_get_helper<bool,Constrains_map>
{
public:

  typedef boost::readable_property_map_tag                    category;
  typedef bool                                                value_type;
  typedef bool                                                reference;
  typedef boost::graph_traits<Surface const>::edge_descriptor key_type;

  Constrains_map() : mConstrains(false) {}

  reference operator[](key_type const& e) const { return e->is_border() || is_constrained(e) ; }
  
  void set_is_constrained ( key_type const& e, bool is )  { mConstrains[e]=is; }
  
  bool is_constrained( key_type const& e ) const { return mConstrains.is_defined(e) ? mConstrains[e] : false ; }
  
private:
  
  CGAL::Unique_hash_map<key_type,bool> mConstrains ;
  
};

int main( int argc, char** argv ) 
{
  Surface surface; 
  
  std::ifstream is(argv[1]) ; is >> surface ;

  int edge_id = 0;
  edge_iterator eb, ee;
  for (boost::tie(eb, ee) = boost::edges(surface); eb != ee; ++eb)
  {
    eb->id() = edge_id++;
  }

  int vertex_id = 0;
  vertex_iterator vb, ve;
  for (boost::tie(vb, ve) = boost::vertices(surface); vb != ve; ++vb)
  {
    vb->id() = vertex_id++;
  }

  // This is a stop predicate (defines when the algorithm terminates).
  // The simplification stops when the length of all edges is greater than the minimum threshold.
  CGAL::internal::Minimum_length_predicate<Surface> stop(0.010);

  std::map<size_t, bool> is_vertex_fixed_map;
  is_vertex_fixed_map.clear();

  for (boost::tie(vb, ve) = boost::vertices(surface); vb != ve; vb++)
  {
    int r = rand() % 10000;
    if (r < 100)
    {
      vertex_descriptor vi = *vb;
      is_vertex_fixed_map[vi->id()] = true;
    }
  }

  Constrains_map constrains_map ;
     
  for (boost::tie(eb, ee) = boost::edges(surface); eb != ee; ++eb)
  {
    vertex_descriptor vi = boost::source(*eb, surface);
    vertex_descriptor vj = boost::target(*eb, surface);
    size_t vi_idx = vi->id();
    size_t vj_idx = vj->id();

    if (is_vertex_fixed_map.find(vi_idx) == is_vertex_fixed_map.end()
     || is_vertex_fixed_map.find(vj_idx) == is_vertex_fixed_map.end())
    {
      continue;
    }
    if (is_vertex_fixed_map[vi_idx] || is_vertex_fixed_map[vj_idx])
    {
      constrains_map.set_is_constrained(*eb, true);
    }
  }

  int cnt = 0;
  for (boost::tie(vb, ve) = boost::vertices(surface); vb != ve; ++vb)
  {
    int id = vb->id();
    if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
    {
      if (is_vertex_fixed_map[id])
      {
        cnt++;
      }
    }
  }
  std::cerr << "before collapse " << cnt << " fixed vertices\n";

  // This the actual call to the simplification algorithm.
  // The surface and stop conditions are mandatory arguments.
  // The index maps are needed because the vertices and edges
  // of this surface lack an "id()" field.
  int r = SMS::edge_collapse
            (surface
            ,stop
            ,CGAL::get_cost(SMS::Edge_length_cost<Surface>())
                      .get_placement(SMS::Midpoint_placement<Surface>())
                      .edge_is_border_map(constrains_map)
            );
  
  cnt = 0;
  for (boost::tie(vb, ve) = boost::vertices(surface); vb != ve; ++vb)
  {
    int id = vb->id();
    if (is_vertex_fixed_map.find(id) != is_vertex_fixed_map.end())
    {
      if (is_vertex_fixed_map[id])
      {
        cnt++;
      }
    }
  }
  std::cerr << "after collapse " << cnt << " fixed vertices\n";

  std::cout << "\nFinished...\n" << r << " edges removed.\n" 
            << (surface.size_of_halfedges()/2) << " final edges.\n" ;
        
  std::ofstream os( argc > 2 ? argv[2] : "out.off" ) ; os << surface ;
  
  return 0 ;      
}

