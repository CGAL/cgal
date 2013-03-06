#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/boost/graph/properties_Polyhedron_3.h>
#include <CGAL/boost/graph/halfedge_graph_traits_Polyhedron_3.h>

/** 
 * Provide simple functionality for iterating over spoke and rim edges
 *   - use get_descriptor() to obtain active edge
 *   - get_iterator() always holds spoke edges */
 /// \code
 /// // how to use Spokes_and_rims_iterator
 /// boost::tie(e_begin, e_end) = boost::out_edges(vertex, polyhedron);
 /// Spokes_and_rims_iterator<Polyhedron> rims_it(e_begin, polyhedron);
 /// 
 /// for ( ; rims_it.get_iterator() != e_end; ++rims_it )
 /// {
 ///   edge_descriptor active_edge = rims_it.get_descriptor();
 ///   // use active_edge as you like
 /// }
 /// \endcode
template<class Polyhedron>
class Spokes_and_rims_iterator
{
public:
  typedef typename boost::graph_traits<Polyhedron>::out_edge_iterator	out_edge_iterator;
  typedef typename boost::graph_traits<Polyhedron>::edge_descriptor		edge_descriptor; 

  Spokes_and_rims_iterator(out_edge_iterator edge_iterator, Polyhedron& polyhedron)
    : iterator(edge_iterator), descriptor(*edge_iterator), polyhedron(polyhedron), is_current_rim(false)
  { }

  /// descriptor will be assigned to next valid edge, note that iterator might not change
  Spokes_and_rims_iterator<Polyhedron>&
  operator++() 
  {
    // loop through one spoke then one rim edge
    if(!is_current_rim && !boost::get(CGAL::edge_is_border, polyhedron, descriptor)) // it is rim edge's turn
    {
	    is_current_rim = true;
	    descriptor = CGAL::next_edge(descriptor, polyhedron);
    }
    else // if current edge is rim OR there is no rim edge (current spoke edge is boudary)
    {    // then iterate to next spoke edge
	    is_current_rim = false;
	    descriptor = *(++iterator);
    }
    return *this;
  }

  out_edge_iterator get_iterator()   { return iterator; }
  edge_descriptor   get_descriptor() { return descriptor; }

private:
  bool is_current_rim;        ///< current descriptor is rim or spoke
  out_edge_iterator iterator; ///< holds spoke edges (i.e. descriptor is not always = *iterator)
  edge_descriptor descriptor; ///< current active edge descriptor for looping
  Polyhedron& polyhedron;
};