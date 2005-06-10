#ifndef _GRAPH_H
#define _GRAPH_H

#include <CGAL/Node.h>

CGAL_BEGIN_NAMESPACE

template < class E >
class Graph {

 public:
  typedef std::list<E> Elements;
  typedef typename Elements::const_iterator Elements_iterator;
  typedef std::map<E, Node<E>*> Nodes_map;
  typedef typename Nodes_map::iterator Nodes_map_iterator;

 private:
  Nodes_map nodes;

 public:
  // Constructor
  Graph() {
    nodes = Nodes_map();
  }

 public:
  Nodes_map get_nodes() const {
    return nodes;
  }
  
  void add_node(const E& e, const std::list<E>& loe) {
    Node<E>* n = new Node<E>(e);
    nodes.insert( std::make_pair(e, n) );

    // add all the connected components of e
    n = ( *( nodes.find(e) ) ).second;
    for (Elements_iterator it = loe.begin();
	 it != loe.end();
	 ++it) {
      Node<E>* succ = new Node<E>(*it);
      nodes.insert( std::make_pair(*it, succ) );      
      succ = ( *( nodes.find(*it) ) ).second;
      n->add_succ(succ);
    }
  }

  void remove_node(const E& e) {
    Node<E>* n = ( *( nodes.find(e) ) ).second;
    delete(n);

    nodes.erase(e);
  }

  // go all over the connected component where n is
  void go_through_connected_component(Node<E>*& n) {
    if ( ! n->is_visited() ) {
      n->set_visited(true);
      Nodes_map succ = n->get_succ();
      for (Nodes_map_iterator nit = succ.begin();
	   nit != succ.end();
	   ++nit) {
	go_through_connected_component( (*nit).second );
      }
    }
  }

  // go over all the connected components and count them
  int nb_connected_components() {
    int count = 0;
    for (Nodes_map_iterator nit = nodes.begin();
	 nit != nodes.end();
	 ++nit) {
      Node<E>* n = (*nit).second;
      if ( ! n->is_visited() ) {
	count++;
	go_through_connected_component(n);
      }
    }
    for (Nodes_map_iterator nit = nodes.begin();
	 nit != nodes.end();
	 ++nit) {
      Node<E>* n = (*nit).second;
      n->set_visited(false);
    }    
    return count;
  }
  
  int size() const {
    return nodes.size();
  }

}; // end Graph

CGAL_END_NAMESPACE

#endif // GRAPH
