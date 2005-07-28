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
  int nb_nodes;

 public:
  // Constructor
  Graph() {
    nodes = Nodes_map();
    nb_nodes = 0;
  }

 public:
  Nodes_map& get_nodes() {
    return nodes;
  }
  
  void add_node(const E& e, const std::list<E>& loe) {
    ++nb_nodes;
    typename Nodes_map::iterator nit = nodes.find(e);
    Node<E>* n;
    if (  nit == nodes.end()) {
      n = new Node<E>(e);
      nodes.insert( std::make_pair(e, n) );
    }
    else n = nit->second;
    for (Elements_iterator it = loe.begin();
	 it != loe.end();
	 ++it) {
      Node<E>* succ;
      nit = nodes.find(*it);
      if (  nit == nodes.end()) {
	Node<E>* succ = new Node<E>(*it);
	nodes.insert( std::make_pair(*it, succ) ); 
      }
      else succ = nit->second;
      n->add_succ(succ);
    }
  }


  void remove_node(const E& e) {
    --nb_nodes;

    Node<E>* n = ( *( nodes.find(e) ) ).second;
    delete(n);

    nodes.erase(e);
  }

  // travels the connected component of n
  int nb_visited;
  void go_through_connected_component(Node<E>*& n, bool status) {
      CGAL_assertion ( n->is_visited() != status );
      ++nb_visited;
      n->set_visited (status);
      Nodes_map& succ = n->get_succ();
      for (Nodes_map_iterator nit = succ.begin();
	   nit != succ.end();
	   ++nit)
	if ((*nit).second->is_visited() != status )
	  go_through_connected_component( (*nit).second, status );
/*     } */
  }

  // says whether there is only one CC
  bool is_graph_connected () {
    CGAL_assertion (!nodes.empty ());
    Node<E>* n = (*(nodes.begin())).second;

    CGAL_assertion ( !n->is_visited ());
    nb_visited = 0;
    go_through_connected_component(n, true);

    CGAL_assertion (nb_visited > 0);
    CGAL_assertion (nb_visited <= size());
    bool result = (nb_visited == size());

    go_through_connected_component(n, false);
/*     for (Nodes_map_iterator nit = nodes.begin(); */
/* 	 nit != nodes.end(); */
/* 	 ++nit) { */
/*       Node<E>* n = (*nit).second; */
/*       n->set_visited(false); */
/*     } */

    return result;
  }
  
  int size() const {
/*     return nodes.size(); */
    return nb_nodes;
  }

}; // end Graph

CGAL_END_NAMESPACE

#endif // GRAPH
