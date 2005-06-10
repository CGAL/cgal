#ifndef _NODE_H
#define _NODE_H

#include <list>
#include <map>

CGAL_BEGIN_NAMESPACE

template < class E >
class Node {

 public:
  typedef std::list<E> Elements;
  typedef std::map<E, Node<E>*> Nodes_map;
  typedef typename Nodes_map::iterator Nodes_map_iterator;

 private:
  bool visited;
  E element;
  Nodes_map succ;
  
 public:
  // Constructor
  Node(const E& e) {
    visited = false;
    element = e;
    succ = Nodes_map();
  }
  
  ~Node() {
    for (Nodes_map_iterator nit = succ.begin();
	 nit != succ.end();
	 ++ nit) {
      Node<E>* succ_node = (*nit).second;
      succ_node->get_succ().erase(element);
    }
  }

 public:
  bool is_visited() const {
    return visited;
  }

  void set_visited(const bool b) {
    visited = b;
  }

  E get_element() const {
    return element;
  }

  void set_element(const E& e) {
    element = e;
  }

  std::map<E, Node<E>*>& get_succ() {
    return succ;
  }

  void add_succ(Node<E>*& n) {
    E e = n->get_element();
    succ.insert( std::make_pair(e, n) );
    n->succ.insert( std::make_pair(element, this) );
  }

  void remove_succ(Node<E>*& n) {
    E e = n->get_element();
    succ.erase(e);
    n->succ.erase(element);
  }

}; // end Node

CGAL_END_NAMESPACE

#endif // NODE
