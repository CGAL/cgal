// ======================================================================
//
// Copyright (c) 2002 The CGAL Consortium
//
// This software and related documentation is part of an INTERNAL release
// of the Computational Geometry Algorithms Library (CGAL). It is not
// intended for general use.
//
// ----------------------------------------------------------------------
//
// release       :
// release_date  :
//
// file          : include/CGAL/Kd_tree.h
// package       : ASPAS
// revision      : 2.4 
// revision_date : 2003/02/01 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_KD_TREE_H
#define CGAL_KD_TREE_H
#include <CGAL/Kd_tree_node.h>
#include <cassert>
#include <CGAL/Compact_container.h>

namespace CGAL {


template <class Traits> 
class Kd_tree {
public:
  
  typedef typename Traits::Item Item;
  typedef typename std::list<Item>::iterator input_iterator;
  typedef typename Traits::NT NT;
  typedef Kd_tree_node<Traits> Node;
  typedef Kd_tree<Traits> Tree;

  typedef typename Compact_container<Node>::iterator Node_handle;
  typedef std::vector<Item*>::iterator Item_iterator;
private:

  Compact_container<Node> nodes;

  Node_handle tree_root;

  Kd_tree_rectangle<NT>* bbox;
  std::list<Item> pts;

  // Instead of storing the points in arrays in the Kd_tree_node
  // we put all the data in a vector in the Kd_tree.
  // and we only store an iterator range in the Kd_tree_node.
  // 
  std::vector<Item*> data;
  Item_iterator data_iterator;
  Traits tr;
  int the_item_number;

  // protected copy constructor
  Kd_tree(const Tree& tree) {};


  // Instead of the recursive construction of the tree in the class Kd_tree_node
  // we do this in the tree class. The advantage is that we then can optimize
  // the allocation of the nodes.

  // The leaf node
  Node_handle 
  create_leaf_node(Point_container<Item>& c)
  {
    Node n;
    Node_handle nh = nodes.insert(n);
    nh->n = c.size();
    nh->the_node_type = Node::LEAF;
    if (c.size()>0) { 
      nh->data = data_iterator;
      data_iterator = std::copy(c.begin(), c.end(), data_iterator);
    }
    return nh;
  }


  // The internal node

  // TODO: Similiar to the leaf_init function above, a part of the code should be
  //       moved to a the class Kd_tree_node.
  //       It is not proper yet, but the goal was to see if there is
  //       a potential performance gain through the Compact_container
  Node_handle 
  create_internal_node_use_extension(Point_container<Item>& c, Traits& t) 
  {
    Node n;
    Node_handle nh = nodes.insert(n);

    nh->the_node_type = Node::EXTENDED_INTERNAL;

    Point_container<Item> 
      c_low = Point_container<Item>(c.dimension());
    Kd_tree_rectangle<NT> bbox(c.bounding_box());

    t.split(nh->sep, c, c_low);
	        
    int cd  = nh->sep.cutting_dimension();

    nh->low_val = bbox.min_coord(cd);
    nh->high_val = bbox.max_coord(cd);

    if (c_low.size() > t.bucket_size())
      nh->lower_ch = create_internal_node_use_extension(c_low,t);
    else
      nh->lower_ch = create_leaf_node(c_low);

    if (c.size() > t.bucket_size())
      nh->upper_ch = create_internal_node_use_extension(c,t);
    else
      nh->upper_ch = create_leaf_node(c);

    return nh;
  }

  
  // Note also that I duplicated the code to get rid if the if's for
  // the boolean use_extension which was constant over the construction
  Node_handle 
  create_internal_node(Point_container<Item>& c, Traits& t) 
  {
    Node n;
    Node_handle nh = nodes.insert(n);

    nh->the_node_type = Node::INTERNAL;

    Point_container<Item> 
      c_low = Point_container<Item>(c.dimension());
    Kd_tree_rectangle<NT> bbox(c.bounding_box());

    t.split(nh->sep, c, c_low);
	        
    if (c_low.size() > t.bucket_size())
      nh->lower_ch = create_internal_node(c_low,t);
    else
      nh->lower_ch = create_leaf_node(c_low);

    if (c.size() > t.bucket_size())
      nh->upper_ch = create_internal_node(c,t);
    else
      nh->upper_ch = create_leaf_node(c);

    return nh;
  }
  




public:

  
  Kd_tree(input_iterator first, input_iterator beyond,
	    Traits t = Traits()) : tr(t) {
    assert(first != beyond);
    int dim = first->dimension();
    
    std::copy(first, beyond, std::back_inserter(pts));

    data = std::vector<Item*>(pts.size()); // guarantees that iterators we store in Kd_tree_nodes stay valid
    data_iterator = data.begin();

    Point_container<Item> c(dim, pts.begin(), pts.end());

    bbox = new Kd_tree_rectangle<NT>(c.bounding_box());
    
    the_item_number=c.size();
    if (c.size() <= t.bucket_size())
      tree_root = create_leaf_node(c);
    else 
		if (t.use_extended_nodes())
		tree_root = create_internal_node_use_extension(c,t);  
		else
		tree_root = create_internal_node(c,t); 
	
  }

  template <class OutputIterator, class Rectangle>
	OutputIterator search(OutputIterator it, Rectangle& r, NT eps=NT(0)) {
                //  Why do you create this on the heap?
		Kd_tree_rectangle<NT>* b = new Kd_tree_rectangle<NT>(*bbox);
		it=tree_root->tree_items_in_rectangle(it,r,b,eps);
		delete b;
		return it;
	}

  template <class OutputIterator>
	OutputIterator search(OutputIterator it, Item& center, NT radius, NT eps=NT(0)) {
		Kd_tree_rectangle<NT>* b = new Kd_tree_rectangle<NT>(*bbox);
		it=tree_root->tree_items_in_sphere(it, center,
		(radius-eps)*(radius-eps), radius*radius,
                (radius+eps)*(radius+eps), b);
		delete b;
		return it;
	}

  template <class OutputIterator>
	OutputIterator report_all_points(OutputIterator it) 
	{it=tree_root->tree_items(it);
	 return it;}

    ~Kd_tree() {
		  delete bbox;
	};


  Traits traits() const {return tr;} // Returns the traits class;

  Node_handle root() { return tree_root; }

  Kd_tree_rectangle<NT>* bounding_box() {return bbox; }

  int item_number() {return the_item_number;}

  // Print statistics of the tree.
  void statistics () {
    std::cout << "Tree statistics:" << std::endl;
    std::cout << "Number of items stored: " 
		  << tree_root->num_items() << std::endl;
    std::cout << " Tree depth: " << tree_root->depth() << std::endl;
  }


};

} // namespace CGAL
#endif // CGAL_KD_TREE_H
