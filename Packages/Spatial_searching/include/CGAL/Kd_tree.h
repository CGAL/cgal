// Copyright (c) 2002  Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_KD_TREE_H
#define CGAL_KD_TREE_H
#include <CGAL/basic.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/point_generators_2.h>
#include <CGAL/algorithm.h>
#include <CGAL/Orthogonal_standard_search.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Kd_tree_node.h>
#include <cassert>
#include <CGAL/Compact_container.h>

namespace CGAL {


template <class TreeTraits> 
class Kd_tree {
public:
  
  typedef typename TreeTraits::Point Point;
  typedef typename TreeTraits::Container Point_container;
  typedef typename std::list<Point>::iterator input_iterator;
  typedef typename TreeTraits::NT NT;
  typedef Kd_tree_node<TreeTraits> Node;
  typedef Kd_tree<TreeTraits> Tree;

  typedef typename Compact_container<Node>::iterator Node_handle;
  typedef typename std::vector<Point*>::iterator Point_iterator;


private:

  Compact_container<Node> nodes;

  Node_handle tree_root;

  Kd_tree_rectangle<NT>* bbox;
  std::list<Point> pts;

  // Instead of storing the points in arrays in the Kd_tree_node
  // we put all the data in a vector in the Kd_tree.
  // and we only store an iterator range in the Kd_tree_node.
  // 
  std::vector<Point*> data;
  Point_iterator data_iterator;
  TreeTraits tr;
  int the_item_number;

  // protected copy constructor
  Kd_tree(const Tree& tree) {};


  // Instead of the recursive construction of the tree in the class Kd_tree_node
  // we do this in the tree class. The advantage is that we then can optimize
  // the allocation of the nodes.

  // The leaf node
  Node_handle 
  create_leaf_node(Point_container& c)
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
  create_internal_node_use_extension(Point_container& c, TreeTraits& t) 
  {
    Node n;
    Node_handle nh = nodes.insert(n);
    
    nh->the_node_type = Node::EXTENDED_INTERNAL;

    Point_container
      c_low = Point_container(c.dimension());
    
    t.split(nh->sep, c, c_low);
	        
    int cd  = nh->sep.cutting_dimension();
    
    nh->low_val = c_low.bounding_box().min_coord(cd);
   

    
    nh->high_val = c.bounding_box().max_coord(cd);

    
    assert(nh->sep.cutting_value() >= nh->low_val);
    assert(nh->sep.cutting_value() <= nh->high_val);

    

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
  create_internal_node(Point_container& c, TreeTraits& t) 
  {
    Node n;
    Node_handle nh = nodes.insert(n);
    
    nh->the_node_type = Node::INTERNAL;

    Point_container
    c_low = Point_container(c.dimension());
    
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

  //introduced for backward compability
  Kd_tree() {}
  
  Kd_tree(input_iterator first, input_iterator beyond,
	    TreeTraits t = TreeTraits()) : tr(t) {
    assert(first != beyond);
    int dim = first->dimension();
    
    std::copy(first, beyond, std::back_inserter(pts));
 
    data = std::vector<Point*>(pts.size()); // guarantees that iterators we store in Kd_tree_nodes stay valid
    data_iterator = data.begin();

    Point_container c(dim, pts.begin(), pts.end());

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

 
  template <class OutputIterator, class FuzzyQueryItem>
	OutputIterator search(OutputIterator it, const FuzzyQueryItem& q) {
		Kd_tree_rectangle<NT>* b = new Kd_tree_rectangle<NT>(*bbox);
		tree_root->search(it,q,b);
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


  TreeTraits traits() const {return tr;} // Returns the traits class;

  Node_handle root() { return tree_root; }

  Kd_tree_rectangle<NT>* bounding_box() {return bbox; }

  int item_number() {return the_item_number;}

  // Print statistics of the tree.
  std::ostream& statistics (std::ostream& s) {
    s << "Tree statistics:" << std::endl;
    s << "Number of items stored: " 
		  << tree_root->num_items() << std::endl;
    s << " Tree depth: " << tree_root->depth() << std::endl;
    return s;
  }


};

} // namespace CGAL
#endif // CGAL_KD_TREE_H
