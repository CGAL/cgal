// Copyright (c) 2002 Utrecht University (The Netherlands).
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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_KD_TREE_H
#define CGAL_KD_TREE_H
#include <CGAL/basic.h>
#include <cassert>
#include<vector>

#include <CGAL/algorithm.h>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Splitters.h>
#include <CGAL/Compact_container.h>


namespace CGAL {


template <class SearchTraits, class Splitter_=Sliding_midpoint<SearchTraits>, class UseExtendedNode = Tag_true >
class Kd_tree {

public:

  typedef Splitter_ Splitter;
  typedef typename SearchTraits::Point_d Point_d;
  typedef typename Splitter::Container Point_container;
  
  typedef typename SearchTraits::FT FT;
  typedef Kd_tree_node<SearchTraits, Splitter, UseExtendedNode > Node;
  typedef Kd_tree<SearchTraits, Splitter> Tree;

  typedef typename Compact_container<Node>::iterator Node_handle;
  typedef typename std::vector<Point_d*>::iterator Point_d_iterator;
  typedef typename Splitter::Separator Separator;
  std::vector<Point_d>::iterator iterator;

private:

  Splitter split;
  Compact_container<Node> nodes;

  Node_handle tree_root;

  Kd_tree_rectangle<SearchTraits>* bbox;
  std::vector<Point_d> pts;

  // Instead of storing the points in arrays in the Kd_tree_node
  // we put all the data in a vector in the Kd_tree.
  // and we only store an iterator range in the Kd_tree_node.
  // 
  std::vector<Point_d*> data;
  SearchTraits tr;

  // protected copy constructor
  Kd_tree(const Tree& tree) 
  {};


  // Instead of the recursive construction of the tree in the class Kd_tree_node
  // we do this in the tree class. The advantage is that we then can optimize
  // the allocation of the nodes.

  // The leaf node
  Node_handle 
  create_leaf_node(Point_container& c)
  {
    Node_handle nh = nodes.construct_insert(c.size(), Node::LEAF);

    nh->data = c.begin();
    return nh;
  }

 
  // The internal node

  Node_handle 
  create_internal_node(Point_container& c, const Tag_true&)
  {
    return create_internal_node_use_extension(c);
  }

  Node_handle 
  create_internal_node(Point_container& c, const Tag_false&)
  {
    return create_internal_node(c);
  }

 
 
  // TODO: Similiar to the leaf_init function above, a part of the code should be
  //       moved to a the class Kd_tree_node.
  //       It is not proper yet, but the goal was to see if there is
  //       a potential performance gain through the Compact_container
  Node_handle 
  create_internal_node_use_extension(Point_container& c) 
  {
    Node_handle nh = nodes.construct_insert(Node::EXTENDED_INTERNAL);
    
    Point_container c_low(c.dimension());
    
    split(nh->separator(), c, c_low);
	        
    int cd  = nh->separator().cutting_dimension();
    
    nh->low_val = c_low.bounding_box().min_coord(cd);
    nh->high_val = c.bounding_box().max_coord(cd);
    
    assert(nh->separator().cutting_value() >= nh->low_val);
    assert(nh->separator().cutting_value() <= nh->high_val);

    if (c_low.size() > split.bucket_size()){
      nh->lower_ch = create_internal_node_use_extension(c_low);
    }else{
      nh->lower_ch = create_leaf_node(c_low);
    }
    if (c.size() > split.bucket_size()){
      nh->upper_ch = create_internal_node_use_extension(c);
    }else{
      nh->upper_ch = create_leaf_node(c);
    }
    
    return nh;
  }

  
  // Note also that I duplicated the code to get rid if the if's for
  // the boolean use_extension which was constant over the construction
  Node_handle 
  create_internal_node(Point_container& c) 
  {
    Node_handle nh = nodes.construct_insert(Node::INTERNAL);
    
    Point_container c_low(c.dimension());
    
    split(nh->separator(), c, c_low);
	        
    if (c_low.size() > split.bucket_size()){
      nh->lower_ch = create_internal_node(c_low);
    }else{
      nh->lower_ch = create_leaf_node(c_low);
    }
    if (c.size() > split.bucket_size()){
      nh->upper_ch = create_internal_node(c);
    }else{
      nh->upper_ch = create_leaf_node(c);
    }
    return nh;
  }
  


public:

  //introduced for backward compability
  Kd_tree() {}
  
  template <class InputIterator>
  Kd_tree(InputIterator first, InputIterator beyond,
	Splitter s = Splitter()) : split(s) 
  {
    assert(first != beyond);
    std::copy(first, beyond, std::back_inserter(pts));
    const Point_d& p = *pts.begin();
    typename SearchTraits::Construct_cartesian_const_iterator_d ccci;
    int dim = std::distance(ccci(p), ccci(p,0)); 

    data.reserve(pts.size());
    for(int i = 0; i < pts.size(); i++){
      data[i] = &pts[i];
    }
    Point_container c(dim, data.begin(), data.end());

    bbox = new Kd_tree_rectangle<SearchTraits>(c.bounding_box());
    
    if (c.size() <= split.bucket_size()){
      tree_root = create_leaf_node(c);
    }else {
      tree_root = create_internal_node(c, UseExtendedNode()); 
    }
  }

 
  template <class OutputIterator, class FuzzyQueryItem>
  OutputIterator 
  search(OutputIterator it, const FuzzyQueryItem& q) 
  {
    Kd_tree_rectangle<SearchTraits> b(*bbox);
    tree_root->search(it,q,b);
    return it;
  }

  ~Kd_tree() {  
    delete bbox;
  }


  SearchTraits 
  traits() const 
  {
    return tr;
  }

  Node_handle 
  root() const 
  { 
    return tree_root; 
  }

  const Kd_tree_rectangle<SearchTraits>&
  bounding_box() const 
  {
    return *bbox; 
  }

  iterator
  begin() const
  {
    return pts.begin();
  }

  iterator
  end() const
  {
    return pts.end();
  }

  int 
  size() const 
  {
    return pts.size();
  }

  // Print statistics of the tree.
  std::ostream& 
  statistics(std::ostream& s) 
  {
    s << "Tree statistics:" << std::endl;
    s << "Number of items stored: " 
      << tree_root->num_items() << std::endl;
    s << " Tree depth: " << tree_root->depth() << std::endl;
    return s;
  }


};

} // namespace CGAL
#endif // CGAL_KD_TREE_H
