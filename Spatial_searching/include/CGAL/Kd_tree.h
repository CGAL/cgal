// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL$
// $Id$
// 
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_KD_TREE_H
#define CGAL_KD_TREE_H

#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <vector>

#include <CGAL/algorithm.h>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Splitters.h>
#include <CGAL/Compact_container.h>

#ifdef CGAL_HAS_THREADS
#include <boost/thread/mutex.hpp>
#endif

namespace CGAL {

  //template <class SearchTraits, class Splitter_=Median_of_rectangle<SearchTraits>, class UseExtendedNode = Tag_true >
template <class SearchTraits, class Splitter_=Sliding_midpoint<SearchTraits>, class UseExtendedNode = Tag_true >
class Kd_tree {

public:
  typedef SearchTraits Traits;
  typedef Splitter_ Splitter;
  typedef typename SearchTraits::Point_d Point_d;
  typedef typename Splitter::Container Point_container;
  
  typedef typename SearchTraits::FT FT;
  typedef Kd_tree_node<SearchTraits, Splitter, UseExtendedNode > Node;
  typedef Kd_tree<SearchTraits, Splitter> Tree;
  typedef Kd_tree<SearchTraits, Splitter,UseExtendedNode> Self;

  typedef typename Compact_container<Node>::iterator Node_handle;
  typedef typename Compact_container<Node>::const_iterator Node_const_handle;
  typedef typename std::vector<const Point_d*>::const_iterator Point_d_iterator;
  typedef typename std::vector<const Point_d*>::const_iterator Point_d_const_iterator;
  typedef typename Splitter::Separator Separator;
  typedef typename std::vector<Point_d>::const_iterator iterator;
  typedef typename std::vector<Point_d>::const_iterator const_iterator;
  
  typedef typename std::vector<Point_d>::size_type size_type;

private:
  SearchTraits traits_;
  Splitter split;
  Compact_container<Node> nodes;

  Node_handle tree_root;

  Kd_tree_rectangle<FT>* bbox;
  std::vector<Point_d> pts;

  // Instead of storing the points in arrays in the Kd_tree_node
  // we put all the data in a vector in the Kd_tree.
  // and we only store an iterator range in the Kd_tree_node.
  // 
  std::vector<const Point_d*> data;


  #ifdef CGAL_HAS_THREADS
  mutable boost::mutex building_mutex;//mutex used to protect const calls inducing build()
  #endif
  bool built_;

  // protected copy constructor
  Kd_tree(const Tree& tree)
    : traits_(tree.traits_),built_(tree.built_)
  {};


  // Instead of the recursive construction of the tree in the class Kd_tree_node
  // we do this in the tree class. The advantage is that we then can optimize
  // the allocation of the nodes.

  // The leaf node
  Node_handle 
  create_leaf_node(Point_container& c)
  {
    Node_handle nh = nodes.emplace(static_cast<unsigned int>(c.size()), Node::LEAF);

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
    Node_handle nh = nodes.emplace(Node::EXTENDED_INTERNAL);
    
    Point_container c_low(c.dimension(),traits_);
    split(nh->separator(), c, c_low);
	        
    int cd  = nh->separator().cutting_dimension();
    
    nh->low_val = c_low.bounding_box().min_coord(cd);
    nh->high_val = c.bounding_box().max_coord(cd);
    
    CGAL_assertion(nh->separator().cutting_value() >= nh->low_val);
    CGAL_assertion(nh->separator().cutting_value() <= nh->high_val);

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
    Node_handle nh = nodes.emplace(Node::INTERNAL);
    
    Point_container c_low(c.dimension(),traits_);
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

  Kd_tree(Splitter s = Splitter(),const SearchTraits traits=SearchTraits())
    : traits_(traits),split(s), built_(false)
  {}
  
  template <class InputIterator>
  Kd_tree(InputIterator first, InputIterator beyond,
	  Splitter s = Splitter(),const SearchTraits traits=SearchTraits()) 
    : traits_(traits),split(s), built_(false) 
  {
    pts.insert(pts.end(), first, beyond);
  }

  bool empty() const {
    return pts.empty();
  }
  
  void 
  build()
  {
    const Point_d& p = *pts.begin();
    typename SearchTraits::Construct_cartesian_const_iterator_d ccci=traits_.construct_cartesian_const_iterator_d_object();
    int dim = static_cast<int>(std::distance(ccci(p), ccci(p,0))); 

    data.reserve(pts.size());
    for(unsigned int i = 0; i < pts.size(); i++){
      data.push_back(&pts[i]);
    }
    Point_container c(dim, data.begin(), data.end(),traits_);
    bbox = new Kd_tree_rectangle<FT>(c.bounding_box());
    if (c.size() <= split.bucket_size()){
      tree_root = create_leaf_node(c);
    }else {
      tree_root = create_internal_node(c, UseExtendedNode()); 
    }
    built_ = true;
  }

private:  
  //any call to this function is for the moment not threadsafe
  void const_build() const {
    #ifdef CGAL_HAS_THREADS
    //this ensure that build() will be called once
    boost::mutex::scoped_lock scoped_lock(building_mutex);
    if(!is_built())
    #endif
      const_cast<Self*>(this)->build(); //THIS IS NOT THREADSAFE
  }
public:
  
  bool is_built() const
  {
    return built_;
  }

  void invalidate_built()
  {
    if(is_built()){
      nodes.clear();
      data.clear();
      delete bbox;
      built_ = false;
    }
  }
  
  void clear()
  {
    invalidate_built();
    pts.clear();
  }
  
  void
  insert(const Point_d& p)
  {
    invalidate_built();
    pts.push_back(p);
  } 
 
  template <class InputIterator>
  void 
  insert(InputIterator first, InputIterator beyond)
  {
    invalidate_built();
    pts.insert(pts.end(),first, beyond);
  }


  template <class OutputIterator, class FuzzyQueryItem>
  OutputIterator 
  search(OutputIterator it, const FuzzyQueryItem& q) const
  {
    if(! pts.empty()){
  
      if(! is_built()){
	const_build();
      }
      Kd_tree_rectangle<FT> b(*bbox);
      tree_root->search(it,q,b);
    }
    return it;
  }

  ~Kd_tree() {
    if(is_built()){
      delete bbox;
    }
  }


  const SearchTraits&
  traits() const 
  {
    return traits_;
  }

  Node_const_handle 
  root() const 
  { 
    if(! is_built()){
      const_build();
    }
    return tree_root; 
  }

  Node_handle 
  root()
  {
    if(! is_built()){
      build();
    }
    return tree_root; 
  }

  void
  print() const
  {
    if(! is_built()){
      const_build();
    }
    root()->print();
  }

  const Kd_tree_rectangle<FT>&
  bounding_box() const 
  {
    if(! is_built()){
      const_build();
    }
    return *bbox; 
  }

  const_iterator
  begin() const
  {
    return pts.begin();
  }

  const_iterator
  end() const
  {
    return pts.end();
  }

  size_type 
  size() const 
  {
    return pts.size();
  }

  // Print statistics of the tree.
  std::ostream& 
  statistics(std::ostream& s) const
  {
    if(! is_built()){
      const_build();
    }
    s << "Tree statistics:" << std::endl;
    s << "Number of items stored: " 
      << root()->num_items() << std::endl;
    s << "Number of nodes: " 
      << root()->num_nodes() << std::endl;
    s << " Tree depth: " << root()->depth() << std::endl;
    return s;
  }


};

} // namespace CGAL

#endif // CGAL_KD_TREE_H
