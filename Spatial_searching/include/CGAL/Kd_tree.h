// Copyright (c) 2002,2011,2014 Utrecht University (The Netherlands), Max-Planck-Institute Saarbruecken (Germany).
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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>),
//               : Waqar Khan <wkhan@mpi-inf.mpg.de>

#ifndef CGAL_KD_TREE_H
#define CGAL_KD_TREE_H

#include <CGAL/license/Spatial_searching.h>


#include <CGAL/basic.h>
#include <CGAL/assertions.h>
#include <vector>

#include <CGAL/algorithm.h>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Splitters.h>
#include <CGAL/internal/Get_dimension_tag.h>

#include <deque>
#include <boost/container/deque.hpp>
#include <boost/optional.hpp>

#ifdef CGAL_HAS_THREADS
#include <CGAL/mutex.h>
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
  typedef Kd_tree_leaf_node<SearchTraits, Splitter, UseExtendedNode > Leaf_node;
  typedef Kd_tree_internal_node<SearchTraits, Splitter, UseExtendedNode > Internal_node;
  typedef Kd_tree<SearchTraits, Splitter> Tree;
  typedef Kd_tree<SearchTraits, Splitter,UseExtendedNode> Self;

  typedef Node* Node_handle;
  typedef const Node* Node_const_handle;
  typedef Leaf_node* Leaf_node_handle;
  typedef const Leaf_node* Leaf_node_const_handle;
  typedef Internal_node* Internal_node_handle;
  typedef const Internal_node* Internal_node_const_handle;
  typedef typename std::vector<const Point_d*>::const_iterator Point_d_iterator;
  typedef typename std::vector<const Point_d*>::const_iterator Point_d_const_iterator;
  typedef typename Splitter::Separator Separator;
  typedef typename std::vector<Point_d>::const_iterator iterator;
  typedef typename std::vector<Point_d>::const_iterator const_iterator;

  typedef typename std::vector<Point_d>::size_type size_type;

  typedef typename internal::Get_dimension_tag<SearchTraits>::Dimension D;

private:
  SearchTraits traits_;
  Splitter split;


  // wokaround for https://svn.boost.org/trac/boost/ticket/9332
#if   (_MSC_VER == 1800) && (BOOST_VERSION == 105500)
  std::deque<Internal_node> internal_nodes;
  std::deque<Leaf_node> leaf_nodes;
#else
  boost::container::deque<Internal_node> internal_nodes;
  boost::container::deque<Leaf_node> leaf_nodes;
#endif

  Node_handle tree_root;

  Kd_tree_rectangle<FT,D>* bbox;
  std::vector<Point_d> pts;

  // Instead of storing the points in arrays in the Kd_tree_node
  // we put all the data in a vector in the Kd_tree.
  // and we only store an iterator range in the Kd_tree_node.
  //
  std::vector<const Point_d*> data;


  #ifdef CGAL_HAS_THREADS
  mutable CGAL_MUTEX building_mutex;//mutex used to protect const calls inducing build()
  #endif
  bool built_;
  bool removed_;

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
    Leaf_node node(true , static_cast<unsigned int>(c.size()));
    std::ptrdiff_t tmp = c.begin() - data.begin();
    node.data = pts.begin() + tmp;

    leaf_nodes.push_back(node);
    Leaf_node_handle nh = &leaf_nodes.back();

   
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
    Internal_node node(false);
    internal_nodes.push_back(node);
    Internal_node_handle nh = &internal_nodes.back();

    Separator sep;
    Point_container c_low(c.dimension(),traits_);
    split(sep, c, c_low);
    nh->set_separator(sep);

    int cd  = nh->cutting_dimension();
    if(!c_low.empty()){
      nh->lower_low_val = c_low.tight_bounding_box().min_coord(cd);
      nh->lower_high_val = c_low.tight_bounding_box().max_coord(cd);
    }
    else{
      nh->lower_low_val = nh->cutting_value();
      nh->lower_high_val = nh->cutting_value();
    }
    if(!c.empty()){
      nh->upper_low_val = c.tight_bounding_box().min_coord(cd);
      nh->upper_high_val = c.tight_bounding_box().max_coord(cd);
    }
    else{
      nh->upper_low_val = nh->cutting_value();
      nh->upper_high_val = nh->cutting_value();
    }

    CGAL_assertion(nh->cutting_value() >= nh->lower_low_val);
    CGAL_assertion(nh->cutting_value() <= nh->upper_high_val);

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
    Internal_node node(false);
    internal_nodes.push_back(node);
    Internal_node_handle nh = &internal_nodes.back();
    Separator sep;

    Point_container c_low(c.dimension(),traits_);
    split(sep, c, c_low);
    nh->set_separator(sep);

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
    : traits_(traits),split(s), built_(false), removed_(false)
  {}

  template <class InputIterator>
  Kd_tree(InputIterator first, InputIterator beyond,
	  Splitter s = Splitter(),const SearchTraits traits=SearchTraits())
    : traits_(traits),split(s), built_(false), removed_(false)
  {
    pts.insert(pts.end(), first, beyond);
  }

  bool empty() const {
    return pts.empty();
  }

  void
  build()
  {
    // This function is not ready to be called when a tree already exists, one
    // must call invalidate_build() first.
    CGAL_assertion(!is_built());
    CGAL_assertion(!pts.empty());
    CGAL_assertion(!removed_);
    const Point_d& p = *pts.begin();
    typename SearchTraits::Construct_cartesian_const_iterator_d ccci=traits_.construct_cartesian_const_iterator_d_object();
    int dim = static_cast<int>(std::distance(ccci(p), ccci(p,0)));

    data.reserve(pts.size());
    for(unsigned int i = 0; i < pts.size(); i++){
      data.push_back(&pts[i]);
    }
    Point_container c(dim, data.begin(), data.end(),traits_);
    bbox = new Kd_tree_rectangle<FT,D>(c.bounding_box());
    if (c.size() <= split.bucket_size()){
      tree_root = create_leaf_node(c);
    }else {
      tree_root = create_internal_node(c, UseExtendedNode());
    }

    //Reorder vector for spatial locality
    std::vector<Point_d> ptstmp;
    ptstmp.resize(pts.size());
    for (std::size_t i = 0; i < pts.size(); ++i){
      ptstmp[i] = *data[i];
    }
    for(std::size_t i = 0; i < leaf_nodes.size(); ++i){
      std::ptrdiff_t tmp = leaf_nodes[i].begin() - pts.begin();
      leaf_nodes[i].data = ptstmp.begin() + tmp;
    }
    pts.swap(ptstmp);

    data.clear();

    built_ = true;
  }

private:
  //any call to this function is for the moment not threadsafe
  void const_build() const {
    #ifdef CGAL_HAS_THREADS
    //this ensure that build() will be called once
    CGAL_SCOPED_LOCK(building_mutex);
    if(!is_built())
    #endif
      const_cast<Self*>(this)->build(); //THIS IS NOT THREADSAFE
  }
public:

  bool is_built() const
  {
    return built_;
  }

  void invalidate_build()
  {
    if(removed_){
      // Walk the tree to collect the remaining points.
      // Writing directly to pts would likely work, but better be safe.
      std::vector<Point_d> ptstmp;
      //ptstmp.resize(root()->num_items());
      root()->tree_items(std::back_inserter(ptstmp));
      pts.swap(ptstmp);
      removed_=false;
      CGAL_assertion(is_built()); // the rest of the cleanup must happen
    }
    if(is_built()){
      internal_nodes.clear();
      leaf_nodes.clear();
      data.clear();
      delete bbox;
      built_ = false;
    }
  }

  void clear()
  {
    invalidate_build();
    pts.clear();
    removed_ = false;
  }

  void
  insert(const Point_d& p)
  {
    invalidate_build();
    pts.push_back(p);
  }

  template <class InputIterator>
  void
  insert(InputIterator first, InputIterator beyond)
  {
    invalidate_build();
    pts.insert(pts.end(),first, beyond);
  }

private:
  struct Equal_by_coordinates {
    SearchTraits const* traits;
    Point_d const* pp;
    bool operator()(Point_d const&q) const {
      typename SearchTraits::Construct_cartesian_const_iterator_d ccci=traits->construct_cartesian_const_iterator_d_object();
      return std::equal(ccci(*pp), ccci(*pp,0), ccci(q));
    }
  };
  Equal_by_coordinates equal_by_coordinates(Point_d const&p){
    Equal_by_coordinates ret = { &traits(), &p };
    return ret;
  }

public:
  void
  remove(const Point_d& p)
  {
    remove(p, equal_by_coordinates(p));
  }

  template<class Equal>
  void
  remove(const Point_d& p, Equal const& equal_to_p)
  {
#if 0
    // This code could have quadratic runtime.
    if (!is_built()) {
      std::vector<Point_d>::iterator pi = std::find_if(pts.begin(), pts.end(), equal_to_p);
      // Precondition: the point must be there.
      CGAL_assertion (pi != pts.end());
      pts.erase(pi);
      return;
    }
#endif
    bool success = remove_(p, 0, false, 0, false, root(), equal_to_p);
    CGAL_assertion(success);

    // Do not set the flag is the tree has been cleared.
    if(is_built())
      removed_ |= success;
  }
private:
  template<class Equal>
  bool remove_(const Point_d& p,
      Internal_node_handle grandparent, bool parent_islower,
      Internal_node_handle parent, bool islower,
      Node_handle node, Equal const& equal_to_p) {
    // Recurse to locate the point
    if (!node->is_leaf()) {
      Internal_node_handle newparent = static_cast<Internal_node_handle>(node);
      // FIXME: This should be if(x<y) remove low; else remove up;
      if (traits().construct_cartesian_const_iterator_d_object()(p)[newparent->cutting_dimension()] <= newparent->cutting_value()) {
	if (remove_(p, parent, islower, newparent, true, newparent->lower(), equal_to_p))
	  return true;
      }
      //if (traits().construct_cartesian_const_iterator_d_object()(p)[newparent->cutting_dimension()] >= newparent->cutting_value())
	return remove_(p, parent, islower, newparent, false, newparent->upper(), equal_to_p);

      CGAL_assertion(false); // Point was not found
    }

    // Actual removal
    Leaf_node_handle lnode = static_cast<Leaf_node_handle>(node);
    if (lnode->size() > 1) {
      iterator pi = std::find_if(lnode->begin(), lnode->end(), equal_to_p);
      // FIXME: we should ensure this never happens
      if (pi == lnode->end()) return false;
      iterator lasti = lnode->end() - 1;
      if (pi != lasti) {
	// Hack to get a non-const iterator
	std::iter_swap(pts.begin()+(pi-pts.begin()), pts.begin()+(lasti-pts.begin()));
      }
      lnode->drop_last_point();
    } else if (!equal_to_p(*lnode->begin())) {
      // FIXME: we should ensure this never happens
      return false;
    } else if (grandparent) {
      Node_handle brother = islower ? parent->upper() : parent->lower();
      if (parent_islower)
	grandparent->set_lower(brother);
      else
	grandparent->set_upper(brother);
    } else if (parent) {
      tree_root = islower ? parent->upper() : parent->lower();
    } else {
      clear();
    }
    return true;
  }

public:
  //For efficiency; reserve the size of the points vectors in advance (if the number of points is already known).
  void reserve(size_t size)
  {
    pts.reserve(size);
  }

  //Get the capacity of the underlying points vector.
  size_t capacity()
  {
    return pts.capacity();
  }


  template <class OutputIterator, class FuzzyQueryItem>
  OutputIterator
  search(OutputIterator it, const FuzzyQueryItem& q) const
  {
    if(! pts.empty()){

      if(! is_built()){
	const_build();
      }
      Kd_tree_rectangle<FT,D> b(*bbox);
      return tree_root->search(it,q,b);
    }
    return it;
  }


  template <class FuzzyQueryItem>
  boost::optional<Point_d>
  search_any_point(const FuzzyQueryItem& q) const
  {
    if(! pts.empty()){

      if(! is_built()){
	const_build();
      }
      Kd_tree_rectangle<FT,D> b(*bbox);
      return tree_root->search_any_point(q,b);
    }
    return boost::none;
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
    if(! pts.empty()){
      if(! is_built()){
	const_build();
      }
      root()->print();
    }else{
      std::cout << "empty tree\n";
    }
  }

  const Kd_tree_rectangle<FT,D>&
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
