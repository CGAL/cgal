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
// release       : $CGAL_Revision: CGAL-2.5-I-99 $
// release_date  : $CGAL_Date: 2003/05/23 $
//
// file          : include/CGAL/Kd_tree.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_KD_TREE_H
#define CGAL_KD_TREE_H
#include <CGAL/basic.h>
#include <cassert>
#include<list>

#include <CGAL/algorithm.h>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Splitters.h>
#include <CGAL/Compact_container.h>


namespace CGAL {


  template <class GeomTraits, class Splitter_=Sliding_midpoint<GeomTraits>, class UseExtendedNode = Tag_true >
  class Kd_tree {

public:

    typedef Splitter_ Splitter;
  typedef typename GeomTraits::Point Point;
  typedef typename Splitter::Container Point_container;

  typedef typename GeomTraits::NT NT;
  typedef Kd_tree_node<GeomTraits, Splitter, UseExtendedNode > Node;
  typedef Kd_tree<GeomTraits, Splitter> Tree;

  typedef typename Compact_container<Node>::iterator Node_handle;
  typedef typename std::vector<Point*>::iterator Point_iterator;
  typedef typename Splitter::Separator Separator;

private:

    Splitter split;
  Compact_container<Node> nodes;

  Node_handle tree_root;

  Kd_tree_rectangle<GeomTraits>* bbox;
  std::list<Point> pts;

  // Instead of storing the points in arrays in the Kd_tree_node
  // we put all the data in a vector in the Kd_tree.
  // and we only store an iterator range in the Kd_tree_node.
  // 
  std::vector<Point*> data;
  Point_iterator data_iterator;
  GeomTraits tr;
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
    Node n;
    Node_handle nh = nodes.insert(n);
    
    nh->the_node_type = Node::EXTENDED_INTERNAL;

    Point_container
      c_low = Point_container(c.dimension());
    
    split(nh->separator(), c, c_low);
	        
    int cd  = nh->separator().cutting_dimension();
    
    nh->low_val = c_low.bounding_box().min_coord(cd);
   

    
    nh->high_val = c.bounding_box().max_coord(cd);

    
    assert(nh->separator().cutting_value() >= nh->low_val);
    assert(nh->separator().cutting_value() <= nh->high_val);

    

    if (c_low.size() > split.bucket_size())
      nh->lower_ch = create_internal_node_use_extension(c_low);
    else
      nh->lower_ch = create_leaf_node(c_low);

    if (c.size() > split.bucket_size())
      nh->upper_ch = create_internal_node_use_extension(c);
    else
      nh->upper_ch = create_leaf_node(c);

    
    return nh;
  }

  
  // Note also that I duplicated the code to get rid if the if's for
  // the boolean use_extension which was constant over the construction
  Node_handle 
  create_internal_node(Point_container& c) 
  {
    Node n;
    Node_handle nh = nodes.insert(n);
    
    nh->the_node_type = Node::INTERNAL;

    Point_container
    c_low = Point_container(c.dimension());
    
    split(nh->separator(), c, c_low);
	        
    if (c_low.size() > split.bucket_size())
      nh->lower_ch = create_internal_node(c_low);
    else
      nh->lower_ch = create_leaf_node(c_low);

    if (c.size() > split.bucket_size())
      nh->upper_ch = create_internal_node(c);
    else
      nh->upper_ch = create_leaf_node(c);

    return nh;
  }
  




public:

  //introduced for backward compability
  Kd_tree() {}
  
template <class InputIterator>
  Kd_tree(InputIterator first, InputIterator beyond,
	    Splitter s = Splitter()) : split(s) {
    assert(first != beyond);
    std::copy(first, beyond, std::back_inserter(pts));
    const Point& p = *pts.begin();
    typename GeomTraits::Construct_cartesian_const_iterator ccci;
    int dim = std::distance(ccci(p), ccci(p,0)); 

    data = std::vector<Point*>(pts.size()); // guarantees that iterators we store in Kd_tree_nodes stay valid
    data_iterator = data.begin();

    Point_container c(dim, pts.begin(), pts.end());

    bbox = new Kd_tree_rectangle<GeomTraits>(c.bounding_box());
    
    the_item_number=c.size();
    if (c.size() <= split.bucket_size())
      tree_root = create_leaf_node(c);
    else {
      tree_root = create_internal_node(c, UseExtendedNode()); 
    }
	
  }

 
  template <class OutputIterator, class FuzzyQueryItem>
	OutputIterator search(OutputIterator it, const FuzzyQueryItem& q) {
		Kd_tree_rectangle<GeomTraits> b(*bbox);
		tree_root->search(it,q,b);
		return it;
	}

  template <class OutputIterator>
	OutputIterator report_all_points(OutputIterator it) 
	{it=tree_root->tree_items(it);
	 return it;}

    ~Kd_tree() {  
		  delete bbox;
	};


  GeomTraits traits() const {return tr;} // Returns the traits class;

  Node_handle root()const { return tree_root; }

  Kd_tree_rectangle<GeomTraits>* bounding_box() const {return bbox; }

  int size() const {return the_item_number;}

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
