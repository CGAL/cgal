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
// file          : include/CGAL/Binary_search_tree.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_BINARY_SEARCH_TREE_H
#define CGAL_BINARY_SEARCH_TREE_H
#include <CGAL/Binary_node.h>
namespace CGAL {


template <class Traits> 
class Binary_search_tree {
public:
  
  typedef typename Traits::Item Item;
  typedef typename std::list<Item>::iterator input_iterator;
  typedef typename Traits::NT NT;
  typedef Binary_node<Traits> Node;
  typedef Binary_search_tree<Traits> Tree;

private:
  Node* tree_root;
  Box<NT>* bbox;
  std::list<Item> pts;
  Traits tr;
  int the_item_number;

  // protected copy constructor
  Binary_search_tree(const Tree& tree) {};

public:

  // remove check_validity?, now commented out.
  Binary_search_tree(input_iterator first, input_iterator beyond,
	    Traits t = Traits(), bool check_validity=false) : tr(t) {
    assert(first != beyond);
    int dim = first->dimension();
    std::copy(first, beyond, std::back_inserter(pts));
    Points_container<Item> c(dim, pts.begin(), pts.end());
    // if (check_validity) { assert(c.is_valid()); }
		// std::cout << "validity of container used to store points:" 
		//	    <<  c.is_valid() << std::endl;}
    bbox = new Box<NT>(c.bounding_box());
    the_item_number=c.size();
    if (c.size() <= t.bucket_size())
      tree_root = new Node(c);
    else {
		if (t.use_extended_nodes())
		{tree_root = new Node(c,t,true); 
		 std::cout << "using extended internal nodes" << std::endl;}
		else
		{tree_root = new Node(c,t,false); 
		 std::cout << "not using extended internal nodes" 
			   << std::endl;}
	}
	// if (check_validity) { assert(is_valid()); 
		// std::cout << "validity of constructed binary tree:" 
		// <<  is_valid() << std::endl;
	//}
  }

    ~Binary_search_tree() {
                  delete tree_root; delete bbox;
	};


  Traits traits() const {return tr;} // Returns the traits class;

  Node* root() { return tree_root; }

  Box<NT>* bounding_box() {return bbox; }

  int item_number() {return the_item_number;}

 
  bool is_valid() { 
     // bool verbose = false, int level = 0) {
     // Perform internal consistency checks to verify the correctness
     // of the tree.
    std::list<Item> pts_1;
    std::back_insert_iterator<std::list<Item> > it(pts_1);
    //    unsigned int nit = items(it);
    root()->tree_items(it); // store the items in pts_1
    // check that pts_1 and pts contain the same stuff
    assert( pts_1.size() == root()->num_items());
    typename std::list<Item>::const_iterator i;
    for (i = pts.begin(); i != pts.end(); ++i) {
      typename std::list<Item>::iterator j = 
	std::find(pts_1.begin(), pts_1.end(), *i);
      assert(j != pts_1.end());
      assert(*j == *i);
    }
    return 1;
  }
  



  // Print statistics of the tree.
  void statistics (bool check_validity=false) {
    std::cout << "Tree statistics:" << std::endl;
    std::cout << "Number of items stored: " 
		  << tree_root->num_items() << std::endl;
    std::cout << " Tree depth: " << tree_root->depth() << std::endl;
    if (check_validity) { assert(is_valid());
        // std::cout << " Calling is_valid: " << is_valid() << std::endl;
    }

  }
};

} // namespace CGAL
#endif // CGAL_BINARY_SEARCH_TREE_H
