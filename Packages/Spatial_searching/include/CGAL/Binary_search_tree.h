// ======================================================================
//
// Copyright (c) 2001 The CGAL Consortium
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
// package       : APSPAS
// revision      : 1.0 
// revision_date : 2001/06/12 
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
//
// ======================================================================

#ifndef CGAL_BINARY_SEARCH_TREE_H
#define CGAL_BINARY_SEARCH_TREE_H

#include <CGAL/Leaf_node.h>
#include <CGAL/Extended_internal_node.h>
#include <CGAL/Cartesian.h>
#include <CGAL/Segment_2.h>
#include <iomanip>
#include <CGAL/IO/PS_Stream.h>

namespace CGAL {

/*  static cast does not work without specififying specialization parameters
template <class Tree>
unsigned int num_items(const Tree* root) {
 if (root->is_leaf()) {
        Leaf_node *L=static_cast<Tree*>(root);
        return L->size();
 }
 else return num_items(root->lower()) + num_items(root->upper());
} */

template <class Tree>
unsigned int num_items(const Tree* root) {
 if (root->is_leaf()) return root->size();
 else return num_items(root->lower()) + num_items(root->upper());
}

template <class Tree>
int depth(const Tree* root, const int current_max_depth) {
  if (root->is_leaf()) return current_max_depth;
  else return std::max(depth(root->lower(), current_max_depth + 1),
		       depth(root->upper(), current_max_depth + 1));
}

template <class Tree>
int depth(const Tree* root) { return depth(root, 1); }

template <class T>
struct get_val : public std::unary_function<const T,
  typename std::iterator_traits<T>::reference> {
  typename std::iterator_traits<T>::reference operator()(const T x)
    {return *x;}
};

template <class Tree, class It>
void tree_items(Tree* root, It& it) {
  typedef Tree::Item Item;
  if (root->is_leaf()) std::transform(root->begin(), root->end(), it,
  			      get_val<Item*>());
  else {
    tree_items(root->lower(), it);
    tree_items(root->upper(), it);
  }
}

// template < class R, class Traits = Search_tree_traits_2d<R> > simplified to
template <class Traits> // = Kd_tree_traits_2d> //class The_node_type = Internal_node >
class Binary_search_tree {
public:
  // Iterator over items.
  // typedef typename Traits::Item_iterator Item_iterator;
  typedef Traits::Item Item;
  typedef Traits::InputIterator InputIterator;
  typedef Item::FT NT;
  typedef Base_node<Traits> Node;
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

  // CREATION
  Binary_search_tree(InputIterator first, InputIterator beyond, int bucket_size=1, bool check_validity=false,
	    Traits t = Traits()) : tr(t) {
    // std::cout << "Binary search tree started" << std::endl;
    assert(first != beyond);
    int dim = first->dimension();
    // std::cout << "dim =" << dim << std::endl;
    std::copy(first, beyond, std::back_inserter(pts));
    Points_container<Item> c(dim, pts.begin(), pts.end());
    // std::cout << "\n Built initial container." << std::endl;
    if (check_validity) {std::cout <<  c.is_valid() << std::endl;}
    // std::cout << c;
    // std::cout << "\n\n Binary_search_tree called ...\n" << std::endl;
    bbox = new Box<NT>(c.bounding_box());
    // std::cout << "constructed bbox using c.bounding_box()" << std::endl;
    the_item_number=c.size();
    if (c.size() <= bucket_size)
      tree_root = new Leaf_node<Traits>(c);
    else {
		if (t.use_extended_nodes())
			tree_root = new Extended_internal_node<Traits>(c,bucket_size);
		else
			tree_root = new Internal_node<Traits>(c,bucket_size);
	}

    // std::cout << "\n Finished tree building." << std::endl;
  }

    ~Binary_search_tree() {// std::cout << "~Boxtree_d called" << std::endl;
                  delete tree_root; delete bbox;
                  // std::cout << "~Boxtree_d ready" << std::endl;
	};
        // Destructor. All nodes in the tree are deleted.

// Member Functions
  template <class OutputIterator>
  unsigned int items(OutputIterator& it) {
    tree_items(&root(), it);
    return num_items(&root());
  }

  void generate_postscript_file(const char* filename, const float width,
	  const int i, const int j) {
	  
	  
	  
	  float box_height = bbox->upper(j) - bbox->lower(j); 
	  float box_width  = bbox->upper(i) - bbox->lower(i);
	  const float height= width * (box_height/box_width);
	  PS_Stream::PS_BBox bb(bbox->lower(i)-0.01*box_width, bbox->lower(j)-0.01*box_height, 
							bbox->upper(i)+0.01*box_width, bbox->upper(j)+0.01*box_height);
	  PS_Stream PS(bb,height,filename,PS_Stream::QUIET_EPS);   
	  // removed CGAL::
	  /* PS.init(bbox->lower(i),bbox->upper(i),bbox->lower(j));
	  cgalize(PS);        // removed CGAL::
	  PS.display(); */
	  // test it
	  PS << point_style(PS_Stream::FDOT);
	  PS << point_size(1);
      PS << line_width(1);
	  
	  tree_root->data_to_postscript(PS, i, j, bbox->lower(i), bbox->upper(i),
		  bbox->lower(j), bbox->upper(j));
  }

  Traits traits() const {return tr;} // Returns the traits class;

  // attention changed on Jan 3, 2001
  Node* root() { return tree_root; }

  // Return
  Box<NT>* bounding_box() {return bbox; }

  int item_number() {return the_item_number;}

  //  Item_iterator items_begin() {};
        // Iterator over all items in the tree.

  //  Item_iterator items_end();
        // Corresponding past-the-end iterator.

  // Jan 3, 2001 removed & before root twice
  bool is_valid(bool verbose = false, int level = 0) {
        // Perform internal consistency checks to verify the correctness
        // of the tree.
    std::list<Item> pts_1;
    std::back_insert_iterator<std::list<Item> > it(pts_1);
    //    unsigned int nit = items(it);
    tree_items(root(), it); // store the items in pts_1
    // check that pts_1 and pts contain the same stuff
    assert( pts_1.size() == num_items(root()));
    std::list<Item>::const_iterator i;
    for (i = pts.begin(); i != pts.end(); ++i) {
      std::list<Item>::iterator j = std::find(pts_1.begin(), pts_1.end(), *i);
      assert(j != pts_1.end());
      assert(*j == *i);
    }
    return 1;
  }

  // Print statistics of the tree.
  void statistics (bool check_validity=false) {
    std::cout << "Tree statistics:" << std::endl;
    std::cout << "Number of items stored: " << num_items(tree_root) << std::endl;
    std::cout << " Tree depth: " << depth(tree_root) << std::endl;
    if (check_validity) {
        std::cout << " Calling is_valid: " << is_valid() << std::endl;
    }

  }
};

} // namespace CGAL
#endif // CGAL_BINARY_SEARCH_TREE_H
