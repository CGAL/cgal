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
#include <CGAL/PS_Stream.h>

namespace CGAL {


template <class Traits> 
class Binary_search_tree {
public:
  
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

  Binary_search_tree(InputIterator first, InputIterator beyond,
	    Traits t = Traits(), bool check_validity=false) : tr(t) {
    assert(first != beyond);
    int dim = first->dimension();
    std::copy(first, beyond, std::back_inserter(pts));
    Points_container<Item> c(dim, pts.begin(), pts.end());
    if (check_validity) {std::cout << "validity of container used to store points:" <<  c.is_valid() << std::endl;}
    bbox = new Box<NT>(c.bounding_box());
    the_item_number=c.size();
    if (c.size() <= t.bucket_size())
      tree_root = new Leaf_node<Traits>(c);
    else {
		if (t.use_extended_nodes())
		{tree_root = new Extended_internal_node<Traits>(c,t); std::cout << "using extended internal nodes" << std::endl;}
		else
		{tree_root = new Internal_node<Traits>(c,t); std::cout << "not using extended internal nodes" << std::endl;}
	}
	if (check_validity) {std::cout << "validity of constructed binary tree:" <<  is_valid() << std::endl;}
  }

    ~Binary_search_tree() {
                  delete tree_root; delete bbox;
	};

  void generate_postscript_file(const char* filename, const float width,
	  const int i, const int j) {

      typedef CGAL::Point_2< CGAL::Cartesian<NT> > Point_2D;
      typedef CGAL::Segment_2< CGAL::Cartesian<NT> > Segment_2D;
	  
	  float box_height = bbox->upper(j) - bbox->lower(j); 
	  float box_width  = bbox->upper(i) - bbox->lower(i);
	  const float height= width * (box_height/box_width);
	  PS_Stream::PS_BBox bb(bbox->lower(i)-0.2*box_width, bbox->lower(j)-0.2*box_height, 
	  					bbox->upper(i)+0.2*box_width, bbox->upper(j)+0.2*box_height);
	  // PS_Stream::PS_BBox bb(bbox->lower(i), bbox->lower(j), 
	  //						bbox->upper(i), bbox->upper(j));
	  
	  PS_Stream PS(bb,height,filename,PS_Stream::QUIET_EPS);   
	  PS << point_style(PS_Stream::FDOT);
	  PS << point_size(1);
          PS << line_width(0.5);
          // PS << border(1);
          
	  Point_2D p00(bbox->lower(i),bbox->lower(j));
          Point_2D p01(bbox->lower(i),bbox->upper(j));
          Point_2D p11(bbox->upper(i),bbox->upper(j));
          Point_2D p10(bbox->upper(i),bbox->lower(j));

	  Segment_2D s0(p00,p01);
          Segment_2D s1(p01,p11);
          Segment_2D s2(p11,p10);
          Segment_2D s3(p10,p00);

	  // PS << border_color(BLACK);   works only for Visual
          PS << s0 << s1 << s2 << s3;
	  tree_root->data_to_postscript(PS, i, j, bbox->lower(i), bbox->upper(i),
		  bbox->lower(j), bbox->upper(j));
  }

  /* previous version
  void generate_postscript_file(const char* filename, const float width,
	  const int i, const int j) {
	  
	  float box_height = bbox->upper(j) - bbox->lower(j); 
	  float box_width  = bbox->upper(i) - bbox->lower(i);
	  const float height= width * (box_height/box_width);
	  PS_Stream::PS_BBox bb(bbox->lower(i)-0.01*box_width, bbox->lower(j)-0.01*box_height, 
							bbox->upper(i)+0.01*box_width, bbox->upper(j)+0.01*box_height);
	  PS_Stream PS(bb,height,filename,PS_Stream::QUIET_EPS);   
	  // removed CGAL::
	  // PS.init(bbox->lower(i),bbox->upper(i),bbox->lower(j));
	  // cgalize(PS);        // removed CGAL::
	  // PS.display(); 
	  // test it
	  PS << point_style(PS_Stream::FDOT);
	  PS << point_size(1);
      PS << line_width(1);
	  
	  tree_root->data_to_postscript(PS, i, j, bbox->lower(i), bbox->upper(i),
		  bbox->lower(j), bbox->upper(j));
  } */

  Traits traits() const {return tr;} // Returns the traits class;

  Node* root() { return tree_root; }

  Box<NT>* bounding_box() {return bbox; }

  int item_number() {return the_item_number;}

  bool is_valid(bool verbose = false, int level = 0) {
     // Perform internal consistency checks to verify the correctness
     // of the tree.
    std::list<Item> pts_1;
    std::back_insert_iterator<std::list<Item> > it(pts_1);
    //    unsigned int nit = items(it);
    root()->tree_items(it); // store the items in pts_1
    // check that pts_1 and pts contain the same stuff
    assert( pts_1.size() == root()->num_items());
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
    std::cout << "Number of items stored: " << tree_root->num_items() << std::endl;
    std::cout << " Tree depth: " << tree_root->depth() << std::endl;
    if (check_validity) {
        std::cout << " Calling is_valid: " << is_valid() << std::endl;
    }

  }
};

} // namespace CGAL
#endif // CGAL_BINARY_SEARCH_TREE_H
