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
// file          : include/CGAL/Leaf_node.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef CGAL_LEAF_NODE_H
#define CGAL_LEAF_NODE_H
#include <CGAL/Base_node.h>
// #include <CGAL/Cartesian.h>
// #include <iomanip>
// #include <CGAL/PS_stream.h>

namespace CGAL {

template <class Traits> // = Kd_tree_traits_Point>
class Leaf_node: public Base_node<Traits>  {
public:
  typedef typename Traits::Item Item;
  typedef typename Traits::Item_iterator Item_iterator;
private:
  int n;
  Item_iterator data;
public:
  // typedef typename Traits::NT NT;
  // typedef CGAL::Point_2< CGAL::Cartesian<NT> > Point_2D;
  
  const bool is_leaf() const { return 1;}
  const int size() const { return n;}
  
  Item_iterator const begin() const  {return data;}
  Item_iterator const end() const {return data + n;}
  Leaf_node(Points_container<Item>& c) :
    n(c.size()), data(new Item*[n]) {
    std::copy(c.begin(), c.end(), data);
    // std::cout << "Boxtree_leaf_node_d called" << std::endl;
  }

/*
  void data_to_postscript(PS_Stream& PS,
	const int i, const int j,
	const NT mini, const NT maxi,
	const NT minj, const NT maxj) {
	  // PS << border_color(RED);  works only for visual
	  for (Item_iterator it=begin(); it != end(); it++) { 
	  Point_2D p ( (*(*it))[i], (*(*it))[j] ); 
	  PS << p;
	}
  } */

  // removed default constructor !!
  ~Leaf_node() {
    // std::cout << "~Leaf_node called" << std::endl;
    delete []data;
  }
};


} // namespace CGAL
#endif // CGAL_LEAF_NODE_H
