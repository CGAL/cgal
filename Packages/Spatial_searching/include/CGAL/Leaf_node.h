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
namespace CGAL {

template <class Traits> // = Kd_tree_traits_Point>
class Leaf_node: public Base_node<Traits>  {
public:
  typedef typename Traits::Item Item;
  typedef typename Traits::Item_iterator Item_iterator;
private:
  unsigned int n;
  Item_iterator data;
public:
  // typedef typename Traits::NT NT;
  // typedef CGAL::Point_2< CGAL::Cartesian<NT> > Point_2D;
  
  const bool is_leaf() const { return 1;}
  const unsigned int size() const { return n;}
  
  Item_iterator const begin() const  {return data;}
  Item_iterator const end() const {return data + n;}
  Leaf_node(Points_container<Item>& c) :
    n(c.size()), data(new Item*[n]) {
    std::copy(c.begin(), c.end(), data);
  }


  ~Leaf_node() {
    delete []data;
  }
};


} // namespace CGAL
#endif // CGAL_LEAF_NODE_H
