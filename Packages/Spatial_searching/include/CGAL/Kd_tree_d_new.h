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
// file          : include/CGAL/Kd_tree_d_new.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 2.4 
// revision_date : 2003/02/01 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_KD_TREE_D_NEW_H
#define CGAL_KD_TREE_D_NEW_H
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_iso_box_d.h>
#include <CGAL/Iso_rectangle_2.h>

namespace CGAL {

  // implementing backward compability to old kdtree

  template <class Traits>
  class Kdtree_d : public Kd_tree<Traits> {

  typedef typename Traits::Iso_box_d Box;
  typedef typename Traits::Item Item;

  private: 

  Kd_tree<Traits>* t;

  public:
  
  // constructor

  Kdtree_d(int k = 2) {}

  void build(list<Item> &l) {
	t = new Kd_tree<Traits>(l.begin(), l.end());}

  void delete_all() {delete t;}

  // not implemented for Kd_tree
  bool is_valid() {return true;}

  // not implemented for Kd_tree
  void dump() {}

  ~Kdtree_d() {}
        

  };

  template <class Item,
	    class Splitter=Sliding_midpoint<Item>, 
	    class Separator=Plane_separator<typename Item::R::FT> > 
  class Kdtree_interface_2d : 
	public Kd_tree_traits_point<Item,Splitter,Separator> {

  public:

  typedef typename Item::R R;
  typedef typename Item::R::FT NT;
  typedef Fuzzy_iso_box_d<Item,Iso_rectangle_2<R> > Iso_box_d;

  // constructor 
  Kdtree_interface_2d(unsigned int bucket_size=1, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		Kd_tree_traits_point<Item>(bucket_size,aspect_ratio,use_extended_nodes);
  }

    	
  // destructor
  ~Kdtree_interface_2d() {}

};


} // namespace CGAL
#endif // CGAL_KD_TREE_D_NEW_H
