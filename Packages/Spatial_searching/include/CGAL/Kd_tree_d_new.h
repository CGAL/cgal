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
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Iso_rectangle_d.h>

namespace CGAL {

  // implementing backward compability to old kdtree

  template <class Traits>
  class Kdtree_d {

  typedef typename Traits::Iso_box_d Box;
  typedef typename Traits::Item Item;
 

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

  template <class OutputIterator, class FuzzyQueryItem>
	OutputIterator search(OutputIterator it, const FuzzyQueryItem& q) {
		it = t->search(it,q);
		return it;
	}

  template <class OutputIterator>
	OutputIterator report_all_points(OutputIterator it) 
	{it=t->report_all_points(it);
	 return it;}

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

  typedef Fuzzy_iso_box_d<Item,Iso_rectangle_2<R> > Iso_box_2;
  // work around, because old kd-tree constructor requires unneeded specification of dim
  class Iso_box_d {

  private:

  Iso_box_2 *b;

  public:

  //constuctor 
  Iso_box_d(const Item& p, const Item&q, int dim) {
	b=new Iso_box_2(p,q);
  }
 
  bool contains(const Item& p) const {	 
	return b->contains(p);
  }

  bool inner_range_intersects(const Kd_tree_rectangle<NT>* rectangle) const {   
	return b->inner_range_intersects(rectangle);
  }

  bool outer_range_is_contained_by(const Kd_tree_rectangle<NT>* rectangle) const { 
	return b->outer_range_is_contained_by(rectangle);
  }

  //destructor
  ~Iso_box_d() { delete b;}
};

  //constructor 
  Kdtree_interface_2d(unsigned int bucket_size=1, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		Kd_tree_traits_point<Item>(bucket_size,aspect_ratio,use_extended_nodes);
  }

    	
  // destructor
  ~Kdtree_interface_2d() {}

};

template <class Item,
	    class Splitter=Sliding_midpoint<Item>, 
	    class Separator=Plane_separator<typename Item::R::FT> > 
  class Kdtree_interface_3d : 
	public Kd_tree_traits_point<Item,Splitter,Separator> {

  public:

  typedef typename Item::R R;
  typedef typename Item::R::FT NT;

  typedef Fuzzy_iso_box_d<Item,Iso_cuboid_3<R> > Iso_box_3;
  // work around, because old kd-tree constructor requires unneeded specification of dim
  class Iso_box_d {

  private:

  Iso_box_3 *b;

  public:

  //constuctor 
  Iso_box_d(const Item& p, const Item&q, int dim) {
	b=new Iso_box_3(p,q);
  }
 
  bool contains(const Item& p) const {	 
	return b->contains(p);
  }

  bool inner_range_intersects(const Kd_tree_rectangle<NT>* rectangle) const {   
	return b->inner_range_intersects(rectangle);
  }

  bool outer_range_is_contained_by(const Kd_tree_rectangle<NT>* rectangle) const { 
	return b->outer_range_is_contained_by(rectangle);
  }

  //destructor
  ~Iso_box_d() { delete b;}
};

  //constructor 
  Kdtree_interface_3d(unsigned int bucket_size=1, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		Kd_tree_traits_point<Item>(bucket_size,aspect_ratio,use_extended_nodes);
  }

    	
  // destructor
  ~Kdtree_interface_3d() {}

};

template <class Item,
	  class Splitter=Sliding_midpoint<Item>, 
	  class Separator=Plane_separator<typename Item::R::FT> > 
  class Kdtree_interface : 
	public Kd_tree_traits_point<Item,Splitter,Separator> {

  public:

  typedef typename Item::R R;
  typedef typename Item::R::FT NT;

  typedef Fuzzy_iso_box_d<Item,Iso_rectangle_d<R> > Iso_box;
  // work around, because old kd-tree constructor requires unneeded specification of dim
  class Iso_box_d {

  private:

  Iso_box *b;

  public:

  //constuctor 
  Iso_box_d(const Item& p, const Item&q, int dim) {
	b=new Iso_box(p,q);
  }
 
  bool contains(const Item& p) const {	 
	return b->contains(p);
  }

  bool inner_range_intersects(const Kd_tree_rectangle<NT>* rectangle) const {   
	return b->inner_range_intersects(rectangle);
  }

  bool outer_range_is_contained_by(const Kd_tree_rectangle<NT>* rectangle) const { 
	return b->outer_range_is_contained_by(rectangle);
  }

  //destructor
  ~Iso_box_d() { delete b;}
};

  //constructor 
  Kdtree_interface(unsigned int bucket_size=1, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		Kd_tree_traits_point<Item>(bucket_size,aspect_ratio,use_extended_nodes);
  }

    	
  // destructor
  ~Kdtree_interface() {}

}; 

} // namespace CGAL
#endif // CGAL_KD_TREE_D_NEW_H
