// Copyright (c) 2002  Utrecht University (The Netherlands).
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
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_KD_TREE_D_NEW_H
#define CGAL_KD_TREE_D_NEW_H
#include <CGAL/Kd_tree.h>
#include <CGAL/Fuzzy_iso_rectangle_d.h>
#include <CGAL/Fuzzy_iso_box_d.h>
#include <CGAL/Iso_rectangle_2.h>
#include <CGAL/Iso_cuboid_3.h>
#include <CGAL/Iso_rectangle_d.h>


namespace CGAL {

  // implementing backward compability to old kdtree

  template <class TreeTraits>
  class Kdtree_d {

  public:
  typedef typename TreeTraits::Iso_box_d Box;

  private:

  typedef typename TreeTraits::Point Point;
  typedef std::list<Point> Point_list; 

  Kd_tree<TreeTraits>* t;

  public:
  
  // constructor

  Kdtree_d(int k = 2) {}

  void build(Point_list &l) {
	t = new Kd_tree<TreeTraits>(l.begin(), l.end());}

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

  template <class Point,
	    class Splitter=Sliding_midpoint<Point> > 
  class Kdtree_interface_2d : 
	public Search_traits<Point,Splitter> {

  public:

  typedef typename Kernel_traits<Point>::Kernel K;
  typedef typename K::FT NT;

  typedef Fuzzy_iso_box_d<Point,Iso_rectangle_2<K> > Iso_box_2;
  // work around, because old kd-tree constructor requires unneeded specification of dim
  class Iso_box_d {

  private:

  Iso_box_2 *b;

  public:

  //constuctor 
  Iso_box_d(const Point& p, const Point&q, int dim) {
	b=new Iso_box_2(p,q);
  }
 
  bool contains(const Point& p) const {	 
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
  Kdtree_interface_2d(unsigned int bucket_size=100, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		Search_traits<Point>(bucket_size,aspect_ratio,use_extended_nodes);
  }

    	
  // destructor
  ~Kdtree_interface_2d() {}

};



template <class Point,
	    class Splitter=Sliding_midpoint<Point> > 
  class Kdtree_interface_3d : 
	public Search_traits<Point,Splitter> {

  public:

  typedef typename Kernel_traits<Point>::Kernel K;
  typedef typename K::FT NT;

  typedef Fuzzy_iso_box_d<Point,Iso_cuboid_3<K> > Iso_box_3;
  // work around, because old kd-tree constructor requires unneeded specification of dim
  class Iso_box_d {

  private:

  Iso_box_3 *b;

  public:

  //constuctor 
  Iso_box_d(const Point& p, const Point&q, int dim) {
	b=new Iso_box_3(p,q);
  }
 
  bool contains(const Point& p) const {	 
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
  Kdtree_interface_3d(unsigned int bucket_size=100, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		Search_traits<Point>(bucket_size,aspect_ratio,use_extended_nodes);
  }

    	
  // destructor
  ~Kdtree_interface_3d() {}

};

template <class Point,
	  class Splitter=Sliding_midpoint<Point> > 
  class Kdtree_interface : 
	public Search_traits<Point,Splitter> {

  public:

  typedef typename Kernel_traits<Point>::Kernel K;
  typedef typename K::FT NT;

  typedef Fuzzy_iso_rectangle_d<Point,Iso_rectangle_d<K> > Iso_box;
  // work around, because old kd-tree constructor requires unneeded specification of dim
  class Iso_box_d {

  private:

  Iso_box *b;

  public:

  //constuctor 
  Iso_box_d(const Point& p, const Point&q, int dim) {
	b=new Iso_box(p,q);
  }
 
  bool contains(const Point& p) const {	 
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
  Kdtree_interface(unsigned int bucket_size=100, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		Search_traits<Point>(bucket_size,aspect_ratio,use_extended_nodes);
  }

    	
  // destructor
  ~Kdtree_interface() {}

}; 

} // namespace CGAL
#endif // CGAL_KD_TREE_D_NEW_H
