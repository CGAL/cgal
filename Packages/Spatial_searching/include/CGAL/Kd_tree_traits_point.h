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
// file          : include/CGAL/Kd_tree_traits_point.h
// package       : ASPAS
// revision      : 2.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef CGAL_KD_TREE_TRAITS_POINT_H
#define CGAL_KD_TREE_TRAITS_POINT_H
#include <CGAL/Splitters.h>

namespace CGAL {

  template <class Item_, 
            class Splitter=Sliding_midpoint<Item_>, class Separator_=Plane_separator<typename Item_::R::FT> >
  class Kd_tree_traits_point {

  public:
    typedef Item_ Item;
    typedef Separator_ Separator;
    typedef Item** Item_iterator;
    // CGAL dependency typedef typename Kernel_traits<Item>::Kernel K;
    // CGAL dependency typedef typename K::FT NT;
    typedef typename Item::R::FT NT;
    
    
  private:

    unsigned int the_bucket_size;
    NT the_aspect_ratio;
    bool use_extended_nodes_option;

  public:

       
        
	Kd_tree_traits_point(unsigned int bucket_size=1, 
			     NT aspect_ratio=NT(3), 
			     bool use_extended_nodes=true) {
		the_bucket_size = bucket_size;
		
		the_aspect_ratio = aspect_ratio;
		use_extended_nodes_option = use_extended_nodes;
	}

    	NT aspect_ratio() const {return the_aspect_ratio;}
	
        unsigned int bucket_size() const {return the_bucket_size;}
	bool use_extended_nodes() const {return use_extended_nodes_option;}

	// split c0 in c0 and c1
    	void split(Separator& sep,
    	           Point_container<Item>& c0, Point_container<Item>& c1) 
	{
		Splitter S;
		S(sep,c0,c1,the_aspect_ratio);
    	}

  };

} // namespace CGAL
#endif // KD_TREE_TRAITS_POINT_H
