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
// file          : include/CGAL/Kd_tree_traits_point.h
// package       : APSPAS
// revision      : 1.0 
// revision_date : 2001/06/12 
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
//
// ======================================================================


#ifndef CGAL_KD_TREE_TRAITS_POINT_H
#define CGAL_KD_TREE_TRAITS_POINT_H
#include <CGAL/Splitting_rules.h>
#include <CGAL/Kernel_d/Point_d.h>

namespace CGAL {

  template <class Separator, class Item>
  class Kd_tree_traits_point {
  public:

	// int bucket_size;
    // Split_rule Selected_split_rule;
	// Shrink_rule Selected_shrink_rule;

    // typedef Point_d<R> Item;
    // typedef Plane_separator<NT> Separator;

    typedef Separator Separator;
    typedef Item Item;
    typedef Item::FT NT;
    typedef std::list<Item>::iterator InputIterator;
    typedef std::pair<Item*,NT> Item_with_distance;

    static NT Aspect_ratio() {return 3.0;}

    // static unsigned int bucket_size() {return 1;}
	// static unsigned int dimension() {return 10;}

	static Split_rule  selected_split_rule() {return SLIDING_MIDPOINT;}
	static Shrink_rule selected_shrink_rule() {return NONE;}

	static bool use_extended_nodes() {return true;}

	// split c0 in c0 and c1
    Separator* split(Points_container<Item>& c0, Points_container<Item>& c1,
                     Split_rule The_split_rule) {
		Separator* sep;

		switch (The_split_rule) {

			case SLIDING_MIDPOINT:
				{Sliding_MidPoint<Item> M;
				sep=M.Rule(c0);
				c0.split_container(c1,sep,true);}
				break;

			case SLIDING_FAIR:
				{Sliding_Fair<Item> M;
				sep=M.Rule(c0,Aspect_ratio());
				c0.split_container(c1,sep,true);}
				break;

			case FAIR:
				{Fair<Item> M;
				sep=M.Rule(c0,Aspect_ratio());
				c0.split_container(c1,sep);}
				break;

			case MEDIAN_OF_MAX_SPREAD:
				{Median_Of_Max_Spread<Item> M;
				sep=M.Rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case MEDIAN_OF_BOX:
				{Median_Of_Box<Item> M;
				sep=M.Rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case MIDPOINT_OF_MAX_SPREAD:
				{MidPoint_Of_Max_Spread<Item> M;
				sep=M.Rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case MIDPOINT_OF_BOX:
				{MidPoint_Of_Box<Item> M;
				sep=M.Rule(c0);
				c0.split_container(c1,sep);}
			    break;

			default: 
				std::cerr << "Split rule corrupted\n";
				break;
		}
                
		return sep;
    }
    
  };
		
} // namespace CGAL
#endif // KD_TREE_TRAITS_POINT_H
