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
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef CGAL_KD_TREE_TRAITS_POINT_H
#define CGAL_KD_TREE_TRAITS_POINT_H
#include <CGAL/Splitting_rules.h>
namespace CGAL {

  template <class Separator_, class Item_>
  class Kd_tree_traits_point {

  public:
    typedef Separator_ Separator;
    typedef Item_ Item;
    typedef Item** Item_iterator;
    // CGAL dependency typedef typename Kernel_traits<Item>::Kernel K;
    // CGaL dependency typedef typename K::FT NT;
    typedef typename Item::R::FT NT;
    typedef std::pair<Item*,NT> Item_with_distance;
    typedef typename Split_rule_enumeration::Split_rule split_rule;

  private:

    unsigned int the_bucket_size;
    split_rule the_selected_split_rule;
    NT the_aspect_ratio;
    bool use_extended_nodes_option;

  public:

        
	Kd_tree_traits_point(unsigned int bucket_size=1, 
			     split_rule My_split_rule=Split_rule_enumeration::SLIDING_MIDPOINT,
			     NT aspect_ratio=3.0, 
			     bool use_extended_nodes=true) {
		the_bucket_size = bucket_size;
		the_selected_split_rule = My_split_rule;
		the_aspect_ratio = aspect_ratio;
		use_extended_nodes_option = use_extended_nodes;
	}

    	NT aspect_ratio() const {return the_aspect_ratio;}
	split_rule  selected_split_rule() const {return the_selected_split_rule;}

        unsigned int bucket_size() const {return the_bucket_size;}
	bool use_extended_nodes() const {return use_extended_nodes_option;}

	// split c0 in c0 and c1
    	Separator* split(Point_container<Item>& c0, Point_container<Item>& c1) 
	{
		Separator* sep;

		switch (the_selected_split_rule) {

			case Split_rule_enumeration::SLIDING_MIDPOINT:
				{Sliding_midpoint<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep,true);}
				break;

			case Split_rule_enumeration::SLIDING_FAIR:
				{Sliding_fair<Item> M;
				sep=M.rule(c0,aspect_ratio());
				c0.split_container(c1,sep,true);}
				break;

			case Split_rule_enumeration::FAIR:
				{Fair<Item> M;
				sep=M.rule(c0,aspect_ratio());
				c0.split_container(c1,sep);}
				break;

			case Split_rule_enumeration::MEDIAN_OF_MAX_SPREAD:
				{Median_of_max_spread<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case Split_rule_enumeration::MEDIAN_OF_BOX:
				{Median_of_box<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case Split_rule_enumeration::MIDPOINT_OF_MAX_SPREAD:
				{Midpoint_of_max_spread<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case Split_rule_enumeration::MIDPOINT_OF_BOX:
				{Midpoint_of_box<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep);}
			    break;

			default:
				std::cerr << "Split rule corrupted\n";
		}

		return sep;
    }

  };

} // namespace CGAL
#endif // KD_TREE_TRAITS_POINT_H
