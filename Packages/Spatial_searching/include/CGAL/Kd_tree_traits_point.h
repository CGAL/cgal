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

namespace CGAL {

  template <class Separator_, class Item_>
  class Kd_tree_traits_point {

  public:
    typedef Separator_ Separator;
    typedef Item_ Item;
    typedef Kernel_traits<Item>::Kernel K;
    typedef K::FT NT;
    typedef Item** Item_iterator;
    typedef typename std::list<Item>::iterator InputIterator;
    typedef std::pair<Item*,NT> Item_with_distance;

  private:

    int Bucket_size;
    Split_rule Selected_split_rule;
    NT Aspect_ratio;
    bool Use_extended_nodes;

  public:

        

	Kd_tree_traits_point(int bucket_size=1, Split_rule My_split_rule=SLIDING_MIDPOINT,
						 NT aspect_ratio=3.0, bool use_extended_nodes=true) {
		Bucket_size = bucket_size;
		Selected_split_rule = My_split_rule;
		Aspect_ratio = aspect_ratio;
		Use_extended_nodes = use_extended_nodes;
	}

    	NT aspect_ratio() {return Aspect_ratio;}
	Split_rule  selected_split_rule() {return Selected_split_rule;}
	Shrink_rule selected_shrink_rule() {return NONE;}
    	unsigned int bucket_size() {return Bucket_size;}
	bool use_extended_nodes() {return Use_extended_nodes;}
        bool nearest_neighbour_search() {return Use_nearest_neighbour_search;}

	// split c0 in c0 and c1
    	Separator* split(Points_container<Item>& c0, Points_container<Item>& c1) {
		Separator* sep;

		switch (Selected_split_rule) {

			case SLIDING_MIDPOINT:
				{Sliding_MidPoint<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep,true);}
				break;

			case SLIDING_FAIR:
				{Sliding_Fair<Item> M;
				sep=M.rule(c0,aspect_ratio());
				c0.split_container(c1,sep,true);}
				break;

			case FAIR:
				{Fair<Item> M;
				sep=M.rule(c0,aspect_ratio());
				c0.split_container(c1,sep);}
				break;

			case MEDIAN_OF_MAX_SPREAD:
				{Median_Of_Max_Spread<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case MEDIAN_OF_BOX:
				{Median_Of_Box<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case MIDPOINT_OF_MAX_SPREAD:
				{MidPoint_Of_Max_Spread<Item> M;
				sep=M.rule(c0);
				c0.split_container(c1,sep);}
			    break;

			case MIDPOINT_OF_BOX:
				{MidPoint_Of_Box<Item> M;
				sep=M.rule(c0);
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
