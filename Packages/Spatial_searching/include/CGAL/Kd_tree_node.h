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
// file          : include/CGAL/Kd_tree_node.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_KD_TREE_NODE_H
#define CGAL_KD_TREE_NODE_H

#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Compact_container.h>
namespace CGAL {

  template < class TreeTraits > 
  class Kd_tree;

	template < class TreeTraits > 
	class Kd_tree_node {

	  friend class Kd_tree<TreeTraits>;

	  typedef typename Kd_tree<TreeTraits>::Node_handle Node_handle;
	enum Node_type {LEAF, INTERNAL, EXTENDED_INTERNAL};
	typedef typename TreeTraits::Point Point;

	typedef typename TreeTraits::NT NT;
	typedef typename TreeTraits::Separator Separator;
	typedef   typename Kd_tree<TreeTraits>::Point_iterator Point_iterator;

        private:

	// node type identifier
	Node_type the_node_type;

     	// private variables for leaf nodes
	unsigned int n; // denotes number of items in a leaf node
	Point_iterator data; // iterator to data in leaf node

    	// private variables for internal nodes

	  Node_handle lower_ch, upper_ch;

  	Separator sep;

	// private variables for extended internal nodes
	NT low_val;
  	NT high_val;
                
	public:
		
	  void *   for_compact_container() const { return lower_ch.for_compact_container(); }
	  void * & for_compact_container()       { return lower_ch.for_compact_container(); }

	// default constructor
	Kd_tree_node() {};

        // members for all nodes
	inline bool is_leaf() const { return (the_node_type==LEAF);}

	// members for leaf nodes only
  	inline unsigned int size() const { return n;}
  
  	inline Point_iterator begin() const  {return data;}
  	inline Point_iterator end() const {return data + n;}

	// members for internal node and extended internal node

	inline Node_handle lower() const { return lower_ch; }
  	inline Node_handle upper() const { return upper_ch; }
  	
  	// inline Separator& separator() {return sep; }
  	// use instead
  	
  	inline NT cutting_value() const 
  	{return sep.cutting_value();}
  	
  	inline int cutting_dimension() const 
  	{return sep.cutting_dimension();}

	// members for extended internal node only
	inline NT low_value() const { return low_val; }
  	inline NT high_value() const { return high_val; }
       
        
	

	unsigned int num_items() {
			if (is_leaf()) return size();
			else 
			return lower()->num_items() + upper()->num_items();
		}

	int depth(const int current_max_depth) {
			if (is_leaf()) return current_max_depth;
			else return 
			std::max( lower()->depth(current_max_depth + 1),
		       	upper()->depth(current_max_depth + 1));
		}

	int depth() { return depth(1); }

	template <class OutputIterator>
	OutputIterator tree_items(OutputIterator it) {
            	if (is_leaf()) 
                        { 
		          if (n>0) 
			  for (Point_iterator i=begin(); i != end(); i++) 
				{*it=**i; ++it;} 
			}
		else {
			it=lower_ch->tree_items(it);  
			it=upper_ch->tree_items(it); 
		};
		return it;
	}

        template <class OutputIterator, class FuzzyQueryItem>
	OutputIterator search(OutputIterator it, const FuzzyQueryItem& q,
			      Kd_tree_rectangle<NT>* b) {
		if (is_leaf()) { 
			if (n>0) 
			for (Point_iterator i=begin(); i != end(); i++) 
				if (q.contains(**i)) 
                                {*it=**i; ++it;}
                }
		else {
                        // after splitting b denotes the lower part of b
			Kd_tree_rectangle<NT>* 
			b_upper=b->split(sep.cutting_dimension(),
					      sep.cutting_value());
                             
			if (q.outer_range_is_contained_by(b)) 	
			   it=lower_ch->tree_items(it);
			else
		           if (q.inner_range_intersects(b)) 
			   it=lower_ch->search(it,q,b);

                        if  (q.outer_range_is_contained_by(b_upper))     
			    it=upper_ch->tree_items(it);
			else
			    if (q.inner_range_intersects(b_upper)) 
			    it=upper_ch->search(it,q,b_upper);
		        delete b_upper;
		};
	        return it;				
	}

        
   };


} // namespace CGAL
#endif // CGAL_KDTREE_NODE_H
