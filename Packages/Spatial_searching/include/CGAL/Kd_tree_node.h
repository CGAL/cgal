// Copyright (c) 2002  Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $Source$
// $Revision$ $Date$
// $Name$
//
// Authors       : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_KD_TREE_NODE_H
#define CGAL_KD_TREE_NODE_H


#include <CGAL/Splitters.h>
#include <CGAL/Compact_container.h>
namespace CGAL {

  template <class SearchTraits, class Splitter, class UseExtendedNode> 
  class Kd_tree;

	template < class TreeTraits, class Splitter, class UseExtendedNode > 
	class Kd_tree_node {

	  friend class Kd_tree<TreeTraits,Splitter,UseExtendedNode>;

	  typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Node_handle Node_handle;
	enum Node_type {LEAF, INTERNAL, EXTENDED_INTERNAL};
	typedef typename TreeTraits::Point_d Point_d;

	typedef typename TreeTraits::FT FT;
	typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Separator Separator;
	typedef typename Kd_tree<TreeTraits,Splitter,UseExtendedNode>::Point_d_iterator Point_d_iterator;

        private:

	// node type identifier
	Node_type the_node_type;

     	// private variables for leaf nodes
	unsigned int n; // denotes number of items in a leaf node
	Point_d_iterator data; // iterator to data in leaf node

    	// private variables for internal nodes

	  Node_handle lower_ch, upper_ch;

  	Separator sep;

	// private variables for extended internal nodes
	FT low_val;
  	FT high_val;
                
	public:
		
	  void *   for_compact_container() const { return lower_ch.for_compact_container(); }
	  void * & for_compact_container()       { return lower_ch.for_compact_container(); }

	// default constructor
	Kd_tree_node() {};

	Kd_tree_node(Node_type t ) : the_node_type(t) {};

	Kd_tree_node(int n_, Node_type t ) : n(n_), the_node_type(t) {};

        // members for all nodes
	inline bool is_leaf() const { return (the_node_type==LEAF);}

	// members for leaf nodes only
  	inline unsigned int size() const { return n;}
  
  	inline Point_d_iterator begin() const  {return data;}
  	inline Point_d_iterator end() const {return data + n;}

	// members for internal node and extended internal node

	inline Node_handle lower() const { return lower_ch; }
  	inline Node_handle upper() const { return upper_ch; }
  	
  	// inline Separator& separator() {return sep; }
  	// use instead
  	
  	inline FT cutting_value() const 
  	{return sep.cutting_value();}
  	
  	inline int cutting_dimension() const 
  	{return sep.cutting_dimension();}

	// members for extended internal node only
	inline FT low_value() const { return low_val; }
  	inline FT high_value() const { return high_val; }
       

	  Separator& separator() {
	    return sep;
	  }
	

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
			  for (Point_d_iterator i=begin(); i != end(); i++) 
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
			      Kd_tree_rectangle<TreeTraits>& b) {
		if (is_leaf()) { 
			if (n>0) 
			for (Point_d_iterator i=begin(); i != end(); i++) 
				if (q.contains(**i)) 
                                {*it=**i; ++it;}
                }
		else {
                        // after splitting b denotes the lower part of b
		  Kd_tree_rectangle<TreeTraits> b_upper(b);
		  b.split(b_upper, sep.cutting_dimension(),
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
		};
	        return it;				
	}

        
   };


} // namespace CGAL
#endif // CGAL_KDTREE_NODE_H
