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
// file          : include/CGAL/Kd_tree_node.h
// package       : ASPAS
// revision      : 2.4 
// revision_date : 2003/02/01 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_KD_TREE_NODE_H
#define CGAL_KD_TREE_NODE_H

#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Compact_container.h>
namespace CGAL {

  template < class Traits > 
  class Kd_tree;

	template < class Traits > 
	class Kd_tree_node {

	  friend class Kd_tree<Traits>;

	  typedef Kd_tree_node<Traits> Node;

	  typedef typename Compact_container<Node>::iterator Node_handle;

	enum Node_type {LEAF, INTERNAL, EXTENDED_INTERNAL};
	typedef typename Traits::Item Item;
	typedef   std::vector<Item*>::iterator Item_iterator;
	typedef typename Traits::NT NT;
	typedef typename Traits::Separator Separator;

        private:

	// node type identifier
	Node_type the_node_type;

     	// private variables for leaf nodes
	unsigned int n; // denotes number of items in a leaf node
	Item_iterator data; // iterator to data in leaf node

    	// private variables for internal nodes

	  Node_handle lower_ch, upper_ch;

  	Separator sep;

	// private variables for extended internal nodes
	NT low_val;
  	NT high_val;
                
	public:
		
	  // TODO: I have to ask Sylvain how to get rid if this:
	  void* p;
	  void *   for_compact_container() const { return p; }
	  void * & for_compact_container()       { return p; }

	// default constructor
	Kd_tree_node(): the_node_type() {};

        // members for all nodes
	inline bool is_leaf() const { return (the_node_type==LEAF);}

	// members for leaf nodes only
  	inline unsigned int size() const { return n;}
  
  	inline Item_iterator begin() const  {return data;}
  	inline Item_iterator end() const {return data + n;}

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
			  for (Item_iterator i=begin(); i != end(); i++) 
				{*it=**i; ++it;} 
			}
		else {
			it=lower_ch->tree_items(it);  
			it=upper_ch->tree_items(it); 
		};
		return it;
	}

        template <class OutputIterator, class Rectangle>
	OutputIterator tree_items_in_rectangle(OutputIterator it, 
	Rectangle& r, Kd_tree_rectangle<NT>* b, NT eps) {
            	if (is_leaf()) { 
			if (n>0) 
			for (Item_iterator i=begin(); i != end(); i++) 
				if (r.has_on_bounded_side(**i)) 
                                {*it=**i; ++it;}
                }
		else {
                        // after splitting b denotes the lower part of b
			Kd_tree_rectangle<NT>* 
			b_upper=b->split(sep.cutting_dimension(),
					      sep.cutting_value());
                             
			if (b->is_enclosed_by_dilated_rectangle(r,eps)) 	
			   it=lower_ch->tree_items(it);
			else
		           if (b->intersects_eroded_rectangle(r,eps)) 
			   it=lower_ch->tree_items_in_rectangle(it,r,b,eps);

                        if (b_upper->is_enclosed_by_dilated_rectangle(r,eps))        
			    it=upper_ch->tree_items(it);
			else
			    if (b_upper->intersects_eroded_rectangle(r,eps)) 
			    it=upper_ch->tree_items_in_rectangle(it,r,b_upper,eps);

		        delete b_upper;
		};
	        return it;				
	}

   template <class OutputIterator>
	OutputIterator tree_items_in_sphere(OutputIterator it, Item& center, 
		NT min_squared_radius, NT squared_radius, 
		NT max_squared_radius, Kd_tree_rectangle<NT>* b) {
            	if (is_leaf()) { 
		    if (n>0) 
		    for (Item_iterator item_it=begin(); 
				       item_it != end(); item_it++) 
                       {// test whether the squared distance 
			// between **item_it and center 
			// is at most the squared_radius
			NT distance=NT(0);
                        int dim=center.dimension();
				  
		        for (int i = 0; 
			(i < dim) && (distance <= squared_radius); ++i) {
				distance += 
				(center[i]-(**item_it)[i]) * 
				(center[i]-(**item_it)[i]);
			}
				  
				  if (distance <= squared_radius) 
                                   {*it=**item_it; ++it;}
			     }   
                }
		else {
                             // after splitting b denotes the lower part of b
			     Kd_tree_rectangle<NT>* 
			     b_upper=b->split(sep.cutting_dimension(),
					      sep.cutting_value());
                             
                             
			     if // maximal range query encloses b
	(b->max_squared_Euclidean_distance_to_point_is_at_most
	    (center,max_squared_radius)) 	
				it=lower_ch->tree_items(it);
			     else 
		           	if  // minimal range query intersects b
	(b->min_squared_Euclidean_distance_to_point_is_at_most
	(center,min_squared_radius)) 
				it=lower_ch->tree_items_in_sphere
	(it,center,min_squared_radius,squared_radius,max_squared_radius,b);

			     // the same for b_upper
                              if // maximal range query encloses upper_b
	(b_upper->max_squared_Euclidean_distance_to_point_is_at_most
	(center,max_squared_radius)) 	
				it=upper_ch->tree_items(it);
			     else 
		           	if  // minimal range query intersects upper_b
	(b_upper->min_squared_Euclidean_distance_to_point_is_at_most
	(center,min_squared_radius)) 
				it=upper_ch->tree_items_in_sphere
	(it,center,min_squared_radius,squared_radius,
		max_squared_radius,b_upper);
				
			     delete b_upper;
		}
	return it;
	}
	

   };


} // namespace CGAL
#endif // CGAL_KDTREE_NODE_H
