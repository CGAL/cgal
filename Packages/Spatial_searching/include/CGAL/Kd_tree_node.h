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
// revision      : 1.4 
// revision_date : 2003/02/01 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef CGAL_KD_TREE_NODE_H
#define CGAL_KD_TREE_NODE_H

#include <CGAL/Kd_tree_traits_point.h>

namespace CGAL {

	template < class Traits > 
	class Kd_tree_node {

	public:

	enum Node_type {LEAF, INTERNAL, EXTENDED_INTERNAL};
	typedef typename Traits::Item Item;
  	typedef typename Traits::Item_iterator Item_iterator;
	typedef typename Traits::NT NT;
	typedef typename Traits::Separator Separator;

        private:

	// node type identifier
	Node_type the_node_type;

     	// private variables for leaf nodes
	unsigned int n; // denotes number of items in a leaf node
  	Item_iterator data; // iterator to data in leaf node

    	// private variables for internal nodes

    	Kd_tree_node* lower_ch;
  	Kd_tree_node* upper_ch;
  	Separator* sep;

	// private variables for extended internal nodes
	NT low_val;
  	NT high_val;
                
	public:
		
	// default constructor
	Kd_tree_node() {};

	// constructor for leaf node
        /*
        Kd_tree_node(Point_container<Item>& c) :
    		n(c.size()), data(new Item*[c.size()]) {
		the_node_type=LEAF;
    		std::copy(c.begin(), c.end(), data);
  	}; */

        Kd_tree_node(Point_container<Item>& c) :
    		n(c.size()) {
                if (n>0) data=new Item*[n];
		the_node_type=LEAF;
    		std::copy(c.begin(), c.end(), data);
  	};
	// constructor for internal node or extended internal node;
	Kd_tree_node(Point_container<Item>& c, Traits& t, bool use_extension) {
		
		
		if (use_extension) 
			the_node_type=EXTENDED_INTERNAL;
		else 
			the_node_type=INTERNAL;

		

    		Point_container<Item> c_low = Point_container<Item>(c.dimension());
			Kd_tree_rectangle<NT> bbox(c.bounding_box());

    		sep = t.split(c, c_low);
	
    		int cd  = sep->cutting_dimension();

		if (use_extension) {
    			low_val = bbox.lower(cd);
    			high_val = bbox.upper(cd);
		};

    		if (c_low.size() > t.bucket_size())
      			lower_ch = new Kd_tree_node<Traits>(c_low,t,use_extension);
    		else
      			lower_ch = new Kd_tree_node<Traits>(c_low);

    		if (c.size() > t.bucket_size())
      			upper_ch = new Kd_tree_node<Traits>(c,t,use_extension);
    		else
      			upper_ch = new Kd_tree_node<Traits>(c);

  	};

        // members for all nodes
	inline bool is_leaf() const { return (the_node_type==LEAF);}

	// members for leaf nodes only
  	inline unsigned int size() const { return n;}
  
  	inline Item_iterator begin() const  {return data;}
  	inline Item_iterator end() const {return data + n;}

	// members for internal node and extended internal node
	inline Kd_tree_node* lower() const { return lower_ch; }
  	inline Kd_tree_node* upper() const { return upper_ch; }
  	inline Separator* separator() const {return sep; }

	// members for extended internal node only
	inline NT low_value() const { return low_val; }
  	inline NT high_value() const { return high_val; }
       
        ~Kd_tree_node() {switch (the_node_type) {

			case LEAF: {
				if (n>0) delete []data;}
			     	break;

			case INTERNAL: { 
				delete sep; delete lower_ch; delete upper_ch;}
				break;
			case EXTENDED_INTERNAL:
				delete sep; delete lower_ch; delete upper_ch;
				break;	
			default:{
				std::cerr << "Node corrupted\n";
				}
		}
	};

	unsigned int num_items() {
			if (is_leaf()) return size();
			else return lower()->num_items() + upper()->num_items();
		}

	int depth(const int current_max_depth) {
			if (is_leaf()) return current_max_depth;
			else return std::max( lower()->depth(current_max_depth + 1),
		       upper()->depth(current_max_depth + 1));
		}

	int depth() { return depth(1); }

	template <class OutputIterator>
	void tree_items(OutputIterator& it) {
            	if (is_leaf()) 
				for (Item_iterator i=begin(); i != end(); i++) *it=**i;
			else {
				lower()->tree_items(it);
				upper()->tree_items(it);
			}
		}

   };


} // namespace CGAL
#endif // CGAL_KDTREE_NODE_H
