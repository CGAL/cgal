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
// file          : include/CGAL/Orthogonal_standard_search.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef  ORTHOGONAL_STANDARD_SEARCH_H
#define  ORTHOGONAL_STANDARD_SEARCH_H
#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Weighted_Minkowski_distance.h>

namespace CGAL {

template <class Traits, class Query_item, class Distance>
class Orthogonal_standard_search {

public:

typedef typename Traits::Item Item;
typedef typename Traits::NT NT;
typedef std::pair<Item*,NT> Item_with_distance;

typedef Kd_tree_node<Traits> Node;
typedef Kd_tree<Traits> Tree;

//private:

typedef Item** Item_iterator;
typedef Kd_tree_rectangle<NT> Rectangle; 

private:

int number_of_internal_nodes_visited;
int number_of_leaf_nodes_visited;
int number_of_items_visited;

bool search_nearest;

NT multiplication_factor;
Query_item* query_object;
int total_item_number;
NT distance_to_root;   
int dim;

typedef std::list<Item_with_distance> NN_list;

NN_list l;
int max_k;
int actual_k;

Distance* distance_instance;

	inline bool branch(NT distance) {
		if (actual_k<max_k) return true;
		else 
			if (search_nearest) return ( distance < l.rbegin()->second * multiplication_factor ); 
			else return ( multiplication_factor * distance > l.begin()->second );
	};

	inline void insert(Item* I, NT dist) {
		bool insert;
		if (actual_k<max_k) insert=true;
		else 
			if (search_nearest) insert=( dist < l.rbegin()->second ); 
			else insert=(dist > l.begin()->second);
        if (insert) {
	   		actual_k++;	 	
			typename NN_list::iterator it=l.begin();
			for (; (it != l.end()); ++it) { if (dist < it->second) break;}
        		Item_with_distance NN_Candidate(I,dist);
        		l.insert(it,NN_Candidate);
        		if (actual_k > max_k) {
				actual_k--;
                if (search_nearest) l.pop_back();
				else l.pop_front();
        	};
  		}

	};

	
	public:

	template<class OutputIterator>  
	OutputIterator  the_k_neighbours(OutputIterator res)
	{   
		typename NN_list::iterator it=l.begin(); 
		for (; it != l.end(); it++) { *res= *it; res++; }
		return res;     
	};


    // constructor
    Orthogonal_standard_search(Tree& tree, Query_item& q,  
    Distance& d, int k, NT Eps, bool Search_nearest=true) {

	distance_instance=&d;

	multiplication_factor=
	distance_instance->transformed_distance(1.0+Eps);
        
	max_k=k;
	actual_k=0;
	search_nearest = Search_nearest; 
		
        if (search_nearest) 
		distance_to_root=
        	distance_instance->min_distance_to_queryitem(q,
						*(tree.bounding_box()));
        else 
		distance_to_root=
        	distance_instance->max_distance_to_queryitem(q,
						*(tree.bounding_box()));

        query_object = &q;
        dim=query_object->dimension();

        total_item_number=tree.item_number();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
       
        compute_neighbours_orthogonally(tree.root(), distance_to_root);
       
    }

    
    // Print statistics of the standard search process.
    void statistics () {
    	std::cout << "Standard search statistics:" << std::endl;
    	std::cout << "Number of internal nodes visited:" << number_of_internal_nodes_visited << std::endl;
    	std::cout << "Number of leaf nodes visited:" << number_of_leaf_nodes_visited << std::endl;
    	std::cout << "Number of items visited:" << number_of_items_visited << std::endl;
    }

    // destructor
    ~Orthogonal_standard_search() { 
		l.clear();  
    };

    private:
   
    void compute_neighbours_orthogonally(Node* N, NT rd) {
		
                if (!(N->is_leaf())) {
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->separator()->cutting_dimension();
                        NT old_off, new_rd;
                        NT new_off =
                        (*query_object)[new_cut_dim] - 
						N->separator()->cutting_value();
                        if ( ((new_off < 0.0) && (search_nearest)) ||
                        (( new_off >= 0.0) && (!search_nearest))  ) {
				compute_neighbours_orthogonally(N->lower(),rd);
                                if (search_nearest) {
                                	old_off= (*query_object)[new_cut_dim]-
								N->low_value();
                                	if (old_off>0.0) old_off=0.0;
                                }
				else 
				{	
                                	old_off= (*query_object)[new_cut_dim] - N->high_value();
					if (old_off<0.0) old_off=0.0;
                                }
                                new_rd=
                                distance_instance->
                                new_distance(rd,old_off,new_off,new_cut_dim);
				if (branch(new_rd)) compute_neighbours_orthogonally(N->upper(),new_rd);                               
                        }
                        else { // compute new distance
                                compute_neighbours_orthogonally(N->upper(),rd); 
				if (search_nearest) {
                                	old_off= N->high_value() - (*query_object)[new_cut_dim];
                                	// if (old_off>0.0) old_off=0.0;
				}
                                else 
                                {       
                                	old_off= N->low_value() - (*query_object)[new_cut_dim];
					// if (old_off<0.0) old_off=0.0;
				}  
                                new_rd=distance_instance->
                                new_distance(rd,old_off,new_off,new_cut_dim);
				if (branch(new_rd)) compute_neighbours_orthogonally(N->lower(),new_rd);       
                        }
                }
                else
				{
                  // n is a leaf
                  number_of_leaf_nodes_visited++;
                  if (N->size() > 0)
                  for (Item_iterator it=N->begin(); it != N->end(); it++) {
                        number_of_items_visited++;
			NT distance_to_query_object=
                        distance_instance->
                        distance(*query_object,**it);
                        insert(*it,distance_to_query_object);
                  }
		}
    }

    
    


    
   

}; // class 



} // namespace CGAL


#endif  // ORTHOGONAL_STANDARD_SEARCH
