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
// file          : include/CGAL/Orthogonal_standard_search.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
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
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Splitters.h>

namespace CGAL {

template <class GeomTraits, 
	  class Distance_=Euclidean_distance<GeomTraits>, 
          class Splitter_=Sliding_midpoint<GeomTraits> ,
	  class Tree_=Kd_tree<GeomTraits, Splitter_, Tag_true> >
class Orthogonal_standard_search {

public:

  typedef Splitter_ Splitter;
  typedef Tree_  Tree;
  typedef Distance_ Distance;
  typedef typename GeomTraits::Point Point;
  typedef typename GeomTraits::Point Query_item;
  typedef typename GeomTraits::NT NT;
  typedef std::pair<Point,NT> Point_with_distance;
  
  typedef typename Tree::Node_handle Node_handle;
  
  typedef typename Tree::Point_iterator Point_iterator;
  // af was: typedef typename Kd_tree<GeomTraits, Splitter>::Point_iterator Point_iterator;
  typedef Kd_tree_rectangle<GeomTraits> Rectangle; 

private:

int number_of_internal_nodes_visited;
int number_of_leaf_nodes_visited;
int number_of_items_visited;

bool search_nearest;

NT multiplication_factor;
Query_item query_object;
int total_item_number;
NT distance_to_root;   

typedef std::list<Point_with_distance> NN_list;

NN_list l;
int max_k;
int actual_k;


Distance* distance_instance;

	inline bool branch(NT distance) {
		if (actual_k<max_k) return true;
		else 
		  if (search_nearest) return 
		  ( distance * multiplication_factor < l.rbegin()->second);
		  else return 
		  ( distance > 
		    l.begin()->second * multiplication_factor);
	};

	inline void insert(Point* I, NT dist) {
		bool insert;
		if (actual_k<max_k) insert=true;
		else 
			if (search_nearest) insert=
			( dist < l.rbegin()->second ); 
			else insert=(dist > l.begin()->second);
        if (insert) {
	   		actual_k++;	 	
			typename NN_list::iterator it=l.begin();
			for (; (it != l.end()); ++it) 
			{ if (dist < it->second) break;}
        		Point_with_distance NN_Candidate(*I,dist);
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
	OutputIterator  the_k_neighbors(OutputIterator res)
	{   
		typename NN_list::iterator it=l.begin(); 
		for (; it != l.end(); it++) { *res= *it; res++; }
		return res;     
	}


    // constructor
    Orthogonal_standard_search(Tree& tree, const Query_item& q,  
    int k=1, NT Eps=NT(0.0), bool Search_nearest=true, const Distance& d=Distance()) {
   
	distance_instance=new Distance(d);

	multiplication_factor=
	d.transformed_distance(1.0+Eps);
        
	max_k=k;
	actual_k=0;
	search_nearest = Search_nearest; 
		
        // if (search_nearest) 
		distance_to_root=
        	d.min_distance_to_queryitem(q,
						*(tree.bounding_box()));
        // else 
	//	distance_to_root=
        //	distance_instance->max_distance_to_queryitem(q,
	//					*(tree.bounding_box()));

        query_object = q;

        total_item_number=tree.size();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
       
        compute_neighbors_orthogonally(tree.root(), distance_to_root);
       
    }

    
    // Print statistics of the standard search process.
    std::ostream& statistics (std::ostream& s) {
    	s << "Standard search statistics:" << std::endl;
    	s << "Number of internal nodes visited:" 
	<< number_of_internal_nodes_visited << std::endl;
    	s << "Number of leaf nodes visited:" 
	<< number_of_leaf_nodes_visited << std::endl;
    	s << "Number of items visited:" 
	<< number_of_items_visited << std::endl;
        return s;
    }

    // destructor
    ~Orthogonal_standard_search() { 
		l.clear();  
		delete distance_instance;
    };

    private:
   
    void compute_neighbors_orthogonally(Node_handle N, NT rd) {
		
      typename GeomTraits::Construct_cartesian_const_iterator construct_it;
      typename GeomTraits::Cartesian_const_iterator query_object_it = construct_it(query_object);
                if (!(N->is_leaf())) {
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->cutting_dimension();
                        NT old_off, new_rd;
                        NT new_off =
                        *(query_object_it + new_cut_dim) - 
					N->cutting_value();
                        if ( ((new_off < NT(0.0)) && (search_nearest)) ||
                        (( new_off >= NT(0.0)) && (!search_nearest))  ) {
				compute_neighbors_orthogonally(N->lower(),rd);
                                if (search_nearest) {
                                	old_off= *(query_object_it + new_cut_dim)-
							N->low_value();
                                	if (old_off>NT(0.0)) old_off=NT(0.0);
                                }
				else 
				{	
                                	old_off= *(query_object_it + new_cut_dim) 
					- N->high_value();
					if (old_off<NT(0.0)) old_off=NT(0.0);
                                }
                                new_rd=
                                distance_instance->
				new_distance(rd,old_off,new_off,new_cut_dim);
				if (branch(new_rd)) 
				compute_neighbors_orthogonally(N->upper(),
								new_rd);                               
                        }
                        else { // compute new distance
                                compute_neighbors_orthogonally(N->upper(),rd); 
				if (search_nearest) {
                                	old_off= N->high_value() - 
					*(query_object_it + new_cut_dim);
                                	if (old_off>NT(0.0)) old_off=NT(0.0);
				}
                                else 
                                {       
                                	old_off= N->low_value() - 
					*(query_object_it + new_cut_dim);
					if (old_off<NT(0.0)) old_off=NT(0.0);
				}  
                                new_rd=
                                distance_instance->
				new_distance(rd,old_off,new_off,new_cut_dim);
				if (branch(new_rd)) 
				compute_neighbors_orthogonally(N->lower(),
								new_rd);       
                        }
                }
                else
				{
                  // n is a leaf
                  number_of_leaf_nodes_visited++;
                  if (N->size() > 0)
                  for (Point_iterator it=N->begin(); it != N->end(); it++) {
                        number_of_items_visited++;
			NT distance_to_query_object=
                        distance_instance->
                        distance(query_object,**it);
                        insert(*it,distance_to_query_object);
                  }
		}
    }

    
    


    
   

}; // class 



} // namespace CGAL


#endif  // ORTHOGONAL_STANDARD_SEARCH
