// Copyright (c) 2002 Utrecht University (The Netherlands).
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
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef  ORTHOGONAL_K_NEIGHBOR_SEARCH_H
#define  ORTHOGONAL_K_NEIGHBOR_SEARCH_H
#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Splitters.h>

namespace CGAL {

template <class SearchTraits, 
  class Distance_=Euclidean_distance<SearchTraits>,
          class Splitter_=Sliding_midpoint<SearchTraits> ,
	  class Tree_=Kd_tree<SearchTraits, Splitter_, Tag_true> >
class Orthogonal_k_neighbor_search {

public:

  typedef Splitter_ Splitter;
  typedef Tree_  Tree;
  typedef Distance_ Distance;
  typedef typename SearchTraits::Point_d Point_d;
  typedef typename Distance::Query_item Query_item;

  typedef typename SearchTraits::FT FT;
  typedef std::pair<Point_d,FT> Point_with_transformed_distance;
  
  typedef typename Tree::Node_handle Node_handle;
  
  typedef typename Tree::Point_d_iterator Point_d_iterator;

private:

int number_of_internal_nodes_visited;
int number_of_leaf_nodes_visited;
int number_of_items_visited;

bool search_nearest;

FT multiplication_factor;
Query_item query_object;
int total_item_number;
FT distance_to_root;   

typedef std::list<Point_with_transformed_distance> NN_list;

public:
  typedef typename NN_list::const_iterator iterator;
private:

NN_list l;
int max_k;
int actual_k;


Distance distance_instance;

	inline bool branch(FT distance) {
		if (actual_k<max_k) return true;
		else 
		  if (search_nearest) return 
		  ( distance * multiplication_factor < l.rbegin()->second);
		  else return 
		  ( distance > 
		    l.begin()->second * multiplication_factor);
	};

	inline void insert(Point_d* I, FT dist) {
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
        		Point_with_transformed_distance NN_Candidate(*I,dist);
        		l.insert(it,NN_Candidate);
        		if (actual_k > max_k) {
				actual_k--;
                if (search_nearest) l.pop_back();
				else l.pop_front();
        	};
  		}

	};

	
	public:
  /*
	template<class OutputIterator>  
	OutputIterator  the_k_neighbors(OutputIterator res)
	{   
		typename NN_list::iterator it=l.begin(); 
		for (; it != l.end(); it++) { *res= *it; res++; }
		return res;     
	}
  */
  iterator begin() const
  {
    return l.begin();
  }

  iterator end() const
  {
    return l.end();
  }


    // constructor
    Orthogonal_k_neighbor_search(Tree& tree, const Query_item& q,  
    int k=1, FT Eps=FT(0.0), bool Search_nearest=true, const Distance& d=Distance())
      : distance_instance(d) {

	multiplication_factor=
	d.transformed_distance(1.0+Eps);
        
	max_k=k;
	actual_k=0;
	search_nearest = Search_nearest; 
		
        // if (search_nearest) 
		distance_to_root=
        	d.min_distance_to_rectangle(q, tree.bounding_box());
        // else 
	//	distance_to_root=
        //	distance_instance->max_distance_to_rectangle(q,
	//					tree.bounding_box());

        query_object = q;

        total_item_number=tree.size();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
       
        compute_neighbors_orthogonally(tree.root(), distance_to_root);
       
    }

    
    // Print statistics of the k_neighbor search process.
    std::ostream& statistics (std::ostream& s) {
    	s << "K_Neighbor search statistics:" << std::endl;
    	s << "Number of internal nodes visited:" 
	<< number_of_internal_nodes_visited << std::endl;
    	s << "Number of leaf nodes visited:" 
	<< number_of_leaf_nodes_visited << std::endl;
    	s << "Number of items visited:" 
	<< number_of_items_visited << std::endl;
        return s;
    }

    // destructor
    ~Orthogonal_k_neighbor_search() { 
		l.clear();  
    };

    private:
   
    void compute_neighbors_orthogonally(Node_handle N, FT rd) {
		
      typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
      typename SearchTraits::Cartesian_const_iterator_d query_object_it = construct_it(query_object);
                if (!(N->is_leaf())) {
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->cutting_dimension();
                        FT old_off, new_rd;
                        FT new_off =
                        *(query_object_it + new_cut_dim) - 
					N->cutting_value();
                        if ( ((new_off < FT(0.0)) && (search_nearest)) ||
                        (( new_off >= FT(0.0)) && (!search_nearest))  ) {
				compute_neighbors_orthogonally(N->lower(),rd);
                                if (search_nearest) {
                                	old_off= *(query_object_it + new_cut_dim)-
							N->low_value();
                                	if (old_off>FT(0.0)) old_off=FT(0.0);
                                }
				else 
				{	
                                	old_off= *(query_object_it + new_cut_dim) 
					- N->high_value();
					if (old_off<FT(0.0)) old_off=FT(0.0);
                                }
                                new_rd=
                                distance_instance.new_distance(rd,old_off,
							       new_off,
							       new_cut_dim);
				if (branch(new_rd)) 
				compute_neighbors_orthogonally(N->upper(),
								new_rd);                               
                        }
                        else { // compute new distance
                                compute_neighbors_orthogonally(N->upper(),rd); 
				if (search_nearest) {
                                	old_off= N->high_value() - 
					*(query_object_it + new_cut_dim);
                                	if (old_off>FT(0.0)) old_off=FT(0.0);
				}
                                else 
                                {       
                                	old_off= N->low_value() - 
					*(query_object_it + new_cut_dim);
					if (old_off<FT(0.0)) old_off=FT(0.0);
				}  
                                new_rd=
                                distance_instance. new_distance(rd,old_off,
								new_off,
								new_cut_dim);
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
                  for (Point_d_iterator it=N->begin(); it != N->end(); it++) {
                        number_of_items_visited++;
			FT distance_to_query_object=
                        distance_instance.transformed_distance(query_object,**it);
                        insert(*it,distance_to_query_object);
                  }
		}
    }

    
    


    
   

}; // class 



} // namespace CGAL


#endif  // ORTHOGONAL_K_NEIGHBOR_SEARCH
