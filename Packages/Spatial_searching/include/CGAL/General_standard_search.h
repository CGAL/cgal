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
// file          : include/CGAL/General_standard_search.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef  CGAL_GENERAL_STANDARD_SEARCH_H
#define  CGAL_GENERAL_STANDARD_SEARCH_H
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
	  class Tree_=Kd_tree<GeomTraits, Splitter_, Tag_false> >

class General_standard_search {

public:
  typedef Splitter_ Splitter;
  typedef Distance_ Distance;
  typedef Tree_ Tree;
  typedef typename GeomTraits::Point Point;
  typedef typename GeomTraits::NT NT;
  typedef std::pair<Point,NT> Point_with_distance;
  
  typedef typename Tree::Node_handle Node_handle;
  
  typedef typename Tree::Point_iterator Point_iterator;
  typedef Kd_tree_rectangle<GeomTraits> Rectangle;
  typedef typename Distance::Query_item Query_item;
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
      else return ( distance >
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
    General_standard_search(Tree& tree, const Query_item& q, 
			    int k=1, NT Eps=NT(0.0), 
			    bool Search_nearest=true,
			    const Distance& d=Distance()) {

	distance_instance=new Distance(d);

	multiplication_factor=
	distance_instance->transformed_distance(NT(1.0)+Eps);
        
	max_k=k;
	actual_k=0;
	search_nearest = Search_nearest; 
       
        query_object = q;

        total_item_number=tree.size();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
       
        
        compute_neighbors_general(tree.root(), tree.bounding_box());

    }

    // Print statistics of the general standard search process.
    std::ostream& statistics (std::ostream& s) {
    	s << "General search statistics:" << std::endl;
    	s << "Number of internal nodes visited:" 
	<< number_of_internal_nodes_visited << std::endl;
    	s << "Number of leaf nodes visited:" 
	<< number_of_leaf_nodes_visited << std::endl;
    	s << "Number of items visited:" 
	<< number_of_items_visited << std::endl;
        return s;
    }


    // destructor
    ~General_standard_search() { 
		l.clear();  
		delete distance_instance;
   };

    private:
   

  void compute_neighbors_general(Node_handle N, Kd_tree_rectangle<GeomTraits>* r) {
		
                if (!(N->is_leaf())) {
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->cutting_dimension();
			NT  new_cut_val=N->cutting_value();

			Kd_tree_rectangle<GeomTraits>* r_lower =
			  new Kd_tree_rectangle<GeomTraits>(*r);

			// modifies also r_lower to lower half
			Kd_tree_rectangle<GeomTraits>* r_upper =
			r_lower->split(new_cut_dim, new_cut_val);

                        NT distance_to_lower_half;
                        NT distance_to_upper_half;

                        if (search_nearest) { 

                        	distance_to_lower_half = 
                        	distance_instance -> 
				min_distance_to_queryitem(query_object, 
							  *r_lower);
				
                        	distance_to_upper_half = 
                        	distance_instance -> 
				min_distance_to_queryitem(query_object, 
							  *r_upper);
			

			} 
			else
			{ 

                        	distance_to_lower_half = 
                        	distance_instance -> 
				max_distance_to_queryitem(query_object, 
							  *r_lower);

                        	distance_to_upper_half = 
                        	distance_instance -> 
				max_distance_to_queryitem(query_object, 
							  *r_upper);

			}
 
			if ( (( search_nearest) && 
			     (distance_to_lower_half < distance_to_upper_half)) 
			     ||
			     ((!search_nearest) && 
			     (distance_to_lower_half >= 
			      distance_to_upper_half))  )
			{
			   compute_neighbors_general(N->lower(), r_lower);
			   if (branch(distance_to_upper_half)) 
			   compute_neighbors_general (N->upper(), r_upper);
			}  
			else
			{	compute_neighbors_general(N->upper(), r_upper);
				if (branch(distance_to_lower_half)) 
				compute_neighbors_general (N->lower(), 
							    r_lower);
			}

			delete r_lower; delete r_upper; 
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


#endif  // STANDARD_SEARCH
