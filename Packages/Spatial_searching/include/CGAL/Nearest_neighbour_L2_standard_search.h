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
// file          : include/CGAL/Nearest_neighbour_L2_standard_search.h
// package       : ASPAS
// revision      : 1.4 
// revision_date : 2002/16/08 
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================

#ifndef  NEAREST_NEIGHBOUR_L2_STANDARD_SEARCH_H
#define  NEAREST_NEIGHBOUR_L2_STANDARD_SEARCH_H
#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree_traits_point.h>
// #include <CGAL/Box.h>
namespace CGAL {

template <class NT, class Point> 
  NT Min_squared_distance_l2_to_box(const Point& p,
					      const Kd_tree_rectangle<NT>& r) {
    NT distance(0.0);
    NT h;
    for (int i = 0; i < r.dimension(); ++i) {
      h=p[i];
      if (h < r.lower(i)) distance += (r.lower(i)-h)*(r.lower(i)-h);
	  if (h > r.upper(i)) distance += (h-r.upper(i))*(h-r.upper(i));
	}
    return distance;
  }

  template <class NT, class Point> 
  NT Max_squared_distance_l2_to_box(const Point& p, const Kd_tree_rectangle<NT>& r) {
    NT distance(0.0);
    NT h;
    for (int i = 0; i < r.dimension(); ++i) {
      h=p[i];
      if (h >= (r.lower(i)+r.upper(i))/2.0) 
		  distance += (h-r.lower(i))*(h-r.lower(i)); 
	  else
		  distance += (r.upper(i)-h)*(r.upper(i)-h);
	}
    return distance;
  }

template <class Tree_traits, class Search_traits> //= Kd_tree_traits_2d>
class Nearest_neighbour_L2_standard_search {

public:

typedef typename Tree_traits::Item Item;
typedef typename Tree_traits::NT NT;
typedef Item** Item_iterator;
typedef std::pair<Item*,NT> Item_with_distance;

typedef Kd_tree_node<Tree_traits> Node;
typedef Kd_tree<Tree_traits> Tree;

private:

int number_of_neighbours_computed;

int number_of_internal_nodes_visited;
int number_of_leaf_nodes_visited;
int number_of_items_visited;

bool search_nearest;

NT multiplication_factor;
Item* query_point;
int total_item_number;
NT distance_to_root;   
int dim;

typedef std::list<Item_with_distance> NN_list;

NN_list l;
int max_k;
int actual_k;



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
	OutputIterator  the_k_nearest_neighbours(OutputIterator res)
	{   
		typename NN_list::iterator it=l.begin(); 
		for (; it != l.end(); it++) { *res= *it; res++; }
		return res;     
	};


    // constructor
    Nearest_neighbour_L2_standard_search(Tree& tree, Item& q, int k, NT Eps=0.0) {

		multiplication_factor=(1.0+Eps) * (1.0+Eps);
        
		max_k=k;
		actual_k=0;
		Search_traits s;
		search_nearest = s.Search_nearest(); 
		
        if (search_nearest) distance_to_root=
        Min_squared_distance_l2_to_box<NT,Item>(q,*(tree.bounding_box()));
        else distance_to_root=
        Max_squared_distance_l2_to_box<NT,Item>(q,*(tree.bounding_box()));
       
        query_point = &q;
        dim=query_point->dimension();

        total_item_number=tree.item_number();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
        number_of_neighbours_computed=0;

        compute_the_k_nearest_neighbours(tree.root(), distance_to_root);
    }

   
    void compute_the_k_nearest_neighbours(Node* N, NT rd) {
		// std::cout << "internal_nodes_visited= " << number_of_internal_nodes_visited << std::endl;
        // std::cout << "leaf_nodes_visited= " << number_of_leaf_nodes_visited << std::endl;
                if (!(N->is_leaf())) {
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->separator()->cutting_dimension();
                        NT old_off, new_rd;
                        NT new_off =
                        (*query_point)[new_cut_dim] - 
						N->separator()->cutting_value();
                        if ( ((new_off < 0.0) && (search_nearest)) ||
                        (( new_off >= 0.0) && (!search_nearest))  ) {
								compute_the_k_nearest_neighbours(N->lower(),rd);
                                old_off= (*query_point)[new_cut_dim]-
								N->low_value();
                                if (old_off>0.0) old_off=0.0;
                                new_rd=rd - old_off*old_off + 
								new_off*new_off;
								if (branch(new_rd)) compute_the_k_nearest_neighbours(N->upper(),new_rd);                               
                        }
                        else { // compute new distance
                                compute_the_k_nearest_neighbours(N->upper(),rd); 
                                old_off= N->high_value() - 
								(*query_point)[new_cut_dim];
                                if (old_off>0.0) old_off=0.0;
                                new_rd=rd - old_off*old_off + 
								new_off*new_off;
								if (branch(new_rd)) compute_the_k_nearest_neighbours(N->lower(),new_rd);       
                        }
                }
                else
				{
                  // n is a leaf
                  number_of_leaf_nodes_visited++;
                  for (Item_iterator it=N->begin(); it != N->end(); it++) {
                        number_of_items_visited++;
                        NT distance_to_query_point=0.0;
                        for (int i=0; i< dim; i++)       {
								NT h=((*query_point)[i]- (*(*it))[i]);
                                distance_to_query_point += h*h;
                        }
						// std::cout << "before insert" << std::endl;
                        insert(*it,distance_to_query_point);
						// std::cout << "after insert" << std::endl;
                  }
		}
}

// destructor
~Nearest_neighbour_L2_standard_search () { 
	l.clear();  
};

}; // class 



} // namespace CGAL


#endif  // NEAREST_NEIGHBOUR_L2_STANDARD_SEARCH_H
