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
// file          : include/CGAL/Nearest_neighbour_Linf.h
// package       : APSPAS
// revision      : 1.0 
// revision_date : 2001/06/12 
// maintainer    : Hans Tangelder (<hanst@cs.uu.nl>)
//
// ======================================================================

#ifndef  NEAREST_NEIGHBOUR_L2_INF_H
#define  NEAREST_NEIGHBOUR_L2_INF_H

#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Extended_internal_node.h>
#include <CGAL/Kd_tree_traits_point.h>
#include <CGAL/Box.h>

using std::list; // to avoid compiler crash on MSVC++

namespace CGAL {

template <class Tree_traits, class Search_traits> //= Kd_tree_traits_2d>
class Nearest_neighbour_Linf {

    private:

    typedef Tree_traits::Item Item;
    typedef Item::FT NT;
    typedef Item** Item_iterator;
    typedef Base_node<Tree_traits> Node;
    typedef Binary_search_tree<Tree_traits> Tree;

    typedef Tree_traits::Item_with_distance Item_with_distance;
    typedef std::pair<Node*,NT> Node_with_distance;

    public:

    class iterator;

    private:

    class Priority_higher
    {
    public:

        bool search_nearest;

        Priority_higher() {
                Search_traits s;
                search_nearest = s.Search_nearest();
        }
        //highest priority is smallest distance
        bool operator() (Node_with_distance* n1, Node_with_distance* n2) const {
                if (search_nearest) { return (n1->second > n2->second);}
                else {return (n2->second > n1->second);}
        }
    };

    class Distance_smaller
    {

    public:

        bool search_nearest;

        Distance_smaller() {
                Search_traits s;
                search_nearest = s.Search_nearest();
        }

        //highest priority is smallest distance
        bool operator() (Item_with_distance* p1, Item_with_distance* p2) const {
		if (search_nearest) {return (p1->second > p2->second);}
                else {return (p2->second > p1->second);}
        }
    };

    typedef std::vector<Node_with_distance*> Node_with_distance_vector;

    typedef std::vector<Item_with_distance*> Item_with_distance_vector;

    typedef std::vector<NT> Distance_vector;

    iterator *start;
    iterator *past_the_end;

    public:

    // constructor
    Nearest_neighbour_Linf(Tree& tree, Item& q, NT Eps=0.0) {
		start = new iterator(tree,q,Eps);
                past_the_end = new iterator();
	};

    // destructor
    ~Nearest_neighbour_Linf() {
		delete start;
                delete past_the_end;
    };

    iterator begin() {
		return *start;
    }

    iterator end() {
		return *past_the_end;
    }

    class iterator {

    public:

    typedef typename std::input_iterator_tag iterator_category;
    typedef typename Item_with_distance value_type;
    typedef int distance_type;

    class Iterator_implementation;
    Iterator_implementation *Ptr_implementation;

    public:

    // default constructor
    class iterator() {Ptr_implementation=0;}

    int the_number_of_items_visited() {
        return Ptr_implementation->number_of_items_visited;
    }

    // constructor
    class iterator(Tree& tree, Item& q, NT eps=0.0) {
        Ptr_implementation = new Iterator_implementation(tree, q, eps);
    }

    // copy constructor
    class iterator(iterator& Iter) {
        Ptr_implementation = Iter.Ptr_implementation;
        if (Ptr_implementation != 0) Ptr_implementation->reference_count++;

    }

    Item_with_distance& operator* () {
                return *(*Ptr_implementation);
    }

    // prefix operator
    iterator& operator++() {
        ++(*Ptr_implementation);
        return *this;
    }

    // postfix operator
    std::auto_ptr<Item_with_distance> operator++(int) {
        std::auto_ptr<Item_with_distance> result = (*Ptr_implementation)++;
        return result;
    }

    bool operator==(const iterator& It) const {

        if (
                ((Ptr_implementation == 0) || Ptr_implementation->Item_PriorityQueue.empty()) &&
                ((It.Ptr_implementation == 0) ||  It.Ptr_implementation->Item_PriorityQueue.empty())
        )
        return true;
        // else
        return (Ptr_implementation == It.Ptr_implementation);
    }

    bool operator!=(const iterator& It) const {
        return !(*this == It);
    }

    ~iterator() {
        // std::cout << "called ~iterator"  << std::endl;
        if (Ptr_implementation != 0) {
               //  std::cout << "reference_count is" << Ptr_implementation->reference_count << std::endl;
                Ptr_implementation->reference_count--;
                if (Ptr_implementation->reference_count==0) {
                        delete Ptr_implementation;
                        Ptr_implementation = 0;
                }
               // if (Ptr_implementation != 0) std::cout << "new reference_count is" << Ptr_implementation->reference_count << std::endl;
        }
        // else std::cout << "Ptr_implementation is null" << std::endl;
    }



    class Iterator_implementation {

    public:

    int number_of_neighbours_computed;
    int number_of_internal_nodes_visited;
    int number_of_leaf_nodes_visited;
    int number_of_items_visited;

    private:

    bool search_nearest;

    NT multiplication_factor;

    Item* query_point;

    int total_item_number;

    NT distance_to_root;

    NT rd;

    std::priority_queue<Node_with_distance*, Node_with_distance_vector,
    Priority_higher> PriorityQueue;

	public:

	int reference_count;

    std::priority_queue<Item_with_distance*, Item_with_distance_vector,
    Distance_smaller> Item_PriorityQueue;

    // constructor
    Iterator_implementation(Tree& tree, Item& q, NT Eps=0.0) {

        reference_count=1;

        multiplication_factor=(1.0+Eps);

        Search_traits s;
        search_nearest = s.Search_nearest();

        if (search_nearest) distance_to_root=
        Min_distance_linf_to_box<NT,Item>(q,*(tree.bounding_box()));
        else distance_to_root=
        Max_distance_linf_to_box<NT,Item>(q,*(tree.bounding_box()));
       
       
        query_point = &q;

        total_item_number=tree.item_number();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
        number_of_neighbours_computed=0;



        Node_with_distance *The_Root = new Node_with_distance(tree.root(),distance_to_root);
        PriorityQueue.push(The_Root);

        // rd is the distance of the top of the priority queue to q
        rd=The_Root->second;
        Compute_the_next_nearest_neighbour();
    }

    // * operator
    Item_with_distance& operator* () {
			return *(Item_PriorityQueue.top());
    }

    // prefix operator
    Iterator_implementation& operator++() {
        // std::cout << "called ++" << std::endl;
        Delete_the_current_item_top();
        Compute_the_next_nearest_neighbour();
        return *this;
    }

    // postfix operator
    std::auto_ptr<Item_with_distance> operator++(int) {
        Item_with_distance Value = *(Item_PriorityQueue.top());
        std::auto_ptr<Item_with_distance>
        result(new Item_with_distance(Value));
        ++*this;
        return result;
    }


    //destructor
    ~Iterator_implementation() {

        // std::cout << "called destructor" << std::endl;
	while (PriorityQueue.size()>0) {
                Node_with_distance* The_top=PriorityQueue.top();
                PriorityQueue.pop();
                delete The_top;
	};
	while (Item_PriorityQueue.size()>0) {
                Item_with_distance* The_top=Item_PriorityQueue.top();
                Item_PriorityQueue.pop();
                delete The_top;
        };
    }

    private:

    void Delete_the_current_item_top() {
        Item_with_distance* The_item_top=Item_PriorityQueue.top();
        Item_PriorityQueue.pop();
        delete The_item_top;
    }

    void Compute_the_next_nearest_neighbour() {

        // compute the next item
        bool next_neighbour_found=false;

        if (!(Item_PriorityQueue.empty())) {
        if (search_nearest)
        	next_neighbour_found=(multiplication_factor*rd > Item_PriorityQueue.top()->second);
        else
		next_neighbour_found=(multiplication_factor*rd < Item_PriorityQueue.top()->second);
        }
        

       
        // otherwise browse the tree further
        while ((!next_neighbour_found) && (!PriorityQueue.empty())) {
                Node_with_distance* The_node_top=PriorityQueue.top();
                Node* N= The_node_top->first;
                PriorityQueue.pop();
                delete The_node_top;

                while (!(N->is_leaf())) { // compute new distance
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->separator()->cutting_dimension();
                        NT new_rd;
                        NT new_off =
                        (*query_point)[new_cut_dim] - N->separator()->cutting_value();
                        if ( ((new_off < 0.0) && (search_nearest)) ||
                        (( new_off >= 0.0) && (!search_nearest))  ) {
				if (-new_off > rd) new_rd=-new_off;
				else new_rd=rd;
                                Node_with_distance *Upper_Child =
                                new Node_with_distance(N->upper(),new_rd);
                                PriorityQueue.push(Upper_Child);
                                N=N->lower();
                        }
                        else { // compute new distance
				if (new_off > rd) new_rd=new_off;
                                else new_rd=rd;
                                Node_with_distance *Lower_Child =
                                new Node_with_distance(N->lower(),new_rd);
                                PriorityQueue.push(Lower_Child);
                                N=N->upper();
                        }
                }
                // n is a leaf
                number_of_leaf_nodes_visited++;
                for (Item_iterator it=N->begin(); it != N->end(); it++) {
                        number_of_items_visited++;
                        NT distance_to_query_point=0.0;
                        for (int i=0; i< (*it)->dimension(); i++)       {
                                if (fabs((*query_point)[i] - (*(*it))[i]) > distance_to_query_point)
                                distance_to_query_point = fabs((*query_point)[i] - (*(*it))[i]);
                        }
                        Item_with_distance *NN_Candidate=
                        new Item_with_distance(*it,distance_to_query_point);
                        Item_PriorityQueue.push(NN_Candidate);
                }
                // old top of PriorityQueue has been processed,
                // hence update rd
                if (!PriorityQueue.empty()) {
                        rd = PriorityQueue.top()->second;
			if (!PriorityQueue.empty()) {
                        rd = PriorityQueue.top()->second;
			if (search_nearest)
				next_neighbour_found =
                		(multiplication_factor*rd > Item_PriorityQueue.top()->second);
                        else
				next_neighbour_found =
                		(multiplication_factor*rd < Item_PriorityQueue.top()->second);
                        
                	}
                }
                else // priority queue empty => last neighbour found
                {
                        next_neighbour_found=true;
                };
                number_of_neighbours_computed++;
        }   // next_neighbour_found or priority queue is empty
        // in the latter case also the item priority quee is empty
    }
}; // class Iterator_implementaion
}; // class iterator
}; // class Nearest neighbour_L2

template <class Tree_traits, class Search_traits>
void swap (Nearest_neighbour_Linf<Tree_traits,Search_traits>::iterator& x,
        Nearest_neighbour_Linf<Tree_traits,Search_traits>::iterator& y) {
        Nearest_neighbour_Linf<Tree_traits,Search_traits>::iterator::Iterator_implementation
        *tmp = x.Ptr_implementation;
        x.Ptr_implementation  = y.Ptr_implementation;
        y.Ptr_implementation = tmp;
    }
} // namespace CGAL


#endif  // NEAREST_NEIGHBOUR_LINF_H
