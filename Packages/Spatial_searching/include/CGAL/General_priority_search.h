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
// file          : include/CGAL/General_priority_search.h
// package       : ASPAS (3.12)
// maintainer    : Hans Tangelder <hanst@cs.uu.nl>
// revision      : 3.0
// revision_date : 2003/07/10
// authors       : Hans Tangelder (<hanst@cs.uu.nl>)
// coordinator   : Utrecht University
//
// ======================================================================


#ifndef  GENERAL_PRIORITY_SEARCH_H
#define  GENERAL_PRIORITY_SEARCH_H
#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree_rectangle.h>
#include <CGAL/Euclidean_distance.h>
namespace CGAL {

template <class TreeTraits, 
          class Distance=Euclidean_distance<typename TreeTraits::Point>, 
	  class QueryItem=typename TreeTraits::Point, 
	  class Tree=Kd_tree<TreeTraits> >
class General_priority_search { 

public:

typedef typename TreeTraits::Point Point;
typedef typename TreeTraits::NT NT;
typedef typename Tree::Point_iterator Point_iterator;
typedef typename Tree::Node_handle Node_handle;

typedef Kd_tree_rectangle<NT> Node_box;   

class Cell 
    {
    private:

    	Node_box* the_box;
    	Node_handle     the_node;

    public:

        // constructor
        Cell (Node_box* Nb, Node_handle N)
	{
		the_box = Nb;
		the_node = N;
	}

        Node_box* box() {return the_box;};
        Node_handle    node() {return the_node;};
        

	~Cell() {}

    };

    

typedef std::pair<Point*,NT> Point_with_distance;
typedef std::pair<Cell*,NT> Cell_with_distance;

// this forward declaration may problems for g++ 
class iterator;


    typedef std::vector<Cell_with_distance*> Cell_with_distance_vector;

    typedef std::vector<Point_with_distance*> Point_with_distance_vector;

    typedef std::vector<NT> Distance_vector;

    iterator *start;
    iterator *past_the_end;

    public:

    // constructor
    General_priority_search(Tree& tree, QueryItem& q, const Distance& tr=Distance(),
    NT Eps=NT(0.0), bool search_nearest=true)
    {
        start = new iterator(tree,q,tr,Eps,search_nearest);
        past_the_end = new iterator();
    };

    // destructor
    ~General_priority_search() {
	delete start;
        delete past_the_end;
    };

    iterator begin() {
		return *start;
    }

    iterator end() {
		return *past_the_end;
    }

    std::ostream&  statistics(std::ostream& s) {
	start->statistics(s);
	return s;
    }

    class iterator {

    public:

    typedef std::input_iterator_tag iterator_category;
    typedef Point_with_distance value_type;
    typedef int distance_type;

    class Iterator_implementation;
    Iterator_implementation *Ptr_implementation;

    public:

    // default constructor
    iterator() {Ptr_implementation=0;}

    int the_number_of_items_visited() {
        return Ptr_implementation->number_of_items_visited;
    }

    // constructor
    iterator(Tree& tree, QueryItem& q, const Distance& tr, NT eps, 
	     bool search_nearest) {
        Ptr_implementation =
        new Iterator_implementation(tree, q, tr, eps, search_nearest);
    }

    // copy constructor
    iterator(const iterator& Iter) {
        Ptr_implementation = Iter.Ptr_implementation;
        if (Ptr_implementation != 0) Ptr_implementation->reference_count++;
    }

    Point_with_distance& operator* () {
                return *(*Ptr_implementation);
    }

    // prefix operator
    iterator& operator++() {
        ++(*Ptr_implementation);
        return *this;
    }

    // postfix operator
    std::auto_ptr<Point_with_distance> operator++(int) {
        std::auto_ptr<Point_with_distance> result = (*Ptr_implementation)++;
        return result;
    }

    bool operator==(const iterator& It) const {

        if (
                ((Ptr_implementation == 0) || 
		Ptr_implementation->Item_PriorityQueue->empty()) &&
                ((It.Ptr_implementation == 0) ||  
		It.Ptr_implementation->Item_PriorityQueue->empty())
        )
        return true;
        // else
        return (Ptr_implementation == It.Ptr_implementation);
    }

    bool operator!=(const iterator& It) const {
        return !(*this == It);
    }

    std::ostream& statistics (std::ostream& s) {
    	Ptr_implementation->statistics(s);
        return s;
    }

    ~iterator() {
        if (Ptr_implementation != 0) {
                Ptr_implementation->reference_count--;
                if (Ptr_implementation->reference_count==0) {
                        delete Ptr_implementation;
                        Ptr_implementation = 0;
                }
        }
    }


    class Iterator_implementation {

    



    private:

    NT multiplication_factor;

    QueryItem* query_point;

    int total_item_number;

    NT distance_to_root;

    bool search_nearest_neighbour;

    NT rd;

    class Priority_higher
    {
    public:

        bool search_nearest;

        Priority_higher(bool search_the_nearest_neighbour) {
                search_nearest = search_the_nearest_neighbour;
        }

        //highest priority is smallest distance
        bool operator() (Cell_with_distance* n1, Cell_with_distance* n2) const 
        {
                if (search_nearest) { return (n1->second > n2->second);}
                else {return (n2->second > n1->second);}
        }
    };

class Distance_smaller
    {

    public:

        bool search_nearest;

        Distance_smaller(bool search_the_nearest_neighbour) {
                // Search_traits s;
                search_nearest = search_the_nearest_neighbour;
        }

        //highest priority is smallest distance
        bool operator() (Point_with_distance* p1, Point_with_distance* p2) const 
        {
		if (search_nearest) {return (p1->second > p2->second);}
                else {return (p2->second > p1->second);}
        }
    };

    std::priority_queue<Cell_with_distance*, Cell_with_distance_vector,
    Priority_higher>* PriorityQueue;

    std::priority_queue<Point_with_distance*, Point_with_distance_vector,
    Distance_smaller>* Item_PriorityQueue;
    

    Distance* Distance_instance;

    public:

    int reference_count;

    int number_of_internal_nodes_visited;
    int number_of_leaf_nodes_visited;
    int number_of_items_visited;
    int number_of_neighbours_computed;

    // constructor
    Iterator_implementation(Tree& tree, QueryItem& q,const Distance& tr,
        NT Eps, bool search_nearest)
    {
        
	
        PriorityQueue = new std::priority_queue<Cell_with_distance*, 
	Cell_with_distance_vector, Priority_higher> 
        (Priority_higher(search_nearest));

        Item_PriorityQueue= new std::priority_queue<Point_with_distance*, 
	Point_with_distance_vector,
    	Distance_smaller>
	(Distance_smaller(search_nearest));

	reference_count=1;
        Distance_instance=new Distance(tr);
        multiplication_factor=
	Distance_instance->transformed_distance(NT(1)+Eps);

        Node_box *bounding_box = new Node_box(*(tree.bounding_box()));
        
        search_nearest_neighbour=search_nearest;

        if (search_nearest) distance_to_root=
        Distance_instance->min_distance_to_queryitem(q,*bounding_box);
        else distance_to_root=
   	Distance_instance->max_distance_to_queryitem(q,*bounding_box);

        

        query_point = &q;

        total_item_number=tree.item_number();

        number_of_leaf_nodes_visited=0;
        number_of_internal_nodes_visited=0;
        number_of_items_visited=0;
        number_of_neighbours_computed=0;


        Cell *Root_Cell = new Cell(bounding_box,tree.root());
        Cell_with_distance  *The_Root = 
	new Cell_with_distance(Root_Cell,distance_to_root);

        PriorityQueue->push(The_Root);

        // rd is the distance of the top of the priority queue to q
        rd=The_Root->second;
        
        Compute_the_next_nearest_neighbour();
    }

    // * operator
    Point_with_distance& operator* () {    
			return *(Item_PriorityQueue->top());
    }

    // prefix operator
    Iterator_implementation& operator++() {
        
        Delete_the_current_item_top();
        Compute_the_next_nearest_neighbour();
        return *this;
    }

    // postfix operator
    std::auto_ptr<Point_with_distance> operator++(int) {     
        
        Point_with_distance Value = *(Item_PriorityQueue->top());        
        std::auto_ptr<Point_with_distance>
        result(new Point_with_distance(Value));
        ++*this;        
        
        return result;
    }


    // Print statistics of the general priority search process.
    std::ostream& statistics (std::ostream& s) {
    	s << "General priority search statistics:" << std::endl;
    	s << "Number of internal nodes visited:" << 
		      number_of_internal_nodes_visited << std::endl;
    	s << "Number of leaf nodes visited:" << 
	number_of_leaf_nodes_visited << std::endl;
    	s << "Number of points visited:" << 
	number_of_items_visited << std::endl;
        s << "Number of neighbours computed:" << 
	number_of_neighbours_computed << std::endl;
        return s;
    }

    //destructor
    ~Iterator_implementation() {

	while (PriorityQueue->size()>0) {
                Cell_with_distance* The_top=PriorityQueue->top();
                PriorityQueue->pop();
                delete The_top->first->box();
		delete The_top->first;
                delete The_top;
	};
	while (Item_PriorityQueue->size()>0) {
                Point_with_distance* The_top=Item_PriorityQueue->top();
                Item_PriorityQueue->pop();
                delete The_top;
        };
	delete Distance_instance;
    }

    private:

    void Delete_the_current_item_top() {
        Point_with_distance* The_item_top=Item_PriorityQueue->top();
        Item_PriorityQueue->pop();
        delete The_item_top;
    }

    void Compute_the_next_nearest_neighbour() {

      
        // compute the next item
        bool next_neighbour_found=false;
        if (!(Item_PriorityQueue->empty())) {
        if (search_nearest_neighbour)
        	next_neighbour_found=
		(multiplication_factor*rd > Item_PriorityQueue->top()->second);
        else
		next_neighbour_found=
		(rd < multiplication_factor*Item_PriorityQueue->top()->second);
        }
        // otherwise browse the tree further
        while ((!next_neighbour_found) && (!PriorityQueue->empty())) {
                
                Cell_with_distance* The_node_top=PriorityQueue->top();
                Node_handle N= The_node_top->first->node();
                Node_box* B= The_node_top->first->box();
                PriorityQueue->pop();
                delete The_node_top->first;
                delete The_node_top;

                while (!(N->is_leaf())) { // compute new distances
                        number_of_internal_nodes_visited++;
                        int new_cut_dim=N->cutting_dimension();
                        NT  new_cut_val=N->cutting_value();
                        
			Node_box* lower_box = new Node_box(*B);
                        Node_box* upper_box = 
			lower_box->split(new_cut_dim, new_cut_val);
			delete B;
			if (search_nearest_neighbour) {
NT distance_to_box_lower =
Distance_instance->min_distance_to_queryitem(*query_point, *lower_box);
NT distance_to_box_upper =
Distance_instance->min_distance_to_queryitem(*query_point, *upper_box);
if (distance_to_box_lower <= distance_to_box_upper) {
	Cell* C_upper = new Cell(upper_box, N->upper());
	Cell_with_distance *Upper_Child =
	new Cell_with_distance(C_upper,distance_to_box_upper);
	PriorityQueue->push(Upper_Child);
	N=N->lower();
	B=lower_box;
}
else {
	Cell* C_lower = new Cell(lower_box, N->lower());
	Cell_with_distance *Lower_Child =
	new Cell_with_distance(C_lower,distance_to_box_lower);
	PriorityQueue->push(Lower_Child);
	N=N->upper();
	B=upper_box;
}
                        }
			else { // search furthest
NT distance_to_box_lower =
Distance_instance->max_distance_to_queryitem(*query_point, *lower_box);
NT distance_to_box_upper =
Distance_instance->max_distance_to_queryitem(*query_point, *upper_box);
if (distance_to_box_lower >= distance_to_box_upper) {
	Cell* C_upper = new Cell(upper_box, N->upper());
	Cell_with_distance *Upper_Child =
	new Cell_with_distance(C_upper,distance_to_box_upper);
	PriorityQueue->push(Upper_Child);
	N=N->lower();
	B=lower_box;
}
else {
	Cell* C_lower = new Cell(lower_box, N->lower());
	Cell_with_distance *Lower_Child =
        new Cell_with_distance(C_lower,distance_to_box_lower);
	PriorityQueue->push(Lower_Child);
        N=N->upper();
	B=upper_box;
}
                        }
                }
                // n is a leaf
                delete B;
                number_of_leaf_nodes_visited++;
                if (N->size() > 0) {
                  for (Point_iterator it=N->begin(); it != N->end(); it++) {
                        number_of_items_visited++;
                        NT distance_to_query_point=
                        Distance_instance->
                        distance(*query_point,**it);
                        Point_with_distance *NN_Candidate=
                        new Point_with_distance(*it,distance_to_query_point);
                        Item_PriorityQueue->push(NN_Candidate);
                  }
                  // old top of PriorityQueue has been processed,
                  // hence update rd
                  if (!PriorityQueue->empty()) {
                        rd = PriorityQueue->top()->second;
                        if (search_nearest_neighbour)
				next_neighbour_found =
                  (multiplication_factor*rd > 
		   Item_PriorityQueue->top()->second);
                        else
				next_neighbour_found =
                  (multiplication_factor*rd < 
		   Item_PriorityQueue->top()->second);
                  }
                  else // priority queue empty => last neighbour found
                  {
                        next_neighbour_found=true;
                  };
                  number_of_neighbours_computed++;
          }
        }   // next_neighbour_found or priority queue is empty
        // in the latter case also the item priority queue is empty
        
    }
}; // class Iterator_implementation
}; // class iterator
}; // class 

} // namespace CGAL


#endif  // GENERAL_PRIORITY_SEARCH_H
