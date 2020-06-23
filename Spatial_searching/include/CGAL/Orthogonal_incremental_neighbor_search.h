// Copyright (c) 2002,2011 Utrecht University (The Netherlands).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
//
// $URL$
// $Id$
// SPDX-License-Identifier: GPL-3.0-or-later OR LicenseRef-Commercial
//
//
// Author(s)     : Hans Tangelder (<hanst@cs.uu.nl>)
//                 Clement Jamin (clement.jamin.pro@gmail.com)

#ifndef CGAL_ORTHOGONAL_INCREMENTAL_NEIGHBOR_SEARCH
#define CGAL_ORTHOGONAL_INCREMENTAL_NEIGHBOR_SEARCH

#include <CGAL/license/Spatial_searching.h>

#include <CGAL/disable_warnings.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/tuple.h>
#include <CGAL/internal/Search_helpers.h>

#include <cstring>
#include <list>
#include <queue>
#include <memory>
#include <iterator> // for std::distance
#include <tuple> // std::get for tuple

namespace CGAL {

  template <class SearchTraits,
            class Distance_= typename internal::Spatial_searching_default_distance<SearchTraits>::type,
            class Splitter_ = Sliding_midpoint<SearchTraits>,
            class Tree_= Kd_tree<SearchTraits, Splitter_, Tag_true, Tag_false> >
  class Orthogonal_incremental_neighbor_search {

  public:
    typedef Splitter_ Splitter;
    typedef Tree_  Tree;
    typedef Distance_ Distance;
    typedef typename SearchTraits::Point_d Point_d;
    typedef typename Distance::Query_item Query_item;
    typedef typename SearchTraits::FT FT;
    typedef typename Tree::Point_d_iterator Point_d_iterator;
    typedef typename Tree::Node_const_handle Node_const_handle;

    typedef std::pair<Point_d,FT> Point_with_transformed_distance;
    typedef std::tuple<Node_const_handle,FT,std::vector<FT> > Node_with_distance;
    typedef std::vector<Node_with_distance*> Node_with_distance_vector;
    typedef std::vector<Point_with_transformed_distance*> Point_with_transformed_distance_vector;

    template<class T>
    struct Object_wrapper
    {
      T object;
      Object_wrapper(const T& t):object(t){}
      const T& operator* () const { return object; }
      const T* operator-> () const { return &object; }
    };

    class Iterator_implementation {
      SearchTraits traits;
    public:

      int number_of_neighbours_computed;
      int number_of_internal_nodes_visited;
      int number_of_leaf_nodes_visited;
      int number_of_items_visited;

    private:

      typedef std::vector<FT> Distance_vector;

      Distance_vector dists;

      Distance orthogonal_distance_instance;
      internal::Distance_helper<Distance, SearchTraits> m_distance_helper;

      FT multiplication_factor;

      Query_item query_point;

      FT distance_to_root;

      bool search_nearest_neighbour;

      FT rd;

      int m_dim;
      Tree const& m_tree;

      class Priority_higher {
      public:

        bool search_nearest;

        Priority_higher(bool search_the_nearest_neighbour)
          : search_nearest(search_the_nearest_neighbour)
        {}

        //highest priority is smallest distance
        bool
        operator() (Node_with_distance* n1, Node_with_distance* n2) const
        {
          return (search_nearest) ? (std::get<1>(*n1) > std::get<1>(*n2)) : (std::get<1>(*n2) > std::get<1>(*n1));
        }
      };

      class Distance_smaller {

      public:

        bool search_nearest;

        Distance_smaller(bool search_the_nearest_neighbour)
          : search_nearest(search_the_nearest_neighbour)
        {}

        //highest priority is smallest distance
        bool operator() (Point_with_transformed_distance* p1, Point_with_transformed_distance* p2) const
        {
          return (search_nearest) ? (p1->second > p2->second) : (p2->second > p1->second);
        }
      };


      std::priority_queue<Node_with_distance*, Node_with_distance_vector,
                          Priority_higher> PriorityQueue;

    public:
      std::priority_queue<Point_with_transformed_distance*, Point_with_transformed_distance_vector,
                          Distance_smaller> Item_PriorityQueue;


    public:

      int reference_count;



      // constructor
      Iterator_implementation(const Tree& tree,const Query_item& q, const Distance& tr,
                              FT Eps=FT(0.0), bool search_nearest=true)
        : traits(tree.traits()),number_of_neighbours_computed(0), number_of_internal_nodes_visited(0),
        number_of_leaf_nodes_visited(0), number_of_items_visited(0),
        orthogonal_distance_instance(tr),
        m_distance_helper(orthogonal_distance_instance, traits),
        multiplication_factor(orthogonal_distance_instance.transformed_distance(FT(1.0)+Eps)),
        query_point(q), search_nearest_neighbour(search_nearest),
        m_tree(tree),
        PriorityQueue(Priority_higher(search_nearest)), Item_PriorityQueue(Distance_smaller(search_nearest)),
        reference_count(1)


      {
        if (m_tree.empty()) return;

        typename SearchTraits::Construct_cartesian_const_iterator_d ccci=traits.construct_cartesian_const_iterator_d_object();
        m_dim = static_cast<int>(std::distance(ccci(q), ccci(q,0)));

        dists.resize(m_dim);
        for(int i=0 ; i<m_dim ; ++i){
          dists[i] = 0;
        }

        if (search_nearest){
          distance_to_root=
            orthogonal_distance_instance.min_distance_to_rectangle(q, m_tree.bounding_box(),dists);
          Node_with_distance *The_Root = new Node_with_distance(m_tree.root(),
                                                                distance_to_root, dists);
          PriorityQueue.push(The_Root);

          // rd is the distance of the top of the priority queue to q
          rd=std::get<1>(*The_Root);
          Compute_the_next_nearest_neighbour();
        }
         else{
           distance_to_root=
            orthogonal_distance_instance.max_distance_to_rectangle(q,
                                                m_tree.bounding_box(), dists);
           Node_with_distance *The_Root = new Node_with_distance(m_tree.root(),
                                                                 distance_to_root, dists);
        PriorityQueue.push(The_Root);

        // rd is the distance of the top of the priority queue to q
        rd=std::get<1>(*The_Root);
        Compute_the_next_furthest_neighbour();
         }


      }

      // * operator
      const Point_with_transformed_distance&
      operator* () const
      {
        return *(Item_PriorityQueue.top());
      }

      // prefix operator
      Iterator_implementation&
      operator++()
      {
        Delete_the_current_item_top();
        if(search_nearest_neighbour)
          Compute_the_next_nearest_neighbour();
        else
          Compute_the_next_furthest_neighbour();
        return *this;
      }

      // postfix operator
      Object_wrapper<Point_with_transformed_distance>
      operator++(int)
      {
        Object_wrapper<Point_with_transformed_distance> result( *(Item_PriorityQueue.top()) );
        ++*this;
        return result;
      }

      // Print statistics of the general priority search process.
      std::ostream&
      statistics (std::ostream& s) const {
            s << "Orthogonal priority search statistics:"
          << std::endl;
            s << "Number of internal nodes visited:"
          << number_of_internal_nodes_visited << std::endl;
            s << "Number of leaf nodes visited:"
          << number_of_leaf_nodes_visited << std::endl;
            s << "Number of items visited:"
          << number_of_items_visited << std::endl;
        s << "Number of neighbours computed:"
          << number_of_neighbours_computed << std::endl;
        return s;
      }


      //destructor
      ~Iterator_implementation()
      {
        while (!PriorityQueue.empty()) {
          Node_with_distance* The_top=PriorityQueue.top();
          PriorityQueue.pop();
          delete The_top;
        }
        while (!Item_PriorityQueue.empty()) {
          Point_with_transformed_distance* The_top=Item_PriorityQueue.top();
          Item_PriorityQueue.pop();
          delete The_top;
        }
      }

    private:

      void
      Delete_the_current_item_top()
      {
        Point_with_transformed_distance* The_item_top=Item_PriorityQueue.top();
        Item_PriorityQueue.pop();
        delete The_item_top;
      }



      // With cache
      bool search_in_leaf(typename Tree::Leaf_node_const_handle node, Tag_true, bool search_furthest)
      {
        typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();
        typename std::vector<FT>::const_iterator cache_point_begin = m_tree.cache_begin() + m_dim*(it_node_point - m_tree.begin());

        for (; it_node_point != it_node_point_end; ++it_node_point)
        {
          number_of_items_visited++;
          FT distance_to_query_point =
            m_distance_helper.transformed_distance_from_coordinates(
              query_point, *it_node_point, cache_point_begin, cache_point_begin + m_dim);

          Point_with_transformed_distance *NN_Candidate =
            new Point_with_transformed_distance(*it_node_point, distance_to_query_point);
          Item_PriorityQueue.push(NN_Candidate);

          cache_point_begin += m_dim;
        }
        // old top of PriorityQueue has been processed,
        // hence update rd

        bool next_neighbour_found;
        if (!(PriorityQueue.empty()))
        {
          rd = std::get<1>(*PriorityQueue.top());
          next_neighbour_found = (search_furthest ?
            multiplication_factor*rd < Item_PriorityQueue.top()->second
            : multiplication_factor*rd > Item_PriorityQueue.top()->second);
        }
        else // priority queue empty => last neighbour found
        {
          next_neighbour_found = true;
        }

        number_of_neighbours_computed++;
        return next_neighbour_found;
      }

      // Without cache
      bool search_in_leaf(typename Tree::Leaf_node_const_handle node, Tag_false, bool search_furthest)
      {
        typename Tree::iterator it_node_point = node->begin(), it_node_point_end = node->end();

        for (; it_node_point != it_node_point_end; ++it_node_point)
        {
          number_of_items_visited++;
          FT distance_to_query_point=
            orthogonal_distance_instance.transformed_distance(query_point, *it_node_point);

          Point_with_transformed_distance *NN_Candidate =
            new Point_with_transformed_distance(*it_node_point, distance_to_query_point);
          Item_PriorityQueue.push(NN_Candidate);
        }
        // old top of PriorityQueue has been processed,
        // hence update rd

        bool next_neighbour_found;
        if (!(PriorityQueue.empty()))
        {
          rd = std::get<1>(*PriorityQueue.top());
          next_neighbour_found = (search_furthest ?
            multiplication_factor*rd < Item_PriorityQueue.top()->second
            : multiplication_factor*rd > Item_PriorityQueue.top()->second);
        }
        else // priority queue empty => last neighbour found
        {
          next_neighbour_found=true;
        }

        number_of_neighbours_computed++;
        return next_neighbour_found;
      }


      void
      Compute_the_next_nearest_neighbour()
      {
        // compute the next item
        bool next_neighbour_found=false;
        if (!(Item_PriorityQueue.empty())) {
            next_neighbour_found=
              (multiplication_factor*rd > Item_PriorityQueue.top()->second);
        }
        typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
        typename SearchTraits::Cartesian_const_iterator_d query_point_it = construct_it(query_point);
        // otherwise browse the tree further
        while ((!next_neighbour_found) && (!PriorityQueue.empty())) {
          Node_with_distance* The_node_top=PriorityQueue.top();
          Node_const_handle N= std::get<0>(*The_node_top);
          dists = std::get<2>(*The_node_top);
          PriorityQueue.pop();
          delete The_node_top;
          FT copy_rd=rd;
          while (!(N->is_leaf())) { // compute new distance
            typename Tree::Internal_node_const_handle node =
            static_cast<typename Tree::Internal_node_const_handle>(N);
            number_of_internal_nodes_visited++;
            int new_cut_dim=node->cutting_dimension();
            FT new_rd,dst = dists[new_cut_dim];
            FT val = *(query_point_it + new_cut_dim);
            FT diff1 = val - node->upper_low_value();
            FT diff2 = val - node->lower_high_value();
            if (diff1 + diff2 < FT(0.0)) {
              new_rd=
                orthogonal_distance_instance.new_distance(copy_rd,dst,diff1,new_cut_dim);

              CGAL_assertion(new_rd >= copy_rd);
                dists[new_cut_dim] = diff1;
                Node_with_distance *Upper_Child =
                  new Node_with_distance(node->upper(), new_rd, dists);
                PriorityQueue.push(Upper_Child);
                dists[new_cut_dim] = dst;
                N=node->lower();

            }
            else { // compute new distance
              new_rd=orthogonal_distance_instance.new_distance(copy_rd,dst,diff2,new_cut_dim);
              CGAL_assertion(new_rd >= copy_rd);
                dists[new_cut_dim] = diff2;
                Node_with_distance *Lower_Child =
                  new Node_with_distance(node->lower(), new_rd, dists);
                PriorityQueue.push(Lower_Child);
                dists[new_cut_dim] = dst;
                N=node->upper();
            }
          }
          // n is a leaf
          typename Tree::Leaf_node_const_handle node =
            static_cast<typename Tree::Leaf_node_const_handle>(N);
          number_of_leaf_nodes_visited++;
    if (node->size() > 0) {
      typename internal::Has_points_cache<Tree, internal::has_Enable_points_cache<Tree>::type::value>::type dummy;
      next_neighbour_found = search_in_leaf(node, dummy, false);
    }
        }   // next_neighbour_found or priority queue is empty
        // in the latter case also the item priority quee is empty
      }


      void
      Compute_the_next_furthest_neighbour()
      {
        // compute the next item
        bool next_neighbour_found=false;
        if (!(Item_PriorityQueue.empty())) {
            next_neighbour_found=
              (rd < multiplication_factor*Item_PriorityQueue.top()->second);
        }
        typename SearchTraits::Construct_cartesian_const_iterator_d construct_it=traits.construct_cartesian_const_iterator_d_object();
        typename SearchTraits::Cartesian_const_iterator_d query_point_it = construct_it(query_point);
        // otherwise browse the tree further
        while ((!next_neighbour_found) && (!PriorityQueue.empty())) {
          Node_with_distance* The_node_top=PriorityQueue.top();
          Node_const_handle N= std::get<0>(*The_node_top);
          dists = std::get<2>(*The_node_top);
          PriorityQueue.pop();
          delete The_node_top;
          FT copy_rd=rd;
          while (!(N->is_leaf())) { // compute new distance
            typename Tree::Internal_node_const_handle node =
              static_cast<typename Tree::Internal_node_const_handle>(N);
            number_of_internal_nodes_visited++;
            int new_cut_dim=node->cutting_dimension();
            FT new_rd,dst = dists[new_cut_dim];
            FT val = *(query_point_it + new_cut_dim);
            FT diff1 = val - node->upper_low_value();
            FT diff2 = val - node->lower_high_value();
            if (diff1 + diff2 < FT(0.0)) {
              diff1 = val - node->upper_high_value();
              new_rd=
                orthogonal_distance_instance.new_distance(copy_rd,dst,diff1,new_cut_dim);
                Node_with_distance *Lower_Child =
                  new Node_with_distance(node->lower(), copy_rd, dists);
                PriorityQueue.push(Lower_Child);
                N=node->upper();
                dists[new_cut_dim] = diff1;
                copy_rd=new_rd;

            }
            else { // compute new distance
              diff2 = val - node->lower_low_value();
              new_rd=orthogonal_distance_instance.new_distance(copy_rd,dst,diff2,new_cut_dim);
                Node_with_distance *Upper_Child =
                  new Node_with_distance(node->upper(), copy_rd, dists);
                PriorityQueue.push(Upper_Child);
                N=node->lower();
                dists[new_cut_dim] = diff2;
                copy_rd=new_rd;
            }
          }
          // n is a leaf
          typename Tree::Leaf_node_const_handle node =
            static_cast<typename Tree::Leaf_node_const_handle>(N);
          number_of_leaf_nodes_visited++;
          if (node->size() > 0) {
            typename internal::Has_points_cache<Tree, internal::has_Enable_points_cache<Tree>::type::value>::type dummy;
            next_neighbour_found = search_in_leaf(node, dummy, true);
          }
        }   // next_neighbour_found or priority queue is empty
        // in the latter case also the item priority quee is empty
      }
    }; // class Iterator_implementaion









  public:
    class iterator;
    typedef iterator const_iterator;

    // constructor
    Orthogonal_incremental_neighbor_search(const Tree& tree,
                                           const Query_item& q, FT Eps = FT(0.0),
                                           bool search_nearest=true, const Distance& tr=Distance())
      : m_tree(tree),m_query(q),m_dist(tr),m_Eps(Eps),m_search_nearest(search_nearest)
    {}

    iterator
    begin() const
    {
      return iterator(m_tree,m_query,m_dist,m_Eps,m_search_nearest);
    }

    iterator
    end() const
    {
      return iterator();
    }

    std::ostream&
    statistics(std::ostream& s)
    {
      begin()->statistics(s);
      return s;
    }




    class iterator {

    public:

      typedef std::input_iterator_tag iterator_category;
      typedef Point_with_transformed_distance       value_type;
      typedef const Point_with_transformed_distance*      pointer;
      typedef const Point_with_transformed_distance&      reference;
      typedef std::size_t               size_type;
      typedef std::ptrdiff_t            difference_type;
      typedef int distance_type;

      //class Iterator_implementation;
      Iterator_implementation *Ptr_implementation;


    public:

      // default constructor
      iterator()
      : Ptr_implementation(0)
      {}

      int
      the_number_of_items_visited()
      {
        return Ptr_implementation->number_of_items_visited;
      }

      // constructor
      iterator(const Tree& tree,const Query_item& q, const Distance& tr=Distance(), FT eps=FT(0.0),
               bool search_nearest=true)
        : Ptr_implementation(new Iterator_implementation(tree, q, tr, eps, search_nearest))
        {}

      // copy constructor
      iterator(const iterator& Iter)
      {
        Ptr_implementation = Iter.Ptr_implementation;
        if (Ptr_implementation != 0) Ptr_implementation->reference_count++;
      }

      iterator& operator=(const iterator& Iter)
      {
        if (Ptr_implementation != Iter.Ptr_implementation){
          if (Ptr_implementation != 0 && --(Ptr_implementation->reference_count)==0) {
              delete Ptr_implementation;
          }
          Ptr_implementation = Iter.Ptr_implementation;
          if (Ptr_implementation != 0) Ptr_implementation->reference_count++;
        }
        return *this;
      }


      const Point_with_transformed_distance&
      operator* () const
      {
        return *(*Ptr_implementation);
      }

      // -> operator
      const Point_with_transformed_distance*
      operator-> () const
      {
        return &*(*Ptr_implementation);
      }

      // prefix operator
      iterator&
      operator++()
      {
        ++(*Ptr_implementation);
        return *this;
      }

      // postfix operator
      Object_wrapper<Point_with_transformed_distance>
      operator++(int)
      {
        return (*Ptr_implementation)++;
      }


      bool
      operator==(const iterator& It) const
      {
        if (
            ((Ptr_implementation == 0) ||
             Ptr_implementation->Item_PriorityQueue.empty()) &&
            ((It.Ptr_implementation == 0) ||
             It.Ptr_implementation->Item_PriorityQueue.empty())
            )
          return true;
        // else
        return (Ptr_implementation == It.Ptr_implementation);
      }

      bool
      operator!=(const iterator& It) const
      {
        return !(*this == It);
      }

      std::ostream&
      statistics (std::ostream& s)
      {
            Ptr_implementation->statistics(s);
        return s;
      }

      ~iterator()
      {
        if (Ptr_implementation != 0) {
          Ptr_implementation->reference_count--;
          if (Ptr_implementation->reference_count==0) {
            delete Ptr_implementation;
            Ptr_implementation = 0;
          }
        }
      }


    }; // class iterator

    //data members
    const Tree& m_tree;
    Query_item m_query;
    Distance m_dist;
    FT m_Eps;
    bool m_search_nearest;
  }; // class

  template <class Traits, class Query_item, class Distance>
  void swap (typename Orthogonal_incremental_neighbor_search<Traits,
             Query_item, Distance>::iterator& x,
             typename Orthogonal_incremental_neighbor_search<Traits,
             Query_item, Distance>::iterator& y)
  {
    typename Orthogonal_incremental_neighbor_search<Traits,
      Query_item, Distance>::iterator::Iterator_implementation
      *tmp = x.Ptr_implementation;
    x.Ptr_implementation  = y.Ptr_implementation;
    y.Ptr_implementation = tmp;
  }

} // namespace CGAL

#include <CGAL/enable_warnings.h>

#endif // CGAL_ORTHOGONAL_INCREMENTAL_NEIGHBOR_SEARCH_H
