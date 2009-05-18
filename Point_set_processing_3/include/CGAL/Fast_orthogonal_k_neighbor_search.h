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
// $URL:$
// $Id: $
//
//
// Author(s)     : Gael Guennebaud (gael.guennebaud@inria.fr), Hans Tangelder (<hanst@cs.uu.nl>)

#ifndef CGAL_FAST_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
#define CGAL_FAST_ORTHOGONAL_K_NEIGHBOR_SEARCH_H

#include <cstring>
#include <list>
#include <queue>
#include <set>
#include <memory>
#include <CGAL/Kd_tree_node.h>
#include <CGAL/Kd_tree.h>
#include <CGAL/Euclidean_distance.h>
#include <CGAL/Splitters.h>

namespace CGAL {


template <typename T, typename Compare = std::less<T> >
class bounded_priority_queue
{
public:

  typedef T value_type;
  typedef typename std::vector<value_type>::const_iterator const_iterator;

  bounded_priority_queue(const Compare& comp = Compare())
    : m_comp(comp)
  {}

  bounded_priority_queue(int size, const Compare& comp = Compare())
    : m_count(0), m_data(size), m_comp(comp)
  {}

  /** Sets the max number of elements in the queue */
  void resize(int new_size)
  {
    if (m_data.size()!=new_size)
      m_data.resize(new_size);
    clear();
  }

  /** \returns the number of elements in the queue */
  inline unsigned int size() const { return m_count; }

  /** Removes all elements of the queue. The max size remains unchanged. */
  inline void clear() { m_count = 0; }

  inline bool full() const { return m_count == m_data.size(); }
  inline bool empty() const { return m_count == 0; }

  /** \returns greatest element */
  inline const value_type& top() const { return m_data[0]; }

  inline void insert(const value_type& x)
  {
    value_type* __restrict data1 = (&m_data[0]-1);
    if (full())
    {
      if (m_comp(x, top()))
      {
        unsigned int j(1), k(2);
        while (k <= m_count)
        {
          value_type* z = &(data1[k]);
          if ((k < m_count) && m_comp(*z, data1[k+1]))
            z = &(data1[++k]);

          if (m_comp(*z, x))
            break;
          data1[j] = *z;
          j = k;
          k = j << 1;
        }
        data1[j] = x;
      }
    }
    else
    {
      int i(++m_count), j;
      while (i >= 2)
      {
        j = i >> 1;
        value_type& y = data1[j];
        if (m_comp(x, y))
          break;
        data1[i] = y;
        i = j;
      }
      data1[i] = x;
    }
  }

  const_iterator begin() const { return m_data.begin(); }
  const_iterator end() const
  {
    const_iterator res = m_data.begin();
    res += m_count;
    return res;
  }

  void sort()
  {
    std::sort(m_data.begin(), m_data.end(), m_comp);
  }

protected:

  unsigned int m_count;
  std::vector<value_type> m_data;
  Compare m_comp;
};


template <class SearchTraits,
          class Distance_= Euclidean_distance<SearchTraits>,
          class Splitter_= Sliding_midpoint<SearchTraits> ,
          class Tree_= Kd_tree<SearchTraits, Splitter_, Tag_true> >
class Fast_orthogonal_k_neighbor_search {

public:

  typedef Splitter_ Splitter;
  typedef Tree_  Tree;
  typedef Distance_ Distance;
  typedef typename SearchTraits::Point_d Point_d;
  typedef typename Distance::Query_item Query_item;

  typedef typename SearchTraits::FT FT;
  typedef std::pair<Point_d*,FT> Point_ptr_with_transformed_distance;
  typedef typename Tree::Node_handle Node_handle;
  typedef typename Tree::Point_d_iterator Point_d_iterator;

private:

  // Comparison functor of the second member of a std::pair
  struct Compare_pair_second
  {
    template<typename T>
    bool operator()(const T& p1, const T& p2) const
    {
      return (p1.second < p2.second);
    }
  };

  typedef bounded_priority_queue<Point_ptr_with_transformed_distance, Compare_pair_second> PQueue;

public:

  typedef typename PQueue::const_iterator iterator;

private:

  // stats:
  int number_of_internal_nodes_visited;
  int number_of_leaf_nodes_visited;
  int number_of_items_visited;

  FT multiplication_factor;
  Query_item query_object;
  int total_item_number;
  FT distance_to_root;

  PQueue pqueue;
  Distance distance_instance;

private:

  // Test if we should continue searching
  inline bool branch(FT distance)
  {
    if (!pqueue.full())
      return true;
    else
      return distance*multiplication_factor < pqueue.top().second;
  }

public:

  iterator begin() const
  {
    return pqueue.begin();
  }

  iterator end() const
  {
    return pqueue.end();
  }


  // constructor
  Fast_orthogonal_k_neighbor_search(Tree& tree, const Query_item& q,
    unsigned int k=1, FT Eps=FT(0.0), bool sorted=false, const Distance& d=Distance())
    : number_of_internal_nodes_visited(0), number_of_leaf_nodes_visited(0), number_of_items_visited(0),
    multiplication_factor(d.transformed_distance(1.0+Eps)), query_object(q),
    total_item_number(tree.size()), pqueue(k), distance_instance(d)

  {
    distance_to_root = d.min_distance_to_rectangle(q, tree.bounding_box());
    compute_neighbors_orthogonally(tree.root(), distance_to_root);
    if (sorted)
      pqueue.sort();
  }

  // Prints statistics of the k_neighbor search process.
  std::ostream& statistics (std::ostream& s)
  {
    s << "K_Neighbor search statistics:" << std::endl;
    s << "Number of internal nodes visited:"
      << number_of_internal_nodes_visited << std::endl;
    s << "Number of leaf nodes visited:"
      << number_of_leaf_nodes_visited << std::endl;
    s << "Number of items visited:"
      << number_of_items_visited << std::endl;
    return s;
  }

private:

  void compute_neighbors_orthogonally(Node_handle N, FT rd)
  {
    typename SearchTraits::Construct_cartesian_const_iterator_d construct_it;
    typename SearchTraits::Cartesian_const_iterator_d query_object_it = construct_it(query_object);
    if (!(N->is_leaf()))
    {
      number_of_internal_nodes_visited++;
      int new_cut_dim=N->cutting_dimension();
      FT old_off, new_rd;
      FT new_off = *(query_object_it + new_cut_dim) - N->cutting_value();
      if (new_off <  FT(0.0))
      {
        compute_neighbors_orthogonally(N->lower(),rd);
        old_off= *(query_object_it + new_cut_dim) - N->low_value();
        if (old_off>FT(0.0))
          old_off=FT(0.0);
        new_rd = distance_instance.new_distance(rd,old_off,new_off,new_cut_dim);
        if (branch(new_rd))
          compute_neighbors_orthogonally(N->upper(), new_rd);
      }
      else // compute new distance
      {
        compute_neighbors_orthogonally(N->upper(),rd);
        old_off= N->high_value() - *(query_object_it + new_cut_dim);
        if (old_off>FT(0.0))
          old_off=FT(0.0);
        new_rd = distance_instance. new_distance(rd,old_off,new_off,new_cut_dim);
        if (branch(new_rd))
          compute_neighbors_orthogonally(N->lower(), new_rd);
      }

    }
    else
    {
      // n is a leaf
      number_of_leaf_nodes_visited++;
      int nb = N->size();
      if (nb > 0)
      {
        number_of_items_visited+=nb;
        Point_d_iterator it = N->begin();
        for (int i=0;i<nb;++i,++it)
        {
          FT distance_to_query_object =
            distance_instance.transformed_distance(query_object,**it);
          pqueue.insert(Point_ptr_with_transformed_distance(*it,distance_to_query_object) );
        }
      }
    }
  }

}; // class

} // namespace CGAL

#endif  // CGAL_FAST_ORTHOGONAL_K_NEIGHBOR_SEARCH_H
