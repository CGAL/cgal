//
//=======================================================================
// Copyright 1997, 1998, 1999, 2000 University of Notre Dame.
// Authors: Andrew Lumsdaine, Lie-Quan Lee, Jeremy G. Siek
//
// Distributed under the Boost Software License, Version 1.0. (See
// accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//=======================================================================
//
// $URL$
// $Id$
//
// NOTE: this file have been taken from boost 1.46.1 for using
//       with Modificable_priority_queue (to enhance the 
//       non-documented mutable_queue).
//       original file is <boost/pending/mutable_queue.hpp>
//
#ifndef CGAL_INTERNAL_BOOST_MUTABLE_QUEUE_HPP
#define CGAL_INTERNAL_BOOST_MUTABLE_QUEUE_HPP

#include <vector>
#include <algorithm>
#include <functional>
#include <CGAL/internal/boost/property_map.h>
#include <CGAL/internal/boost/mutable_heap.hpp>
#include <CGAL/internal/boost/array_binary_tree.hpp>
#include <iterator>

namespace CGAL { namespace internal {
namespace boost_ {

  // The mutable queue whose elements are indexed
  //
  // This adaptor provides a special kind of priority queue that has
  // and update operation. This allows the ordering of the items to
  // change. After the ordering criteria for item x changes, one must
  // call the Q.update(x)
  //
  // In order to efficiently find x in the queue, a functor must be
  // provided to map value_type to a unique ID, which the
  // mutable_queue will then use to map to the location of the
  // item. The ID's generated must be between 0 and N, where N is the
  // value passed to the constructor of mutable_queue

  template <class IndexedType,
            class RandomAccessContainer = std::vector<IndexedType>,
            class Comp = std::less<typename RandomAccessContainer::value_type>,
            class ID = ::boost::identity_property_map >
  class mutable_queue {
  public:
    typedef IndexedType value_type;
    typedef typename RandomAccessContainer::size_type size_type;
  protected:
    typedef typename RandomAccessContainer::iterator iterator;
#if !defined BOOST_NO_STD_ITERATOR_TRAITS
    typedef array_binary_tree_node<iterator, ID> Node;
#else
    typedef array_binary_tree_node<iterator, value_type, ID> Node;
#endif
    typedef compare_array_node<RandomAccessContainer,Comp> Compare;
    typedef std::vector<size_type> IndexArray;
  public:
    typedef Compare value_compare;
    typedef ID id_generator;

    mutable_queue(size_type n, const Comp& x, const ID& _id)
      : index_array(n), comp(x), id(_id) {
      c.reserve(n);
    }
//SL: added this constructor so that index_array is filled with
//    indices equals to n. Maintaining this property in pop allows
//    to have a method to detect if an element is in the queue
    mutable_queue(size_type n, const Comp& x, const ID& _id,bool)
      : index_array(n,n), comp(x), id(_id) {
      c.reserve(n);
    }
    
    template <class ForwardIterator>
    mutable_queue(ForwardIterator first, ForwardIterator last,
                  const Comp& x, const ID& _id)
      : index_array(std::distance(first, last)), comp(x), id(_id)
    {
      while( first != last ) {
        push(*first);
        ++first;
      }
    }

    bool empty() const { return c.empty(); }
//SL: modified this function so that the element popped from the queue
//    has its index set to the max number of element. That way we know
//    that an element is not in the tree if its index is the max number
//    of elements.
    void pop() {
      value_type tmp = c.back();
      c.back() = c.front();
      c.front() = tmp;

      size_type id_f = get(id, c.back());
      size_type id_b = get(id, tmp);
      //SL was: size_type i = index_array[ id_b ];
      index_array[ id_b ] = index_array[ id_f ];
      //SL was: index_array[ id_f ] = i; 
      index_array[ id_f ] = index_array.size(); /*SL added*/
      c.pop_back();
      Node node(c.begin(), c.end(), c.begin(), id);
      down_heap(node, comp, index_array);
    }
    void push(const IndexedType& x) {
      c.push_back(x);
      /*set index-array*/
      index_array[ get(id, x) ] = c.size()-1;
      Node node(c.begin(), c.end(), c.end() - 1, id);
      up_heap(node, comp, index_array);
    }

    void update(const IndexedType& x) {
      size_type current_pos = index_array[ get(id, x) ];
      c[current_pos] = x;

      Node node(c.begin(), c.end(), c.begin()+current_pos, id);
      update_heap(node, comp, index_array);
    }

    value_type& front() { return c.front(); }
    value_type& top() { return c.front(); }

    const value_type& front() const { return c.front(); }
    const value_type& top() const { return c.front(); }

    size_type size() const { return c.size(); }

    void clear() { c.clear(); }

#if 0
        // dwa 2003/7/11 - I don't know what compiler is supposed to
        // be able to compile this, but is_heap is not standard!!
    bool test() {
      return std::is_heap(c.begin(), c.end(), Comp());
    }
#endif

   protected:
    IndexArray index_array;
    Compare comp;
    RandomAccessContainer c;
    ID id;
  };


} } } //namespace CGAL::internal::boost_
#endif // CGAL_INTERNAL_BOOST_MUTABLE_QUEUE_HPP
